<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/51b52963-19de-4f67-8977-072b409dae19" />

# meta_standards_converter

Convert biological study metadata among GEO MINiML, parsed JSON, ArrayExpress MAGE-TAB, and AnnData/H5AD.

## Description

`meta_standards_converter` is a Python package and command-line toolkit for moving study metadata between GEO and ArrayExpress-compatible representations and for attaching that metadata to expression data. It can fetch and parse GEO MINiML, enrich packages with PubMed and SRA/ENA records, read and write MAGE-TAB IDF/SDRF files, normalize processed matrices into H5AD, and process raw FASTQs through pinned nf-core pipelines.

The five primary workflows are:

- `geo2ae`: GEO Series accession to MAGE-TAB IDF and SDRF.
- `geo2json`: GEO Series accession to parsed MINiML-compatible JSON.
- `json2ae`: parsed JSON to MAGE-TAB IDF and SDRF.
- `ae2json`: local, HTTP(S), or BioStudies MAGE-TAB to parsed JSON.
- `json2h5ad`: parsed JSON plus H5AD, matrix, or FASTQ assets to normalized H5AD.

## Installation

Install the base package from GitHub:

```bash
python -m pip install "git+https://github.com/jaychowcl/meta_standards_converter.git"
```

Install locally for development:

```bash
git clone https://github.com/jaychowcl/meta_standards_converter
cd meta_standards_converter
python -m pip install -e .
```

Include AnnData/H5AD support when using `json2h5ad`:

```bash
python -m pip install -e '.[h5ad]'
```

Build the project image, which includes the H5AD extra, Java 21, Nextflow, `gffread`, and the Docker CLI:

```bash
docker build -t meta-standards-converter .
```

### Requirements

- Python `>=3.10`.
- Base dependencies: `requests>=2.31.0` and `python-dateutil>=2.8.2`.
- H5AD dependencies: AnnData, h5py, NumPy, pandas, Scanpy, and SciPy; install the `h5ad` extra.
- Network access for live GEO, BioStudies, PubMed, NCBI SRA, and ENA lookups.
- Host-side raw FASTQ processing: Java, Nextflow, and a supported Nextflow runtime/profile such as Docker or Apptainer.
- GFF/GFF3 annotation conversion: `gffread`.
- Rootless Compose processing: Linux, Docker Engine rootless extras, subordinate UID/GID support, ACL tools, and user-level systemd.

The Python metadata converters do not require Docker. The project image supplies the scientific and workflow dependencies needed by `json2h5ad`, but raw Docker-profile processing also requires access to a Docker daemon.

## Quickstart

### CLI quickstart

Install the package, then run any of its five commands. This example creates parsed JSON and then normalized H5AD. See the [CLI guide](#cli).

```bash
geo2json GSE234602 --out output
json2h5ad output/GSE234602.json --out output
```

### Python API quickstart

Import a converter and call `convert()`. See the [Python API guide](#python-api).

```python
from meta_standards_converter.converters.geo2json import geo2json

packages = geo2json().convert("GSE234602", out="output")
```

### Docker quickstart

Build the image and mount a writable output directory. See the [Docker guide](#docker).

```bash
docker build -t meta-standards-converter .
mkdir -p output
docker run --rm -v "$PWD/output:/out" \
  meta-standards-converter geo2ae GSE234602 --out /out
```

### Rootless Docker Compose quickstart

Provision the dedicated runner once, then build through its rootless daemon. See the [Rootless Docker Compose guide](#rootless-docker-compose).

```bash
sudo "$PWD/scripts/provision-rootless-json2h5ad.sh" "$PWD" "$PWD/.out/json2h5ad"
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" build converter
```

### Inputs & Outputs

| Workflow | Expected input | Output |
| --- | --- | --- |
| `geo2ae` | One or more `GSE...` accessions | `{accession}.idf.txt` and `{accession}.sdrf.txt`; Python returns MAGE-TAB row payloads |
| `geo2json` | One or more `GSE...` accessions | `{GSE}.json`; Python returns a package `list[dict]` |
| `json2ae` | JSON containing one package object or a non-empty package list | IDF/SDRF files; Python returns ordered MAGE-TAB payloads |
| `ae2json` | IDF path, HTTP(S) IDF URL, or BioStudies/ArrayExpress accession; optional SDRF overrides | `{accession}.json`; Python returns a one-package list with a `mage_tab` extension |
| `json2h5ad` | Non-empty parsed package-list JSON plus discovered or explicit H5AD, matrix, or FASTQ assets | Per-sample H5ADs, optional compatible combined H5AD, provenance JSON, and optional nf-core results |

GEO JSON packages contain Series metadata plus the referenced samples, platforms, contributors, organizations, and databases. MAGE-TAB-origin JSON uses the same public package shape and adds `mage_tab.model`, warnings, unmapped data, and lossless round-trip metadata. H5AD outputs retain expression values, normalized `msc_*` observation metadata, flattened MINiML metadata in `uns["msc_miniml"]`, and conversion provenance.

## Guide

### CLI

The package installs `geo2ae`, `geo2json`, `json2ae`, `ae2json`, and `json2h5ad`. Run `<command> --help` for generated usage text.

All commands process multiple positional inputs in order. A failed input is logged, later inputs continue, and the final exit status is `1`; a fully successful invocation returns `0`. Logging defaults to `WARNING`. `-v` selects `INFO`, `-vv` selects `DEBUG`, and `-q` selects `ERROR`.

#### `geo2ae`

Fetch one or more GEO Series and write MAGE-TAB IDF/SDRF files.

```bash
geo2ae GSE234602 --out output
geo2ae GSE234602 GSE34779 --related --keep-empty --out output
```

| Argument | Behavior |
| --- | --- |
| `gse` | One or more GEO Series accessions; required. |
| `-h`, `--help` | Display generated help and exit. |
| `--related`, `--related-series`, `--get-related-series` | Include transitively related GEO super/subseries; disabled by default. |
| `--remove-empty` | Remove empty parsed fields; this is the default. |
| `--keep-empty` | Preserve empty parsed fields; mutually exclusive with `--remove-empty`. |
| `--out` `OUT` | Output directory; default `.`. |
| `-v`, `--verbose` | Increase verbosity; repeat as `-vv` for DEBUG. |
| `-q`, `--quiet` | Emit ERROR logs only; mutually exclusive with verbosity. |
| `--log-file` `LOG_FILE` | Also write logs to this file, replacing an existing file. |

#### `geo2json`

Fetch one or more GEO Series and write parsed MINiML-compatible JSON package lists.

```bash
geo2json GSE234602 --out output
geo2json GSE234602 --no-enrich --keep-empty --out output
```

| Argument | Behavior |
| --- | --- |
| `gse` | One or more GEO Series accessions; required. |
| `-h`, `--help` | Display generated help and exit. |
| `--related`, `--related-series`, `--get-related-series` | Include transitively related GEO super/subseries; disabled by default. |
| `--remove-empty` | Remove empty parsed fields; this is the default. |
| `--keep-empty` | Preserve empty parsed fields; mutually exclusive with `--remove-empty`. |
| `--no-enrich` | Skip PubMed and SRA/ENA enrichment; enrichment is enabled by default. |
| `--out` `OUT` | Output directory; default `.`. |
| `-v`, `--verbose` | Increase verbosity; repeat as `-vv` for DEBUG. |
| `-q`, `--quiet` | Emit ERROR logs only; mutually exclusive with verbosity. |
| `--log-file` `LOG_FILE` | Also write logs to this file, replacing an existing file. |

#### `json2ae`

Read parsed JSON and write MAGE-TAB IDF/SDRF files.

```bash
json2ae output/GSE234602.json --out output
json2ae primary.json related.json --no-enrich --out output
```

| Argument | Behavior |
| --- | --- |
| `json_path` | One or more JSON paths; each must contain one package object or a non-empty package list. |
| `-h`, `--help` | Display generated help and exit. |
| `--no-enrich` | Convert supplied metadata without PubMed/SRA enrichment; enrichment is enabled by default. |
| `--out` `OUT` | Output directory; default `.`. |
| `-v`, `--verbose` | Increase verbosity; repeat as `-vv` for DEBUG. |
| `-q`, `--quiet` | Emit ERROR logs only; mutually exclusive with verbosity. |
| `--log-file` `LOG_FILE` | Also write logs to this file, replacing an existing file. |

`json2ae` validates all packages before converting any of them. If the input came from `ae2json`, an unchanged single-SDRF package can reproduce its original tables exactly; edits to the typed `mage_tab.model` or mapped core fields regenerate the relevant MAGE-TAB content.

#### `ae2json`

Resolve an IDF and its SDRFs, then write a MINiML-compatible JSON package.

```bash
ae2json study.idf.txt --out output
ae2json https://example.org/study.idf.txt --out output
ae2json E-MTAB-1990 --out output
ae2json study.idf.txt --sdrf first.sdrf.txt --sdrf second.sdrf.txt --out output
```

| Argument | Behavior |
| --- | --- |
| `source` | One or more IDF paths, HTTP(S) IDF URLs, or BioStudies/ArrayExpress accessions. |
| `-h`, `--help` | Display generated help and exit. |
| `--sdrf` `PATH_OR_URL` | Override IDF SDRF references; repeat for multiple SDRFs. Requires exactly one `source` and cannot accompany an accession source. |
| `--out` `OUT` | Output directory; default `.`. |
| `-v`, `--verbose` | Increase verbosity; repeat as `-vv` for DEBUG. |
| `-q`, `--quiet` | Emit ERROR logs only; mutually exclusive with verbosity. |
| `--log-file` `LOG_FILE` | Also write logs to this file, replacing an existing file. |

Remote IDF/SDRF text remains in memory. Accession mode uses BioStudies to discover exactly one IDF and at least one SDRF; assay data files are not downloaded.

#### `json2h5ad`

Select the best available expression source for every sample, normalize it into AnnData, and write H5AD outputs. Explicit manifest assets outrank `--asset` entries, which outrank JSON-discovered assets; within a source tier the order is H5AD, matrix, then raw FASTQ.

```bash
json2h5ad output/GSE234602.json --out output
json2h5ad output/GSE234602.json \
  --asset GSM9651991=local.h5ad --out output
json2h5ad output/GSE234602.json \
  --force-reprocess --pipeline rnaseq --genome GRCh38 \
  --gtf references/current.gtf.gz --profile docker --out output
```

| Argument | Behavior |
| --- | --- |
| `json_path` | One or more paths containing a non-empty package list. |
| `-h`, `--help` | Display generated help and exit. |
| `--out` `OUT` | Output directory; default `.`. |
| `--asset-manifest` `ASSET_MANIFEST` | CSV/TSV mapping with required `scope_id` and `path` columns and optional kind, role, read/lane, matrix, checksum, and orientation metadata. |
| `--asset` `ACCESSION=PATH` | Explicit local or remote H5AD, matrix, or FASTQ; repeat as needed. |
| `--force-reprocess` | Ignore processed sources and require raw FASTQ for every sample. |
| `--pipeline` `{auto,scrnaseq,rnaseq}` | Raw-input pipeline; default `auto`, which groups samples by detected modality. |
| `--genome` `GENOME` | nf-core catalogue genome key, optionally combined with `--gtf` or `--gff`. |
| `--fasta` `FASTA` | Local custom FASTA; requires exactly one of `--gtf` or `--gff`. |
| `--gtf` `GTF` | Local GTF or GTF.GZ annotation; mutually exclusive with `--gff`. |
| `--gff` `GFF` | Local GFF/GFF3 annotation; mutually exclusive with `--gtf` and converted to GTF with `gffread`. |
| `--accept-inferred-reference` | Accept supported human/mouse reference inference when no explicit reference is supplied. |
| `--profile` `PROFILE` | Nextflow profile; default `docker`. |
| `--revision` `REVISION` | Override the pinned nf-core revision; defaults are `scrnaseq` 4.2.0 and `rnaseq` 3.26.0. |
| `--params-file` `PARAMS_FILE` | Additional nf-core JSON parameters; converter-owned input, output, and reference values take precedence. |
| `--nextflow-config` `NEXTFLOW_CONFIG` | Additional Nextflow resource/infrastructure config. |
| `--work-dir` `WORK_DIR` | Nextflow work directory; defaults below the study/pipeline output tree. |
| `--resume` | Add `-resume` to the Nextflow invocation. |
| `--overwrite` | Replace normalized H5AD and manifest outputs; existing outputs are protected by default. |
| `--matrix-orientation` `{auto,genes-by-observations,observations-by-genes}` | Delimited matrix orientation; default `auto`, which rejects ambiguous generic matrices. |
| `-v`, `--verbose` | Increase verbosity; repeat as `-vv` for DEBUG. |
| `-q`, `--quiet` | Emit ERROR logs only; mutually exclusive with verbosity. |
| `--log-file` `LOG_FILE` | Also write logs to this file, replacing an existing file. |

Processed assets may be local or HTTP(S)/FTP and may include `.h5ad`, `.h5ad.gz`, 10x HDF5, 10x MTX directories, CSV, TSV, or TXT matrices. Remote processed assets are cached under the output directory and an available MD5 is verified. Raw processing upgrades known ENA/NCBI FTP FASTQ links to HTTPS before writing nf-core samplesheets.

Each successful sample produces `{GSM}.h5ad`. Compatible samples are outer-joined into `{GSE}.h5ad`; incompatible organisms, references, modalities, or feature namespaces leave the sample files intact, omit the combined file, record a partial failure, and cause CLI status `1`. Every run writes `{GSE}.json2h5ad.json` provenance unless output protection rejects an existing file.

### Python API

The converters accept injectable collaborators for testing and integration, but default construction is sufficient for normal use.

Convert GEO to MAGE-TAB:

```python
from meta_standards_converter.converters.geo2ae import geo2ae

magetabs = geo2ae().convert(
    gse="GSE234602",
    related_series=False,
    remove_empty=True,
    out="output",
)
```

`geo2ae.convert(gse, related_series=False, remove_empty=True, out=None)` returns a list of in-memory MAGE-TAB payloads. `out=None` suppresses file writes.

Convert GEO to JSON:

```python
from meta_standards_converter.converters.geo2json import geo2json

packages = geo2json().convert(
    gse="GSE234602",
    related_series=False,
    remove_empty=True,
    enrich=True,
    out="output",
)
```

`geo2json.convert(gse, related_series=False, remove_empty=True, enrich=True, out=None)` returns `list[dict]`; `out` writes `{gse}.json`.

Convert parsed JSON to MAGE-TAB:

```python
from meta_standards_converter.converters.json2ae import json2ae

magetabs = json2ae().convert(
    json_path="output/GSE234602.json",
    out="output",
    enrich=True,
)
```

`json2ae.convert(json_path, out=None, enrich=True)` accepts a package object or package list and returns ordered MAGE-TAB payloads.

Convert MAGE-TAB to parsed JSON:

```python
from meta_standards_converter.converters.ae2json import ae2json

packages = ae2json().convert(
    source="E-MTAB-1990",
    out="output",
    sdrf_sources=None,
)
```

`ae2json.convert(source, out=None, sdrf_sources=None)` returns a one-package list. `sdrf_sources` is a list of explicit local paths or HTTP(S) URLs and follows the same constraints as repeated CLI `--sdrf` values.

Convert parsed JSON and expression assets to H5AD:

```python
from meta_standards_converter.converters.json2h5ad import Asset, json2h5ad

result = json2h5ad().convert(
    json_path="output/GSE234602.json",
    out="output",
    explicit_assets=[Asset("GSM9651991", "local.h5ad", "h5ad")],
    asset_manifest=None,
    asset_specs=None,
    force_reprocess=False,
    matrix_orientation="auto",
    overwrite=False,
    pipeline="auto",
    genome="GRCh38",
    fasta=None,
    gtf="references/current.gtf.gz",
    gff=None,
    accept_inferred_reference=False,
    profile="docker",
    revision=None,
    params_file=None,
    nextflow_config=None,
    work_dir=None,
    resume=False,
)
```

`JSON2H5ADConverter.convert()` returns `ConversionResult`, whose principal fields are `study_accession`, `sample_h5ads`, `combined_h5ad`, `retained_h5ads`, `pipeline_runs`, `manifest_path`, `warnings`, `failures`, `primary_h5ad`, and `partial`. Paths returned in memory are absolute; persisted provenance paths are relative to their artifact parent where possible.

### Docker

Build the image:

```bash
docker build -t meta-standards-converter .
```

With no command, the image displays `geo2ae --help`. Any installed CLI can be supplied after the image name:

```bash
docker run --rm meta-standards-converter geo2json --help
```

Mount host paths for inputs and outputs. Use matching container paths in CLI arguments:

```bash
mkdir -p output
docker run --rm \
  -v "$PWD/output:/work" \
  meta-standards-converter \
  geo2json GSE234602 --out /work

docker run --rm \
  -v "$PWD/output:/work" \
  meta-standards-converter \
  json2ae /work/GSE234602.json --out /work
```

The standard image contains no Docker daemon. Metadata conversion and processed-asset H5AD conversion work without a nested runtime. Raw FASTQ processing with the Docker profile requires a deliberately supplied daemon; use the hardened rootless Compose workflow below.

### Rootless Docker Compose

The rootless workflow is intended for raw `json2h5ad` processing. It creates a locked `nfcore-runner` account, gives it read access to the project and read/write access only to `.out/json2h5ad`, and connects the converter to that account's rootless Docker socket.

Provision once as root:

```bash
sudo "$PWD/scripts/provision-rootless-json2h5ad.sh" \
  "$PWD" "$PWD/.out/json2h5ad"
```

Build and verify the image as the runner:

```bash
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" build converter
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" \
  run --rm converter docker info --format '{{json .SecurityOptions}}'
```

Generate JSON and process raw data. All mounted inputs, outputs, caches, and Nextflow work must remain under `.out/json2h5ad`:

```bash
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" \
  run --rm converter geo2json GSE104830 \
  --out "$PWD/.out/json2h5ad/json" -vv

sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" \
  run --rm converter json2h5ad \
  "$PWD/.out/json2h5ad/json/GSE104830.json" \
  --out "$PWD/.out/json2h5ad/bulk" \
  --force-reprocess --pipeline rnaseq \
  --accept-inferred-reference --profile docker -vv
```

The helper refuses non-rootless daemons. Compose drops all capabilities, enables `no-new-privileges`, makes the root filesystem read-only, and mounts only the dedicated output tree and rootless socket. Final H5AD files use mode `0660`; provisioning establishes and verifies the project-owner and runner ACLs.

### Code flow

```text
CLI or Python API
  |
  +-- GEO accession
  |     -> GEOWebFetcher -> GEOParser -> [MINiMLEnricher]
  |          |                                  |
  |          +-> geo2json: JSON packages        +-> PubMed / NCBI SRA / ENA
  |          `-> geo2ae: AEConstructor -> IDF + SDRF
  |
  +-- parsed JSON
  |     +-> json2ae: validate -> [enrich] -> AEConstructor -> IDF + SDRF
  |     `-> json2h5ad: plan assets -> [nf-core for FASTQ]
  |                         -> normalize AnnData -> sample/combined H5AD + manifest
  |
  `-- IDF path, URL, or BioStudies accession
        -> AEWebFetcher -> AEParser -> JSON package + mage_tab sidecar
```

Network requests pass through service-specific rate limiting, timeouts, and retries. CLI entrypoints catch failures per top-level input, while programmatic converter calls raise errors to their caller. For class-level call graphs, branches, external API operations, and data shapes, use the canonical codebase documentation below.

## Docs

- [Docs index](docs/index.md): routing index with stable anchors, section purposes, and keywords.
- [Codebase docs](docs/codebase.md): canonical architecture, workflow, callable, test, and maintenance handoff.

## Authors

Created by [jaychowcl](https://github.com/jaychowcl) @ [Saez-Rodriguez Group](https://saezlab.org) & [EMBL-EBI Functional Genomics Team](https://www.ebi.ac.uk/about/teams/functional-genomics/) on May 2026
