<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/51b52963-19de-4f67-8977-072b409dae19" />

# meta_standards_converter

Convert GEO Series MINiML metadata into parsed JSON, AnnData/H5AD, and ArrayExpress/MAGE-TAB files.

## Description

`meta_standards_converter` is a Python package for converting biological repository metadata between GEO MINiML, parsed JSON, and ArrayExpress-style MAGE-TAB outputs. The implemented workflow fetches GEO Series MINiML XML, parses each Series into a self-contained metadata package, optionally enriches it with PubMed and SRA/ENA metadata, and renders either JSON or IDF/SDRF tables.

Main features:

- `geo2ae`: GEO Series accession to MAGE-TAB IDF and SDRF TSV files.
- `geo2json`: GEO Series accession to parsed MINiML JSON packages, enriched by default.
- `json2h5ad`: reuse supplied H5ADs, build AnnData from count matrices, or run pinned nf-core RNA-seq pipelines for raw FASTQs.
- Injectable fetchers, parsers, enrichers, and constructors for testing or downstream integration.
- Related GEO super/subseries traversal when requested.

## Installation

Install from GitHub:

```bash
python -m pip install "git+https://github.com/jaychowcl/meta_standards_converter.git"
```

Install locally for development:

```bash
git clone https://github.com/jaychowcl/meta_standards_converter
cd meta_standards_converter
python -m pip install -e .
python -m pip install -e '.[h5ad]'  # include AnnData conversion support
```

Build the Docker image, which includes Java 17, Nextflow, the Docker CLI, and the
H5AD Python dependencies:

```bash
docker build -t meta-standards-converter .
```

### Requirements

- Python `>=3.10`
- Runtime packages: `requests>=2.31.0`, `python-dateutil>=2.8.2`
- H5AD extra: `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, and `h5py`
- Raw FASTQ processing outside the project image: Nextflow, Java, and a supported container runtime such as Docker or Apptainer
- Network access for live workflows:
  - GEO FTP for MINiML tarballs
  - NCBI E-utilities for PubMed and SRA XML
  - ENA Portal API for FASTQ file reports
- Docker is optional for the Python CLI. The rootless Compose workflow requires
  Docker Engine with rootless extras, subordinate UID/GID support, and user-level
  systemd.

## Quickstart

CLI MAGE-TAB conversion, see [CLI guide](#cli):

```bash
geo2ae GSE234602 --out output
```

CLI JSON conversion, see [CLI guide](#cli):

```bash
geo2json GSE234602 --out output
```

Python API, see [Python guide](#python):

```python
from meta_standards_converter.converters.geo2json import geo2json

packages = geo2json().convert("GSE234602", out="output")
```

Docker CLI, see [Docker guide](#docker):

```bash
docker run --rm meta-standards-converter geo2ae --help
```

### Inputs & Outputs

| Interface | Input | Output |
| --- | --- | --- |
| `geo2ae` | One or more `GSE...` accessions | `{accession}.idf.txt` and `{accession}.sdrf.txt` files; Python API returns MAGE-TAB row payloads |
| `geo2json` | One or more `GSE...` accessions | `{GSE}.json` files; Python API returns `list[dict]` parsed packages |
| `json2h5ad` | Parsed MINiML JSON plus discovered or explicit H5AD, matrix, or FASTQ assets | Normalized per-sample H5ADs, compatible combined study H5AD, provenance manifest, and optional nf-core results |
| `GEOParser` | MINiML XML string | One self-contained package per Series |
| `MINiMLEnricher` | Parsed package dict | Same dict with PubMed and SRA/ENA fields added where available |
| `AEConstructor` | Enriched parsed package dict | In-memory MAGE-TAB rows or IDF/SDRF TSV files |

## Guide

### CLI

The package installs three console scripts: `geo2ae`, `geo2json`, and `json2h5ad`.

`geo2ae` writes MAGE-TAB IDF/SDRF files:

```bash
geo2ae GSE234602 --out output
geo2ae GSE234602 GSE34779 --related --out output
```

Options:

- `GSE...`: one or more GEO Series accessions.
- `--out DIR`: output directory, default `.`.
- `--related`, `--related-series`, `--get-related-series`: include related GEO super/subseries.
- `--remove-empty`: remove empty parsed MINiML fields, default behavior.
- `--keep-empty`: preserve empty parsed MINiML fields.
- `-v`, `--verbose`: log INFO messages.
- `-vv`: log DEBUG messages.
- `-q`, `--quiet`: log only errors.
- `--log-file PATH`: also write logs to a file.

`geo2json` writes parsed package JSON:

```bash
geo2json GSE234602 --out output
geo2json GSE234602 --no-enrich --keep-empty --out output
```

Options are the same as `geo2ae`, plus:

- `--no-enrich`: skip PubMed and SRA/ENA enrichment.

`json2h5ad` selects the best source per sample (`H5AD > matrix > raw`):

```bash
json2h5ad output/GSE234602.json --out output
json2h5ad output/GSE234602.json --asset GSM9651991=local.h5ad --out output
json2h5ad output/GSE234602.json --force-reprocess --genome GRCh38 --profile docker --out output
```

Options:

- `JSON...`: one or more parsed MINiML JSON files.
- `--out DIR`: output directory, default `.`.
- `--asset ACCESSION=PATH`: explicit local or remote asset; repeat as needed.
- `--asset-manifest CSV`: detailed sample/study asset mapping, including matrix bundles and FASTQ read/lane metadata.
- `--matrix-orientation`: required for ambiguous delimited matrices.
- `--force-reprocess`: ignore processed sources and rebuild from raw FASTQs.
- `--pipeline auto|scrnaseq|rnaseq`: automatic per-sample routing or an explicit nf-core pipeline.
- `--genome` or `--fasta` with `--gtf`: reference selection. Inference requires `--accept-inferred-reference`.
- `--profile`, `--revision`, `--params-file`, `--nextflow-config`, `--work-dir`, and `--resume`: nf-core/Nextflow controls.
- `--overwrite`: replace normalized outputs; outputs are protected by default.
- `-v`, `-vv`, `-q`, `--log-file PATH`: shared logging options.

CLI behavior:

- Multiple inputs are processed in order.
- A failed input is logged, later inputs continue, and the command exits with status `1`.
- Fully successful runs exit with status `0`.
- Default CLI logging emits warnings and errors only. Library callers can
  capture structured INFO/DEBUG telemetry for request attempts and durations,
  GEO fetches, MINiML structural counts, related-series progress, and
  enrichment totals; payloads and credentials are not logged.

### Python

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

Convert GEO to parsed JSON:

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

Convert parsed JSON to H5AD:

```python
from meta_standards_converter.converters.json2h5ad import json2h5ad

result = json2h5ad().convert(
    "output/GSE234602.json",
    out="output",
    genome="GRCh38",
)
print(result.combined_h5ad)
print(result.sample_h5ads)
```

Install the `h5ad` extra before using processed-matrix conversion. Raw workflows
run directly on the host additionally require Nextflow, Java, and the selected
execution profile. The project image already contains Java, Nextflow, and the
Docker CLI.

### Docker

Build and inspect the CLI:

```bash
docker build -t meta-standards-converter .
docker run --rm meta-standards-converter
```

Run a conversion and mount an output directory:

```bash
mkdir -p output
docker run --rm -v "$PWD/output:/out" meta-standards-converter geo2ae GSE234602 --out /out
```

#### Rootless nf-core execution with Compose

For raw FASTQs, provision the dedicated `nfcore-runner` account and its rootless
Docker daemon once. This administrative step installs rootless prerequisites,
allocates subordinate IDs, enables the user service, and grants the runner
read/write access only to `.out/json2h5ad`:

```bash
sudo "$PWD/scripts/provision-rootless-json2h5ad.sh" "$PWD" "$PWD/.out/json2h5ad"
```

Build and verify the converter with the rootless daemon:

```bash
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" build converter
sudo -u nfcore-runner -H "$PWD/scripts/json2h5ad-compose.sh" \
  run --rm converter docker info --format '{{json .SecurityOptions}}'
```

Generate enriched JSON and process raw data. All inputs, outputs, Nextflow work,
and caches supplied to the container must remain below `.out/json2h5ad`:

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

The Compose service mounts the dedicated rootless socket, never the system
`/var/run/docker.sock`. `META_STANDARDS_REQUIRE_ROOTLESS_DOCKER=1` makes raw
Docker-profile conversion fail before Nextflow starts if the daemon is
unreachable or not rootless. The account still controls its own daemon, so do
not grant it access to secrets or unrelated host directories.

### Main Code Flow

`geo2ae.convert` orchestrates live GEO to MAGE-TAB conversion:

```text
geo2ae.convert(gse, related_series, remove_empty, out)
  miniml = GEOWebFetcher.fetch_gse_miniml(gse)
  packages = GEOParser.parse(miniml, remove_empty, related_series)
  for package in packages:
    enriched = MINiMLEnricher.enrich(package)
    magetab = AEConstructor.miniml2magetab(enriched)
  if out:
    AEConstructor.magetab2file(magetab, out)
  return magetab payloads
```

`geo2json.convert` shares the fetch and parse stages, optionally calls `MINiMLEnricher.enrich`, writes `{gse}.json` when `out` is set, and returns parsed packages.

`json2h5ad.convert` plans and normalizes expression assets:

```text
load parsed package list and explicit asset mappings
for each sample, select H5AD > matrix > raw FASTQ
if raw:
  group samples into nf-core/scrnaseq or nf-core/rnaseq
  confirm the reference, run pinned Nextflow, and discover outputs
normalize each sample into sparse AnnData with GEO metadata and provenance
combine compatible samples with an outer sparse gene join
write per-sample H5ADs, optional study H5AD, and a JSON provenance manifest
```

`GEOParser.parse` converts MINiML XML into per-Series packages:

```text
parse(miniml)
  root = ElementTree.fromstring(miniml)
  parse top-level Organization, Contributor, Database, Platform, Sample, Series nodes
  build iid indexes
  for each Series:
    include referenced Samples
    include Platforms referenced by included Samples
    include referenced Contributors, Databases, and Organizations
  optionally fetch related GSE accessions from Series relations
  optionally remove empty fields
  return package list
```

`MINiMLEnricher.enrich` mutates one parsed package:

```text
enrich(data)
  enrich_pubmed(data)
    PubmedWebFetcher.pubmed_summary(pubmed_id)
      RateLimitedRequester.get(NCBI PubMed ESummary)
      Harmonizer.pubstatus2efo(status)
  enrich_sra(data)
    extract SRA accessions from sample relations
    INSDCWebfetcher.fetch_sra_runs(accession)
      RateLimitedRequester.get(NCBI SRA EFetch)
      RateLimitedRequester.get(ENA Portal filereport)
  return data
```

`AEConstructor.miniml2magetab` builds MAGE-TAB:

```text
miniml2magetab(data)
  protocol_registry = ProtocolRegistry(series_accession)
  technology_type = _detect_ae_technology(data)
  sdrf = SDRFConstructor._miniml2sdrf(data, protocol_registry, technology_type)
  idf = IDFConstructor.miniml2idf(data, protocol_registry, technology_type)
  replace IDF "SDRF File" placeholder with in-memory SDRF table
  return normalized MAGE-TAB rows
```

Important classes:

- `GEOWebFetcher`: validates `GSE...`, builds the GEO FTP tarball URL, downloads it, and extracts `{GSE}_family.xml`.
- `RateLimitedRequester`: wraps `requests.get` with service-specific timeout, request delay, retry status handling, and backoff.
- `PubmedWebFetcher`: calls NCBI PubMed ESummary and returns DOI, authors, title, and harmonized publication status.
- `INSDCWebfetcher`: calls NCBI SRA EFetch and ENA Portal file reports, then returns run-level FASTQ metadata.
- `GEOParser`: parses XML generically, resolves references, traverses related Series, and removes empty values.
- `MINiMLEnricher`: adds `series.pubmed_publication`, `sample.sra_accession`, `sample.sra_run`, and `sample.ena_accession`.
- `ProtocolRegistry`: assigns stable protocol refs shared by IDF and SDRF.
- `IDFConstructor`: emits investigation, design, person, date, publication, protocol, term-source, and platform-specific IDF rows.
- `SDRFConstructor`: selects a technology handler and renders source/sample/extract/file rows.
- `JSONHandler`: reads dotted JSON paths with list indexes and `*` expansion.
- `Harmonizer`: maps GEO protocol and PubMed status labels to ontology terms where known.
- `MetaStore`: placeholder validation class; structural validation is not implemented.
- `SourcePlanner`, `AssetManifest`, and `AssetDownloader`: select, map, cache, and checksum expression assets.
- `NFCoreRunner` and `ReferenceResolver`: prepare pinned nf-core runs and guard reference selection.
- `JSON2H5ADConverter`: normalize processed/pipeline results and combine compatible AnnData objects.

## Docs

- [Docs index](docs/index.md): read first to select focused `docs/codebase.md` sections.
- [Codebase docs](docs/codebase.md): architecture, workflows, pseudocode, callable references, and test notes.

## Authors

Created by [jaychowcl](https://github.com/jaychowcl) @ [Saez-Rodriguez Group](https://saezlab.org) & [EMBL-EBI Functional Genomics Team](https://www.ebi.ac.uk/about/teams/functional-genomics/) on May 2026
