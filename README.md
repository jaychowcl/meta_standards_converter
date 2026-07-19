<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/51b52963-19de-4f67-8977-072b409dae19" />

# meta_standards_converter

Convert GEO Series MINiML metadata into parsed JSON, AnnData/H5AD, and ArrayExpress/MAGE-TAB files.

## Description

`meta_standards_converter` is a Python package for converting biological repository metadata between GEO MINiML, parsed JSON, and ArrayExpress-style MAGE-TAB outputs. The implemented workflow fetches GEO Series MINiML XML, parses each Series into a self-contained metadata package, optionally enriches it with PubMed and SRA/ENA metadata, and renders either JSON or IDF/SDRF tables.

Main features:

- `geo2ae`: GEO Series accession to MAGE-TAB IDF and SDRF TSV files.
- `geo2json`: GEO Series accession to parsed MINiML JSON packages, enriched by default.
- `json2ae`: parsed MINiML JSON package or package list to MAGE-TAB IDF and SDRF TSV files.
- `ae2json`: local, HTTP(S), or BioStudies/ArrayExpress MAGE-TAB IDF and SDRF metadata to MINiML-compatible JSON.
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

Build the Docker image, which includes Java 21, Nextflow, the Docker CLI, and the
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
  - BioStudies API and HTTP file service for accession-based MAGE-TAB metadata
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

CLI JSON-to-MAGE-TAB conversion, see [CLI guide](#cli):

```bash
json2ae output/GSE234602.json --out output
```

CLI MAGE-TAB-to-JSON conversion, see [CLI guide](#cli):

```bash
ae2json E-MTAB-1990 --out output
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
| `json2ae` | One or more parsed MINiML JSON files containing a package object or non-empty package list | `{accession}.idf.txt` and `{accession}.sdrf.txt` files; Python API returns MAGE-TAB row payloads |
| `ae2json` | IDF path, HTTP(S) IDF URL, or BioStudies/ArrayExpress accession; optional explicit SDRF paths/URLs | `{accession}.json`; Python API returns `list[dict]` MINiML-compatible packages |
| `json2h5ad` | Parsed MINiML JSON plus discovered or explicit H5AD, matrix, or FASTQ assets | Normalized per-sample H5ADs, compatible combined study H5AD, provenance manifest, and optional nf-core results |
| `GEOParser` | MINiML XML string | One self-contained package per Series |
| `MINiMLEnricher` | Parsed package dict | Same dict with PubMed and SRA/ENA fields added where available |
| `AEConstructor` | Enriched parsed package dict | In-memory MAGE-TAB rows or IDF/SDRF TSV files |

## Guide

### CLI

The package installs five console scripts: `geo2ae`, `geo2json`, `json2ae`, `ae2json`, and `json2h5ad`.

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

`json2ae` reads parsed MINiML JSON and writes MAGE-TAB IDF/SDRF files:

```bash
json2ae output/GSE234602.json --out output
json2ae parsed-primary.json parsed-related.json --no-enrich --out output
```

Options:

- `JSON...`: one or more parsed MINiML JSON files. Each file may contain one package object or a non-empty package list.
- `--out DIR`: output directory, default `.`.
- `--no-enrich`: skip PubMed/SRA enrichment and convert the supplied metadata exactly as provided.
- `-v`, `-vv`, `-q`, `--log-file PATH`: shared logging options.

`ae2json` reads an IDF plus its referenced SDRF files and writes parsed JSON:

```bash
ae2json study.idf.txt --out output
ae2json https://example.org/study.idf.txt --out output
ae2json E-MTAB-1990 --out output
ae2json study.idf.txt --sdrf first.sdrf.txt --sdrf second.sdrf.txt --out output
```

Options:

- `SOURCE...`: one or more IDF paths, HTTP(S) IDF URLs, or BioStudies/ArrayExpress accessions.
- `--sdrf PATH_OR_URL`: override IDF SDRF references; repeat for multiple SDRFs and use with exactly one source.
- `--out DIR`: output directory, default `.`.
- `-v`, `-vv`, `-q`, `--log-file PATH`: shared logging options.

Remote IDF/SDRF text is parsed in memory. Accession lookup downloads metadata only; it does not download referenced assay data files. The `mage_tab.roundtrip` sidecar retains source tables and a semantic fingerprint, while unmapped rows and columns remain available with warnings.
The core JSON keeps the fixed GEO/MINiML fields for samples, channels, runs,
files, factors, databases, platforms, contributors, and supported protocol
slots. Arbitrary MAGE-TAB protocol graphs, assay identities, performers,
protocol hardware/software, QC/replicate declarations, per-value units and
ontology annotations, and custom rows or columns remain sidecar-only; removing
`mage_tab` is therefore not a lossless operation.

`json2h5ad` selects the best source per sample (`H5AD > matrix > raw`):

Local and remote H5AD assets may be supplied either uncompressed (`.h5ad`) or
gzip-compressed (`.h5ad.gz`), as commonly published by GEO.

```bash
json2h5ad output/GSE234602.json --out output
json2h5ad output/GSE234602.json --asset GSM9651991=local.h5ad --out output
json2h5ad output/GSE234602.json --force-reprocess --genome GRCh38 --profile docker --out output
json2h5ad output/GSE234602.json --force-reprocess --genome GRCh38 --gtf current.gtf.gz --out output
json2h5ad output/GSE234602.json --force-reprocess --fasta genome.fa.gz --gff genes.gff3.gz --out output
```

Options:

- `JSON...`: one or more parsed MINiML JSON files.
- `--out DIR`: output directory, default `.`.
- `--asset ACCESSION=PATH`: explicit local or remote asset; repeat as needed.
- `--asset-manifest CSV`: detailed sample/study asset mapping, including matrix bundles and FASTQ read/lane metadata.
- `--matrix-orientation`: required for ambiguous delimited matrices.
- `--force-reprocess`: ignore processed sources and rebuild from raw FASTQs.
- `--pipeline auto|scrnaseq|rnaseq`: automatic per-sample routing or an explicit nf-core pipeline.
- `--genome`: select a catalogue reference; combine it with `--gtf` or `--gff` to override the catalogue annotation.
- `--fasta` with exactly one of `--gtf` or `--gff`: fully user-supplied reference sequence and annotation. Inference requires `--accept-inferred-reference`.
- `--profile`, `--revision`, `--params-file`, `--nextflow-config`, `--work-dir`, and `--resume`: nf-core/Nextflow controls.
- `--overwrite`: replace normalized outputs; outputs are protected by default.
- `-v`, `-vv`, `-q`, `--log-file PATH`: shared logging options.

Known ENA and NCBI `ftp://` FASTQ URLs are written to nf-core samplesheets
using the archives' equivalent `https://` endpoints, which are more reliable
through rootless container networking.

For `nf-core/scrnaseq`, H5AD discovery covers the complete results tree and
selects CellBender-filtered, then filtered (including QCATCH
`*_filtered_quants.h5ad`), then raw output for each sample.

Final H5ADs expose converter-owned observation metadata through both the
established `msc_*` names and additive dotted aliases grouped by section, such
as `msc.sample.accession`, `msc.sample.channel.organism.value`,
`msc.library.strategy`, and `msc.database.identifier`. Arbitrary MINiML
characteristics appear as both `msc_characteristic_*` and
`msc.characteristics.*`; this includes harmonized tags such as
`msc.characteristics.hz_cell_type`, its ID, and its ontology. Existing source
H5AD observation columns remain intact. Single-cell values are sample
annotations repeated across the sample's cells, and combined studies expose
both `msc_batch` and `msc.combination.batch`.

Canonical organism labels prefer each channel's harmonized `hz_organism`
value, including the harmonizer's container-list and scalar forms, and fall
back to the channel's original `organism` value. Channels without either value
remain empty; the converter does not infer an organism from unrelated fields.

The permitted remaining MINiML metadata is stored as a typed flattened table
in `uns["msc_miniml"]["fields"]`. GSM files contain the relevant sample and
transitively referenced platform, contributor, and database records; GSE files
contain all study records. Dictionary-shaped MINiML references such as
`{"ref": "GPL..."}` are resolved directly. Repository identity
is taken from the MINiML `database` record rather than hard-coded. Publication
storage is restricted to citation metadata (PubMed ID, DOI, title, authors,
status, and ontology identifiers); abstracts and publication full text are
never embedded.

Persisted H5AD and manifest provenance uses paths relative to the containing
artifact and marks each path as internal, external, or remote. Python
`ConversionResult` paths remain absolute for immediate programmatic use.
Nextflow warnings are retained in both pipeline-run and top-level manifest
diagnostics without turning a successful pipeline run into a failure.

User annotations must be local `.gtf[.gz]`, `.gff[.gz]`, or `.gff3[.gz]`
files. GFF3 is converted once with `gffread` to a checksum-addressed GTF shared
by bulk and single-cell runs. Annotation source, format, SHA-256, and effective
GTF are recorded in H5AD and JSON provenance. FASTA/annotation build and
chromosome-name compatibility remain the caller's responsibility.

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

Convert parsed JSON to MAGE-TAB:

```python
from meta_standards_converter.converters.json2ae import json2ae

magetabs = json2ae().convert(
    "output/GSE234602.json",
    enrich=True,
    out="output",
)
```

Convert MAGE-TAB to parsed JSON:

```python
from meta_standards_converter.converters.ae2json import ae2json

packages = ae2json().convert(
    "E-MTAB-1990",
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
    gtf="references/current.gtf.gz",
)
print(result.combined_h5ad)
print(result.sample_h5ads)
```

Install the `h5ad` extra before using processed-matrix conversion. Raw workflows
run directly on the host additionally require Nextflow, Java, and the selected
execution profile. GFF3 input additionally requires `gffread`. The project
image already contains Java, Nextflow, `gffread`, and the Docker CLI.

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
not grant it access to secrets or unrelated host directories. Generated H5ADs
use mode `0660`, allowing the provisioned output ACL to give the invoking
project user read/write access without making results world-readable. The
provisioning and Compose helpers verify effective ACL access before and after
runs. Filesystems that map rootless writers to `nobody:nogroup` remain
supported when the project owner and runner retain effective access; no
privileged ownership repair is performed. General
temporary files remain on a `noexec` `/tmp`; only Nextflow's separate 2 GiB
`/nextflow-tmp` is executable so its AWS/S3 native client can load extracted
libraries needed for iGenomes references.

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

`json2ae.convert` resumes the MAGE-TAB path from persisted parsed JSON:

```text
load one package object or a non-empty package list
validate every package has a usable study accession
for each package:
  optionally enrich PubMed and SRA/ENA metadata (default)
  build MAGE-TAB with AEConstructor.miniml2magetab
optionally write IDF/SDRF files with AEConstructor.magetab2file
return MAGE-TAB payloads in input order
```

`ae2json.convert` performs the reverse metadata conversion:

```text
resolve an IDF path, HTTP(S) URL, or BioStudies accession
resolve referenced SDRFs, or use explicit SDRF overrides
parse known IDF rows and SDRF graph columns into the existing package shape
merge samples and platforms across SDRFs in first-seen order
preserve unmapped rows/columns and conversion warnings under mage_tab
retain source IDF/SDRF rows and a semantic fingerprint under mage_tab.roundtrip
optionally write a one-package JSON list as {accession}.json
return the one-package list
```

`json2h5ad.convert` plans and normalizes expression assets. Merged nf-core/rnaseq
tables may include the standard text `gene_name` column; it is retained in
`adata.var`, while per-sample counts and TPM values remain numeric:

```text
load parsed package list and explicit asset mappings
for each sample, select H5AD > matrix > raw FASTQ
if raw:
  group samples into nf-core/scrnaseq or nf-core/rnaseq
  resolve a catalogue or user FASTA reference and optional user annotation
  normalize GFF3 to a checksum-addressed GTF with gffread
  run pinned Nextflow and discover outputs
normalize each sample into sparse AnnData with msc_* metadata and provenance
flatten permitted MINiML fields into uns["msc_miniml"]
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
- `AEWebFetcher`: resolves local and HTTP(S) IDFs plus SDRFs, or discovers an accession's metadata through BioStudies; remote metadata stays in memory.
- `AEParser`: maps general MAGE-TAB IDF/SDRF metadata into the MINiML-compatible package shape and preserves unmapped content under `mage_tab`.
- `mage_tab.roundtrip`: optional versioned sidecar used to reproduce unchanged AE-origin tables; edited JSON wins and unsupported metadata is restored where it can be matched safely.
- `RateLimitedRequester`: wraps `requests.get` with service-specific timeout, request delay, retry status handling, and backoff.
- `PubmedWebFetcher`: calls NCBI PubMed ESummary and returns DOI, authors, title, and harmonized publication status.
- `INSDCWebfetcher`: calls NCBI SRA EFetch and ENA Portal file reports, then returns run-level FASTQ metadata.
- `GEOParser`: parses XML generically, resolves references, traverses related Series, and removes empty values.
- `MINiMLEnricher`: adds `series.pubmed_publication`, `sample.sra_accession`, `sample.sra_run`, and `sample.ena_accession`.
- `json2ae`: validates persisted parsed packages, optionally enriches them, and delegates IDF/SDRF construction to `AEConstructor`.
- `ae2json`: resolves MAGE-TAB input, delegates parsing to `AEParser`, and optionally writes an accession-named JSON package list.
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
