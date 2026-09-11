<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/51b52963-19de-4f67-8977-072b409dae19" />

# meta_standards_converter

Convert biological study metadata among GEO MINiML, JSON, MAGE-TAB, sample tables, and AnnData/H5AD.

## Description

MSC is a Python library and command-line toolkit for fetching study metadata,
converting metadata formats, and attaching metadata to expression data. It can
read Atlas documents without installing ThematicAtlases.

| Command | Use it to… |
| --- | --- |
| `geo2json` | Fetch a GEO Series and produce MSC MINiML JSON |
| `geo2ae` | Fetch a GEO Series and produce MAGE-TAB IDF/SDRF files |
| `ae2json` | Read local or remote MAGE-TAB and produce MSC MINiML JSON |
| `json2ae` | Convert MSC MINiML or Atlas JSON to MAGE-TAB |
| `json2tsv` | Export sample-level metadata as TSV or CSV |
| `json2h5ad` | Produce a per-sample H5AD catalogue from metadata and expression assets |
| `json2obs` | Export cell observation metadata and optional `var`/`uns` components |
| `miniml-migrate` | Import supported unversioned/1.0 legacy source JSON as v3 |

## Installation

Install from a source checkout in a virtual environment:

```bash
git clone https://github.com/jaychowcl/meta_standards_converter.git
cd meta_standards_converter
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
```

For `json2h5ad` or `json2obs`, install the scientific dependencies:

```bash
python -m pip install '.[h5ad]'
```

For development, use `python -m pip install -e '.[test]'`.

### Requirements

- Python **3.10 or newer**.
- Network access for GEO/BioStudies retrieval and optional PubMed/SRA/ENA enrichment.
- The `h5ad` extra for expression and AnnData workflows.
- For raw FASTQ processing: Java, Nextflow, a supported execution runtime,
  and a reference genome/annotation. GFF conversion also requires `gffread`.

Metadata conversion does not require Docker. See [runtime requirements](docs/codebase.md#runtime-behavior)
for dependency bounds and [raw processing](docs/codebase.md#reference-annotation-flow)
for reference configuration.

## Quickstart

### CLI quickstart

Fetch a study and export its sample metadata (retrieval requires network access):

```bash
geo2json GSE234602 --out output
json2tsv output/GSE234602.json --out tables
```

For an existing package, `json2tsv` performs metadata-only conversion without
fetching expression assets. See the [CLI guide](#cli).

### Python API quickstart

```python
from meta_standards_converter.converters import JSON2TSVConverter

result = JSON2TSVConverter().export_manifest(
    "output/GSE234602.json", outdir="tables"
)
print(result.output_path)
```

See the [Python API guide](#python-api).

### Docker quickstart

```bash
docker build -t meta-standards-converter .
docker run --rm meta-standards-converter geo2json --help
```

Mount an output directory when converting data; see the [Docker guide](#docker).

### Rootless Docker Compose quickstart

Raw FASTQ workflows can use the dedicated rootless runner. Start with the
[Rootless Docker Compose guide](#rootless-docker-compose) for provisioning,
mounts, and runtime prerequisites.

### Inputs & Outputs

MSC 8 requires **MSC MINiML 3.0** at JSON conversion boundaries. Fresh GEO and
MAGE-TAB ingestion produces v3 directly. Saved v2 packages are rejected;
regenerate them from source or use the preceding release to convert them in a
separate environment. `miniml-migrate` accepts unversioned/1.0 input, **not v2**.
See [migration guidance](docs/codebase.md#miniml-v3-only-cutover).

JSON converters accept a native package, a package list, or an Atlas document
(schema 1.0). Atlas conversion selects harmonized datasets. Metadata converters
produce JSON, IDF/SDRF, or one row per sample in TSV/CSV. H5AD and OBS workflows
also need expression assets; their catalogues do not combine expression matrices.

Harmonized `hz_*` values are exported alongside source values by default.
[Replacement profiles](docs/codebase.md#harmonization-overrides) optionally use
harmonized values in ordinary destination fields while preserving the canonical
input. See [data contracts](docs/codebase.md#data-contracts) for provenance and
round-trip limits.

## Guide

### Configuration

There is no mandatory application configuration file. Use CLI options or Python
arguments. Common controls include output paths, optional enrichment,
`--platform-handler`, resource profiles, and expression asset selection.

Use `json2ae ... --no-enrich` to skip the enrichment stage; MAGE-TAB construction
may still resolve missing publication or sequencing evidence. Use
`--replacement-profile-file policy.json` to activate an export replacement profile.
Expression outputs are protected unless `--overwrite` is supplied; output
behavior varies by command, so consult its reference before reusing a destination.

See [configuration and precedence](docs/codebase.md#configuration).

### CLI

```bash
# Reconstruct MAGE-TAB from an existing v3 package.
json2ae output/GSE234602.json --no-enrich --out mage

# Choose CSV rather than TSV.
json2tsv output/GSE234602.json --format csv --out csv-tables

# Inspect all options for the installed command.
json2h5ad --help
```

The [complete CLI reference](docs/codebase.md#cli) lists all arguments, defaults,
and failure behavior. [Platform handlers](docs/codebase.md#configuration)
explain automatic detection and explicit overrides.

### Python API

Import converters from `meta_standards_converter.converters`. GEO and MAGE-TAB
JSON ingestion returns typed `MINiMLPackage` objects; call `to_mapping()` when
you need a serializable mapping. Python metadata converters generally return
in-memory results when no output path is supplied; manifest and expression
entrypoints have their own output contracts.

See [Python examples](docs/codebase.md#python-api-guide),
[public interfaces](docs/codebase.md#public-api-reference), and
[class relationships and injection points](docs/codebase.md#oop-design).

### Docker

The image includes the H5AD dependencies and raw-processing tools. Supply an
installed command after the image name and mount input/output paths explicitly.
The image does not include a Docker daemon. See [Docker usage](docs/codebase.md#docker-guide).

### Rootless Docker Compose

The supported raw-processing setup uses a dedicated runner, restricted mounts,
and its rootless Docker socket. Provisioning changes the host; follow the
[full setup guide](docs/codebase.md#rootless-json2h5ad-runtime).

### Code flow

Source retrieval → parsing into canonical MINiML → optional enrichment or export
replacement → destination construction/projection → result and optional files.
Expression workflows add asset planning, reading or nf-core processing, and
per-sample checkpointing. See [architecture](docs/codebase.md#architecture) and
[execution flows](docs/codebase.md#principal-workflows).

## Testing

With `.[test]` installed:

```bash
python -m pytest tests/e2e -q
python -m pytest tests/policy tests/test_documented_imports.py -q
python -m pytest -q
```

The default suite blocks external network/process effects; live provider checks
are opt-in. See [test coverage and commands](docs/codebase.md#test-plan),
[fixture provenance](tests/fixtures/README.md), and the historical
[rootless acceptance report](docs/rootless-acceptance-2026-07-31.md).

## Docs

- [Docs index](docs/index.md): find topics by purpose, type, command, or keyword.
- [Codebase docs](docs/codebase.md): architecture, OOP design, workflows, API reference, and maintenance guidance.

## Authors

Created by [jaychowcl](https://github.com/jaychowcl) @ [Saez-Rodriguez Group](https://saezlab.org) & [EMBL-EBI Functional Genomics Team](https://www.ebi.ac.uk/about/teams/functional-genomics/) on May 2026
