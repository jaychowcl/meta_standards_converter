<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/51b52963-19de-4f67-8977-072b409dae19" />

# meta_standards_converter

Convert biological study metadata among GEO MINiML, JSON, MAGE-TAB, sample tables, and AnnData/H5AD.

## Installation

Create a Conda environment and install MSC:

1. `conda create -n msconverter -c conda-forge python=3.12 pip`
2. `conda activate msconverter`
3. `pip install "git+https://github.com/jaychowcl/meta_standards_converter.git"`

Activate the environment with `conda activate msconverter` in each new terminal.

## Quickstart

Convert a GEO study to MAGE-TAB IDF and SDRF files:

`msc-convert GSE234602 --out-type magetab --out magetabs`

Replace the accession with your own. The generated files are written under
`magetabs/`; retrieving public studies requires network access.

Show all terminal options:

`msc-convert --help`

`msc-convert` is the CLI for `Converter.convert()`. It accepts accessions,
supported files, URLs and directories. Select a required output format with
`--out-type`: `json`, `magetab`, `tsv`, `csv`, `h5ad` or `obs`.

Defaults:

- Input detection is automatic; use `--in-type` to select a reader.
- Enrichment is `standard`; select `curators` or `off` with `--enrichment`.
- Verified related-study expansion is enabled, independently of enrichment;
  use `--no-expand-studies` to disable it.
- Output goes to the current directory unless `--out` or `--outfile` is supplied.

For other repositories and platform selection, follow the
[curator guide](docs/curators-guide.md). For H5AD/OBS dependencies and processing
requirements, see the [runtime guide](docs/codebase.md#runtime-behavior).

## Docs

- [Curator guide](docs/curators-guide.md): installation, repository examples,
  main options and platform selection.
- [Full CLI options](docs/codebase.md#unified-cli): all flags, summaries,
  JSON reports and exit codes.
- [Python API](docs/codebase.md#unified-converter): `Converter.convert()`,
  supported inputs and outputs, and advanced settings.
- [Docs index](docs/index.md): guides and references by topic.
- [Codebase docs](docs/codebase.md): canonical architecture, legacy commands,
  Python interfaces, Docker setup and testing instructions.

## Authors

Created by [jaychowcl](https://github.com/jaychowcl) @ [Saez-Rodriguez Group](https://saezlab.org) & [EMBL-EBI Functional Genomics Team](https://www.ebi.ac.uk/about/teams/functional-genomics/) on May 2026
