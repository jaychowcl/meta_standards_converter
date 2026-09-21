<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# Guide to meta_standards_converter

[meta_standards_converter on GitHub](https://github.com/jaychowcl/meta_standards_converter)

Converter between different metadata standards and formats using `msc-convert`, the CLI for `Converter.convert()`.

## Installation

Install MSC:

1. `conda create -n msconverter -c conda-forge python=3.12 pip`
2. `conda activate msconverter`
3. `pip install "git+https://github.com/jaychowcl/meta_standards_converter.git"`

## Converting to MAGE-TAB

Activate MSC before conversion: `conda activate msconverter`.

Replace the example accessions with your own. These commands retrieve public
study metadata and write MAGE-TAB IDF and SDRF files under `magetabs/`.

**GEO → MAGE-TAB**

`msc-convert GSE12345 GSE54321 --in-type geo --out-type magetab --out magetabs`

**ArrayExpress → MAGE-TAB**

`msc-convert E-MTAB-12345 --in-type ae_accession --out-type magetab --out magetabs`

**SRA → MAGE-TAB**

`msc-convert SRP123456 --in-type sra --out-type magetab --out magetabs`

**ENA → MAGE-TAB**

`msc-convert ERP123456 --in-type ena --out-type magetab --out magetabs`

**DDBJ → MAGE-TAB**

Retrieve DDBJ accessions through ENA:

`msc-convert DRP123456 --in-type ena --out-type magetab --out magetabs`

## Main options

These correspond to the main arguments of `Converter.convert()`.

| CLI argument | Purpose | Python argument |
| --- | --- | --- |
| Inputs after `msc-convert` | One or more accessions, supported files, URLs or directories. | `input` |
| `--in-type TYPE` | Select the reader; defaults to automatic detection. Alias: `--force-in-type`. | `in_type` |
| `--out-type TYPE` | Required output format: `json`, `magetab`, `tsv`, `csv`, `h5ad` or `obs`. | `out_type` |
| `--out DIRECTORY` | Output directory; the CLI defaults to the current directory. Aliases: `--outdir`, `-o`. | `outdir` |
| `--enrichment MODE` | Enrichment preset: `standard`, `curators` or `off`. Defaults to `standard`. | `enrichment` |

Enrichment presets:

- `standard`: Normal preparation and applicable linked-repository enrichment.
- `curators`: Native/publication preparation without enrichment from other study repositories.
- `off`: Disable supplementary enrichment.

Verified related-study expansion is enabled by default, independently of enrichment.
Use `--no-expand-studies` to disable it. The guide's curator audience does not change
the default `standard` enrichment preset.

## Forcing a platform

List available handlers:

`msc-convert --list-platform-handlers`

Force a platform while allowing automatic input detection:

`msc-convert GSE12345 --out-type magetab --platform-handler tenx_v3_droplet_single_cell_sequencing --out magetabs`

Available handlers:

- `plate_single_cell_sequencing`
- `droplet_single_cell_sequencing`
- `tenx_v2_droplet_single_cell_sequencing`
- `tenx_v3_droplet_single_cell_sequencing`
- `single_cell_sequencing`
- `spatial_sequencing`
- `bulk_sequencing`
- `sequencing`
- `array`
- `generic`

Choose the most specific handler matching the study. The list order does not indicate
specificity. Omitting `--platform-handler` enables automatic platform detection.

## Additional options

The [complete CLI reference](codebase.md#unified-cli) covers study expansion,
platforms, output settings, expression processing, runtime controls, logs and reports.
Most conversion settings have individual flags and map into the Python
`Converter.convert(options={...})` mapping; logging and result reports are CLI controls.
See the [Python options reference](codebase.md#unified-routes) for keys and applicability.

The terminal summary lists each input's status, output paths and diagnostic codes.
`--report-json PATH` saves a machine-readable result report without scientific payloads;
`--report-json -` emits JSON-only stdout. `-v` enables informational logs on stderr.

[Documentation index](index.md)
