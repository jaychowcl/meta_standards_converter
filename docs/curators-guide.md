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

Install MSC with Python 3.12 or newer:

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

<a id="enrichment-presets"></a>
Enrichment presets:

| Mode | Canonical selection | Conflicts and retrieval |
| --- | --- | --- |
| `standard` | Explicit ArrayExpress > GEO > SRA/ENA, independent of the entry route. | Retrieve applicable linked sources; retain alternatives and warn with the selected source and reason. |
| `curators` | Preserve the primary input source's existing populated value. | Native/publication preparation; inspect supplied conflict evidence without fetching other study repositories for comparison. |
| `off` | Preserve supplied metadata. | No supplementary reconciliation. |

Native means the primary input source: an AE input is AE-native. Between SRA
and ENA, the primary provider keeps its populated values and its peer fills gaps.
Only statements joined to the same sample, experiment/library, run and field are
compared. Distinct explicit biological characteristic values remain multivalued;
repository ranking does not replace library-scoped modality interpretation.

Scientific disagreement produces warnings and **partial completeness with valid
output**. Conflict evidence retains source accession, field path, original values,
selected value and reason. This is separate from invalid metadata: an invalid
optional candidate is rejected and the last valid package is retained. `curators`
does not check repositories it has not retrieved; absence of a warning does not
establish agreement with those repositories.

Conflict evidence travels with MAGE-TAB in versioned `Comment[MSCConflict]`
entries and a readable summary. TSV/CSV and AnnData observations include the
applicable ledger in `msc.provenance.conflicts`; AnnData also retains the package
ledger. Reimporting these conflicts keeps the warnings, even with enrichment off,
and does not select values again. Unknown comment versions remain recoverable.

Detailed modality combines compatible library and bound preparation evidence.
Broad single-cell metadata can be refined by nuclei isolation; contradictory
explicit preparations remain unknown. Generated acquisition paths are reconciled
only when verified identities establish that a retained continuation covers their
evidence. Authored archive row occurrences and distinct processing branches remain.

One run can produce several SDRF rows for submitted, archive-generated or processed
files. Different file locations, checksums and sizes remain evidence; similar
basenames or compressed/uncompressed names do not prove equivalence. A sample-only
processed-file branch may have no assay factors. Compare supported identities and
relationships before interpreting extra rows as duplicate samples or missing data.

MAGE-TAB comments stay with their original node, field, protocol or unit occurrence.
Factor values describe the complete experimental row and are rendered after the
workflow columns, together with their qualifiers. Ambiguous duplicate scalar fields
remain explicit residual evidence with diagnostics rather than replacing a value.

Possible source-text encoding corruption is retained and flagged. Standard
enrichment may use a clean equivalent only when the same scoped field corroborates
a reversible decoding transformation; both statements and normalization evidence
remain available. Curators preserves the primary source's text.

Python calls can optionally bound metadata retrieval, for example
`options={"metadata_timeout_seconds": 900, "metadata_max_requests": 1000}`.
Both limits default to `None`. Each top-level input shares its allowance across
native, peer, linked and publication retrieval; retries and redirects count as
HTTP attempts, cache hits do not. On exhaustion, MSC returns the latest strictly
valid metadata with partial coverage diagnostics. If no valid candidate exists,
it fails without publication. These limits do not cover subsequent expression
processing and are currently Python runtime options rather than dedicated CLI flags.

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

The terminal summary lists status totals and output paths; the result report retains detailed diagnostics.
`--report-json PATH` saves a machine-readable result report without scientific payloads;
`--report-json -` emits JSON-only stdout. Human summaries, logs and progress always use stderr.

Interactive terminals show an automatic progress bar. Use `--progress plain` for
stable stage markers, or `--progress off` to disable it. Unknown totals show activity
rather than an estimated percentage. `-v` adds stage statistics; `-vv` adds diagnostic
detail. `--quiet` suppresses progress and human summaries but keeps errors and JSON results.

To keep complete logs with a quiet terminal:

`msc-convert GSE100 --out-type json --out results --quiet --log-file run.jsonl --log-format jsonl`

Log files capture DEBUG by default independently of the console. Use
`--log-file-level INFO` only when a reduced file is desired. The default file format
is readable text; JSONL is suitable for automated inspection. See the
[logging and progress reference](codebase.md#conversion-logging).

## Reading conversion results

MSC checks strict MINiML compatibility before conversion and before publication.
An invalid source produces diagnostics and no output. Optional enrichment that
introduces invalid metadata is rejected; the last valid package can still be
published with partial completeness. Read the diagnostic operation and path to
identify what could not be accepted.

Default MAGE-TAB preserves deposited archive metadata. Array-design identifiers
keep their actual namespace; protocol, factor and vocabulary references are
checked before writing. This does not guarantee an Atlas processing layout.
Release/update labels follow the source repository. The conversion date and
submission-date fallback for Date of Experiment retain their existing behavior.

Tables and AnnData add `msc.expression.modality_detail`: bulk, single_cell,
single_nucleus, spatial, mixed or unknown. The legacy modality field is unchanged.
Unknown means the available library evidence is insufficient or conflicting;
forcing a renderer does not change biological evidence. Canonical tissue uses
explicit tissue/organism-part characteristics, not a source description such as
total RNA. Distinct explicit values remain visible for review.

See [fidelity contracts](codebase.md#fidelity-contracts) for qualification, evidence
scope, rollback and the recorded multi-technology acceptance panel.

[Documentation index](index.md)
