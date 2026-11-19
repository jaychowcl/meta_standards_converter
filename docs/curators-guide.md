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

GEO exports retain source annotations and leave missing curator fields blank.
They do not infer organism part, factors or barcode recipes. Derived files remain
in JSON and are excluded from MAGE-TAB for every source. See the [source fidelity
contract](codebase.md#geo-source-fidelity), including guidance for older JSON.

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
| `standard` | Preserve populated GEO primary-input values; otherwise explicit ArrayExpress > GEO > SRA/ENA. | Retrieve applicable linked sources; retain alternatives and warn with the selected source and reason. |
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
New library conflicts use a stable field order within each run, so repeated
conversions retain the same ledger order. Imported evidence keeps its supplied order.

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

<a id="curation-warnings"></a>
## Finding metadata that needs curation

Use `-v` (or `--verbose`) to review every category of missing metadata after
enrichment. Each category gives the study, number of affected records, up to three
example locations and a suggested curator action. Use `-vv` for every location,
or `--report-json report.json` for the full structured findings:

Add `-v --report-json report.json` to one of the conversion commands above.

For a completely local assessment of supplied metadata, add
`--enrichment off --no-expand-studies`. Assessment itself never makes a lookup.
The same checklist runs for JSON, imported MAGE-TAB and other conversion routes
that have canonical metadata, regardless of output format or verbosity.

For example, a grouped warning for 18 samples looks like:

    WARNING [curation_missing] GSE123: organism part missing in 18 samples.
      Examples: GSM1 (package 0 /sample/0/channel/0/characteristics); GSM2 (package 0 /sample/1/channel/0/characteristics); GSM3 (package 0 /sample/2/channel/0/characteristics)
      Action: Supply verified organism part metadata or confirm it is not applicable.
      Use -vv or --report-json for all affected locations.

| Area | Checklist |
| --- | --- |
| Study | Title, description/summary, overall design, experimental type, contributors and contacts. A reference to a person without supplied information does not satisfy the check. |
| Publications | Citation presence; a usable PMID, DOI or PMCID; title, authors, DOI and status; completeness of supplied status ontology references. An absent citation produces one finding instead of a warning for every citation field. Sample, run and relation citations keep their own scope. |
| Samples | Title, source, organism, organism part (including the tissue alias), developmental stage, disease, genotype and molecular material; empty declared characteristics and factor values. Explicit sample-scoped factor declarations apply only to those samples. |
| Protocols | Descriptions and types of declared protocols; applicable collection, extraction, library preparation, sequencing, hybridization, scan and labeling procedures. Processed data also requires a processing description. Descriptions must be bound to the relevant sample or assay. |
| Assays and runs | Applicable library strategy/source/selection/layout, instrument, array design, run identifiers and raw file references. Explicit MAGE-TAB assay declarations and existing technology/chemistry resolution determine expectations; array assays do not require sequencing library fields. |
| Lookups | Existing retrieval failures plus successful lookups with no match or incomplete returned metadata. PubMed and run-enrichment outcomes retain provider and requested identifier. A missing field is linked to a lookup only when explicit identity evidence connects them. |

Absent values, `null`, whitespace-only text and empty collections count as
missing. Explicit values such as `unknown` and `not provided` are reported
separately as `curation_reported_unknown`. Zero, `false`, `none`, `N/A` and
`not applicable` remain legitimate supplied values. MSC does not infer a negative
biological annotation from missing data. Ontology/unit companions are checked
when a group is supplied; free text does not require ontology annotation.

**Curation findings are advisory.** `curation_*` warnings do not change execution,
completeness, validation or the exit code. A successful conversion with only these
findings still exits 0. Existing failed lookups, invalid metadata and conversion
errors retain their partial/failure behavior and exit codes. Successful lookups
with no match are distinguishable from failed requests; disabled enrichment is
recorded as policy, not a failed request. A field can still need curation even
when a request completed successfully. Scientific JSON, IDF and SDRF content is
unchanged by the assessment.

The ordinary console and `--quiet` suppress new advisories. Existing warning/error
visibility is unchanged. `--log-file run.log` or
`--log-file run.jsonl --log-format jsonl` retains findings independently of console
verbosity (subject to any explicit file-level filter). Reports and Python
`result.items[].diagnostics` always retain all findings. Each advisory includes
`curation.dataset_id`, `package_index`, `entity_type`, `entity_id`, `field` and
`reason`; `location` is a JSON pointer within the prepared package and `action`
describes the curator step. For an absent field, the pointer identifies its
expected field or containing metadata group. Available lookup context is additive.
`--report-json -` keeps stdout JSON-only; human output stays on stderr.

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

Each completed input gets a compact summary with usable output paths and up to
five diagnostic groups, followed by command totals. Each group states the problem
and a recovery action. Partial output written is distinct from no output produced.
Use `-v` for all groups; the result report and DEBUG log retain detailed evidence.
Related warnings from one provider failure share a correlation ID. For example,
a verified ENA fallback reports BioStudies access denial together with the ENA
substitution and the missing ArrayExpress-specific metadata.
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

For `missing_file`, check the named input and companions; for `permission_denied`
or `disk_full`, correct the destination before retrying. Provider retry eligibility
and cooldown are separate from warning severity: confirmed BioStudies access
denial is terminal for automatic retry, while transient optional failures may be
retried. For `conversion_failed`, retain the correlation ID and capture a DEBUG
JSONL log; it includes safe stack locations, never raw exception payloads.
If a log or report fails after conversion, the command returns nonzero and already
written conversion artifacts remain available. These controls also apply to
`miniml-migrate`, which retains its JSON summary on stdout.

## Reading conversion results

MSC checks strict MINiML compatibility before conversion and before publication.
An invalid source produces diagnostics and no output. Optional enrichment that
introduces invalid metadata is rejected; the last valid package can still be
published with partial completeness. Read the diagnostic operation and path to
identify what could not be accepted.

For ArrayExpress accessions, MSC checks BioStudies accessibility before discovering
IDF/SDRF files. Confirmed inaccessible/not-found studies and accessible studies
missing required files automatically try ENA. MSC requires one public ENA study
with an explicit link to the requested ArrayExpress accession; an alias or similar
title is insufficient. Ambiguous files, malformed responses and parsing errors
remain failures. Local IDF/SDRF inputs are never substituted automatically.

A verified fallback writes valid **partial** output because ArrayExpress-specific
annotations are unavailable. The report and `Comment[MSCSourceFallback]` identify
the ENA study, mapping evidence and original problem. JSON preserves that provenance
and remains partial on later export. This loading fallback also applies with
`--enrichment off --no-expand-studies`. Supply the study's original IDF and SDRF
files to recover ArrayExpress metadata. A confirmed access denial is not
automatically retried; temporary provider failures retain their own retry status.

The summary distinguishes `partial output written` from `no output produced`, with
the failure stage and suggested action. Exit zero means complete; partial and
failed conversions return nonzero, so inspect the report and output paths before
discarding results. See the [source-loading call graph](codebase.md#workflow-ae2json).

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
