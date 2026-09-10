# Structured chemistry and mixed technology evidence

`GSM5388031.xml` and `GSM9254695.xml` retain sample identity, complete channels
(including protocols and characteristics), description and library fields from
public GEO MINiML responses acquired on 2026-09-10. `manifest.json` records each
request URL, full response checksum, reduced fixture checksum and reduction.
The original full responses are acquisition evidence in the ignored development
workspace; the fixtures themselves are sufficient for every offline test.

The Series envelopes are constructed GSE900011/GSE900012 studies. They are not
biological records. Removed provider relations, processing and files keep these
checks focused on metadata interpretation without additional provider calls.
`mixed.xml` is entirely constructed, including its explicit index cycle recipe.

## Reviewed expected outputs

Complete IDF/SDRF expectations and operation audits are checked through both
GEO2AE and JSON2AE, followed by real MAGE-TAB reparsing. No fields are normalized
in the artifact comparisons. Protocol whitespace normalization is the existing
IDF registration contract, and repeated comments (including empty cells) must
survive reparsing in order.

- GSM5388031: the channel's `singlecell_type=SC3Pv2` identifies 10x 3-prime v2.
  The protocol's general v2-or-v3 statement cannot remove the sample assignment.
  The code supplies no read geometry or indexing configuration.
- GSM9254695: the sample title and library description identify scRNA sequencing.
  Its shared extraction protocol describes both snRNA and Visium preparation.
  The applicable Chromium 3-prime v3.1 preparation supports the scRNA library;
  a Visium mention cannot select the spatial handler or its geometry.
- Mixed: two constructed samples use one shared protocol registry and one SDRF
  column plan. The first has 3-prime v3.1 chemistry and explicitly reported
  index1/index2 lengths of 10. The second is explicitly Visium and retains the
  existing spatial rendering. The IDF uses the generic sequencing summary.

The real-record expectations were compared with isolated source at `6170907`.
GSM5388031 gains only end bias and library construction columns. GSM9254695
replaces Visium construction and geometry with source-supported single-cell
annotations and the existing single-cell IDF classification/empty comment rows.
All common SDRF columns preserve values and ordering; original protocols remain
unchanged. Existing earlier chemistry expectations were not replaced.

The identifier mapping is grounded in the [official Cell Ranger chemistry
options](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-multi-config-csv-opts),
reviewed 2026-09-10. Only the ten exact identifiers listed in the implementation
are supported here. Unsupported codes remain raw evidence with a diagnostic;
`auto` supplies no chemistry identity. Passing this corpus is not an estimate of
population-wide accuracy.
