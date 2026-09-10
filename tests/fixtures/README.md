<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# Converter contract fixtures

These are offline, independently invoked input-to-artifact examples for all seven MSC converters. Run `.venv/bin/python -m pytest tests/e2e -q` from the repository root. Scientific dependencies must be installed in the project environment; the required expression examples must not silently skip. Installed MSC metadata must match the checkout (MSC 6.0.0 for this reviewed corpus).

## Corpus and provenance

| Case | Origin and reduction | Workflows |
| --- | --- | --- |
| `studies/GSE328265` | Existing public GEO XML reduced to GSM9651991 and its sample reference; other study/platform evidence retained. PubMed, SRA and ENA responses captured on 2026-09-10. | GEO2JSON, GEO2AE, JSON2AE, JSON2TSV |
| `studies/E-MTAB-6486` | Unmodified public IDF/SDRF captured on 2026-09-10: twelve rows, six samples, paired FASTQs. | AE2JSON |
| `studies/pbmc3k` | Public PBMC3k raw counts: first four cells and first six source-order genes expressed in those cells. No numerical transformation; transpose to a readable TSV. Constructed MINiML wrapper, clearly separate from the GEO study. | JSON2H5AD, JSON2OBS |
| `edge_cases/magetab-normalization` | Existing constructed one-row material-type contract, not an original E-MTAB-6486 download. | MAGE-TAB semantic regression |
| `edge_cases/geo-relations` | Constructed reciprocal GSE1/GSE2 relationship; reuses a PubMed response solely as test evidence. | Traversal, inheritance, absent evidence, provider failure |
| `edge_cases/atlas-groups` | Existing synthetic Atlas envelope with explicit native MINiML 3.0 metadata. | Grouped CSV output and failed-dataset diagnostics |

Each study has `provenance.json`. `contracts.sha256.json` seals every input, provider response, provenance record, and expected output in this corpus. The original GEO capture date is unknown and is not invented. The full PBMC download is disposable development material; its checksum, source shape, axis selection, and selected counts remain in tracked provenance. Normal tests download nothing.

## Expected outputs and review

`expected/<converter>/` contains actual filenames and complete expected text/JSON outputs. An expected `.h5ad.json` is a decoded semantic representation of the H5AD: matrix values, sparse format, dtype, ordered axes, categorical metadata, retained MINiML, other metadata, and auxiliary matrix slots. Binary HDF5 bytes are not a portable oracle. CSV parsing intentionally produces float64 counts from the float32 source dataset; values remain exact.

Initial candidates were produced with MSC 6.0.0 and reviewed against independent evidence before acceptance:

- GEO identities, titles, protocol text, organism, molecule, relations and sample references were checked against XML; DOI/authors/title and run/library identifiers against captured PubMed/SRA XML.
- MAGE-TAB sample identity, source values, run references, protocol definitions/references and ordering were checked against source evidence and documented construction contracts. Use explicit `single_cell_sequencing` for the GEO-derived MAGE-TAB examples; see the known auto-detection defect below.
- Public AE source rows were checked for six first-seen sample groups, twelve retained FASTQ URIs, organism and library layout. All five semantic warning messages are stored separately in `diagnostics.json` and asserted.
- PBMC matrix values were checked against the selected raw source slice and readable TSV. The four-by-six values are also explicitly asserted independently in `test_fixture_evidence.py`; expected identities follow MSC's documented sample suffix rule. MINiML retention, absence of expression integration, artifact paths and diagnostics are checked in complete snapshots and focused tests.
- The synthetic Atlas expectation contains exactly the completed GSE100 sample; its manifest retains the GSE200 collection failure diagnostic.

Expected outputs are never regenerated during tests. There is no automatic snapshot-accept command. A future update must explain the changed contract, compare source evidence, inspect full diffs, and refresh the checksum inventory only after review. Historical output alone is insufficient justification.

## Allowed execution differences

Only the temporary workspace root and JSON2OBS staging-directory token are replaced. MAGE-TAB submission date is first required to equal the actual execution date, then represented as `<RUN_DATE>`. Provider/content dates remain exact. Machine available memory is fixed at the external measurement boundary; MSC's estimator and admission algorithm still execute. Release provenance remains exact, not normalized away.

No scientific values, dtypes, ordering, identities, missing values or semantic diagnostics are dropped. Corruption tests change counts, axes, identities, dtypes, retained metadata and diagnostics and require the same comparison helper to reject them.

## Known defect MSC-TEST-001

The retained GSM9651991 extraction protocol explicitly describes Chromium Single Cell **5′** chemistry. Automatic technology detection currently renders **10xV2 / 3 prime tag**. This discrepancy predates the test work. `tests/e2e/test_known_defects.py` records the contradiction as a strict expected failure. The generic single-cell handler is explicitly selected in approved MAGE-TAB goldens to avoid adopting that contradiction as intended output. Production code is unchanged; a future scientific fix must resolve the assertion and remove the expected-failure marker with source-backed review.

## Test boundaries

HTTP replay replaces only `RateLimitedRequester.get`, validates requested resources and fails on unexpected requests. Real archive extraction, parsers, enrichment, converters and writers execute. Raw processing uses an injected bounded command runner that emits a valid tiny count matrix; it verifies MSC dispatch and downstream reading, not nf-core's biological computation. Existing focused tests retain overwrite/refusal, failure recovery, memory limits, partial results, protocol variants and model/patch contracts.
