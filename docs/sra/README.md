<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# SRA converter reference snapshot

This directory is the source-faithful, offline provider reference bundle used
by the implemented `sra2json` converter and its tests. The directory itself is
reference evidence, not runtime code or a promise that every provider field can
be mapped into scientifically complete MAGE-TAB. Provider artifacts are stored
byte-for-byte; mapping and limitations are documented in
[expected-fields.md](expected-fields.md).

The snapshot was taken on 2026-09-01. [manifest.json](manifest.json) records each source URL, retrieval time, content type, available HTTP validators, SHA-256 digest, local path, and fixture accession. A JSON `null` ETag or Last-Modified value means the response did not supply a value retained by the snapshot process.

## Contents

| Path | Contents | Intended use |
|---|---|---|
| `schemas/SRA.*.xsd` | Complete INSDC SRA 1.5 common, study, sample, experiment, run, analysis, submission, and receipt schemas from ENA | Element, attribute, occurrence, and controlled-enumeration contracts |
| `schemas/ncbi/` | NCBI BioSample and BioProject submission schemas | Linked-record structure and converter dependency reference |
| `catalogues/einfo-sra.xml` | NCBI E-utilities `EInfo` response for the SRA database | Searchable SRA field names and term metadata |
| `catalogues/biosample-packages.xml` | NCBI BioSample package definitions | Package-specific required and optional attributes |
| `catalogues/biosample-attributes.xml` | NCBI recognized BioSample attributes | Attribute names, harmonized names, descriptions, and formats |
| `fixtures/SRX017289/` | GEO- and PubMed-linked C. elegans chain | Study/sample/experiment/run, BioProject, BioSample, and PubMed example |
| `fixtures/SRX7812918/` | Non-GEO human-gut amplicon chain | Paired reads, geographic attributes, and deliberately absent publication example |

## Snapshotted NCBI endpoints

| Service | Request represented here | Response role |
|---|---|---|
| E-utilities EFetch | `efetch.fcgi?db=sra&id=<SRX>&retmode=xml` | NCBI composite `EXPERIMENT_PACKAGE_SET` containing experiment, study, sample, and runs |
| E-utilities EFetch | `efetch.fcgi?db=biosample&id=<SAMN>&retmode=xml` | BioSample identifiers, owner/status, package, organism, and free-form or recognized attributes |
| E-utilities EFetch | `efetch.fcgi?db=bioproject&id=<PRJNA>&retmode=xml` | BioProject description, scope, target, publications, and linked identifiers |
| E-utilities ESummary | `esummary.fcgi?db=pubmed&id=<PMID>` | Publication title, authors, DOI, source, dates, and publication status |
| E-utilities EInfo | `einfo.fcgi?db=sra` | SRA database description and query-field catalogue |
| BioSample documentation | `/biosample/docs/packages/?format=xml` | Machine-readable package definitions |
| BioSample documentation | `/biosample/docs/attributes/?format=xml` | Machine-readable recognized-attribute catalogue |

These are public, unrestricted records. No credentials, API keys, or protected dbGaP material are present.

## Fixture chains

| Fixture | Linked chain | Cross-record evidence |
|---|---|---|
| `SRX017289` | `SRP002056` → `SRS011830` → `SRX017289` → `SRR037073` | `PRJNA123835`, `SAMN00009557`, GEO `GSE18729`/`GSM465245`, PMID `20133686`; RNA-Seq, transcriptomic, size fractionation, single-end, Illumina Genome Analyzer II |
| `SRX7812918` | `SRP250911` → `SRS6225446` → `SRX7812918` → `SRR11192680` | `PRJNA609050`, `SAMN14218700`; amplicon, genomic, PCR, paired-end, MiSeq; BioSample collection date, `Canada: London`, and `43.01 N 81.27 W`; no PubMed response because no PMID is linked |

The fixture directory intentionally stores separate BioSample and BioProject responses even though the composite SRA XML also carries their accessions. `sra2json` joins these records by accession; physical nesting in an API response is not ownership.

## Schema use and validation limits

The `SRA.*.xsd` files are official INSDC submission schemas, while NCBI EFetch returns an archive-oriented composite. NCBI's `EXPERIMENT_PACKAGE_SET` contains injected structures such as `Organization`, `Pool`, `SRAFiles`, `CloudFiles`, run `Statistics`, and `Bases`, and can carry historical schema-location declarations. Therefore, well-formed XML parsing and focused field-contract tests are appropriate for these fixtures, but strict validation of the complete NCBI composite against only the SRA 1.5 submission XSDs is not expected to succeed. Do not weaken or rewrite the provider XML to force validation.

## Integrity and refresh policy

Provider files must not be reformatted, normalized, or given project headers. Refreshes should replace a provider artifact only with bytes retrieved from its manifest URL, then update its digest and retrieval metadata. Project-authored Markdown remains subject to the repository author-header policy. The offline tests in `tests/test_provider_reference_material.py` verify coverage, XML syntax, fixture linkage, and every recorded hash.
