<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# Expected ENA endpoint fields

This inventory describes information available for future `ena2json` and MAGE-TAB conversion. It is not a final JSON schema. Portal field catalogues are exhaustive for the snapshot; this document highlights structural and conversion-critical contracts.

## Record hierarchy and accession classes

| Level | Primary forms | Core links and values | Cardinality |
|---|---|---|---:|
| project/study | `PRJ(E|D|N)...`, `(E|D|S)RP...` | title, abstract/description, study type, center, publications, sample/experiment associations | one study to `0..*` experiments |
| sample/BioSample | `SAM(E|D|N)...`, `(E|D|S)RS...` | title/description, taxon, organism, checklist, repeating attributes and cross-references | one sample to `0..*` experiments |
| experiment | `(E|D|S)RX...` | one study reference, one sample descriptor, library descriptor, platform/instrument, design/protocol text | normally one study and one sample |
| run | `(E|D|S)RR...` | experiment reference, files, statistics, dates, center/alias | one experiment to `1..*` runs |
| analysis | `(E|D|S)RZ...` | analysis type, study/sample/run references, process and result files | `0..*` |
| assembly | `GCA_...`, `GCF_...` and ENA assembly forms | assembly name/type/level, genome representation, taxon, study/sample and files | `0..*` |
| taxon | integer tax ID | scientific name, rank, lineage, division, genetic codes and flags | one resolved taxon object |

The initial letter on SRA-style accessions indicates the submitting INSDC partner (`E` ENA, `D` DDBJ, `S` NCBI); it does not change the object semantics. Preserve both project/BioSample primary accessions and secondary study/sample accessions.

## Endpoint contracts

| Endpoint | Field/path contract | Primitive and cardinality | Field class |
|---|---|---|---|
| Portal `/results` | `resultId`, `description`, `primaryAccessionType`, `recordCount`, `lastUpdated` | one TSV row per result type; counts/dates are calculated snapshot values | calculated/archive-injected |
| Portal `/returnFields` | `columnId`, `description`, `type` | one TSV row per returnable indexed column | endpoint schema |
| Portal `/searchFields` | `columnId`, `description`, `type` | one TSV row per queryable indexed column | endpoint schema |
| Portal `/search` | requested/default column names | JSON array or headered TSV; zero rows is valid; values can be absent/blank | indexed, mixed required/optional |
| Portal `/controlledVocab` | `value`, `description` | JSON array; use exact case and spelling | controlled |
| Portal `/filereport` | accession plus file-family columns | JSON array or TSV; fields may be empty strings | archive-injected/calculated |
| Browser `/xml/{accession}` | SRA/ENA XML object roots and descendants | full XML, including repeating free-form attributes and links | required/optional/free-form |
| Xref `/json/search` | service-specific cross-reference objects | JSON array; empty array is a successful no-match response | archive-injected cross-reference |
| Taxonomy `/tax-id/{id}` | `taxId`, `scientificName`, `rank`, `lineage`, `division`, genetic codes, names/flags | one JSON object; most values are strings, `otherNames` is an array | controlled/archive-injected |

The vendored return/search field files cover `study`, `read_study`, `sample`, `read_experiment`, `read_run`, `analysis`, `analysis_study`, `assembly`, and `taxon`. Use the TSV `type` column (`text`, `number`, `date`, `boolean`, and provider variants) as the endpoint primitive contract; do not infer type from one fixture value.

## XML conversion-critical paths

| Object | XPath | Cardinality/value contract | MAGE-TAB relevance |
|---|---|---|---|
| study | `STUDY/@accession`, `IDENTIFIERS/*` | primary plus `0..*` external/submitter IDs | investigation accession and database links |
| study | `DESCRIPTOR/STUDY_TITLE`, `STUDY_ABSTRACT`, `STUDY_TYPE` | title/abstract text; controlled study type | IDF title, description, experiment type |
| study | `STUDY_LINKS/STUDY_LINK/{XREF_LINK,URL_LINK}` | `0..*` links | publication/project/GEO discovery |
| sample | `SAMPLE/@accession`, `SAMPLE_NAME/{TAXON_ID,SCIENTIFIC_NAME}` | accession and taxonomy | source name and organism |
| sample | `SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/{TAG,VALUE,UNITS}` | `0..*` repeating tuples | SDRF characteristics and factors with provenance |
| experiment | `STUDY_REF`, `DESIGN/SAMPLE_DESCRIPTOR` | normally one each | join keys |
| experiment | `LIBRARY_DESCRIPTOR/LIBRARY_*`, `LIBRARY_LAYOUT/*` | controlled strategy/source/selection and `SINGLE`/`PAIRED` | extraction/library protocol comments |
| experiment | `PLATFORM/*/INSTRUMENT_MODEL` | platform-specific model | assay hardware comment |
| run | `RUN/@accession`, `EXPERIMENT_REF` | one accession/reference | assay/run identity |
| run | file and statistics descendants | `0..*` files/read statistics | scan and raw-data files |
| checklist | `CHECKLIST/DESCRIPTOR/FIELD_GROUP/FIELD` | label/name/description, field type, mandatory state, multiplicity, units/restrictions | validation and typed characteristic contract |

Checklist `MANDATORY` values distinguish mandatory, recommended, and optional fields; field groups can impose additional restriction logic. `FIELD_TYPE` can be text, integer/decimal, date, ontology, taxonomy, or controlled-value oriented depending on the definition. Preserve the checklist accession and original attribute tag because aliases and field definitions change over time.

## Controlled values and missing values

The 18 Portal snapshots define the accepted values returned for analysis type, assembly level, broker, category, checklist, datahub, genome representation, host sex, instrument model/platform, library layout/selection/source/strategy, NCBI reporting standard, sex, tag, and taxonomic division. `library_layout.json`, for example, contains exactly `SINGLE` and `PAIRED`.

INSDC missing-value vocabulary is semantically different from an empty field. Top-level forms include `not applicable`, `missing: not collected`, `missing: not provided`, and `missing: restricted access`; current granular reporting forms include reasons such as `control sample`, `sample group`, `synthetic construct`, `lab stock`, `third party data`, `data agreement-established pre-2023`, `endangered species`, and `human-identifiable`. Preserve the complete literal and do not coerce it to JSON null. Optional fields with no value are normally omitted rather than filled with a missing term.

## Semicolon-aligned file reports

ENA represents multiple files within a run by semicolon-delimited sibling columns. Split `fastq_ftp`, `fastq_md5`, and `fastq_bytes` independently but align items by zero-based position. For `SRR11192680`, index 0 is read 1 across all three columns and index 1 is read 2. A shorter hash/byte column yields an unknown value for that file; it must not shift later values. Empty strings mean that file family is not supplied. Scheme-less `ftp.sra.ebi.ac.uk/...` values can be normalized to `ftp://...` while retaining the provider value in provenance.

## Cross-reference behavior

Portal primary/secondary accession columns and Browser XML references establish most joins. Xref results are supplementary and can be empty, as both fixtures demonstrate. Publication discovery must inspect study XML links and, when present, fetch a publication provider such as PubMed; ENA does not guarantee a PMID for every study. Taxonomy resolution enriches the sample tax ID but must not replace the submitted sample attribute values.

## Compact GEO enrichment subset

The compact GEO enrichment route calls ENA Portal `/filereport` with `result=read_run` and fields `run_accession,fastq_ftp,fastq_md5,fastq_bytes`. It groups non-empty files by run and prefers that ENA file list over NCBI `SRAFile` entries. Native `ena2json` additionally uses Portal search/catalogues, Browser XML, Xref, taxonomy, analysis and assembly retrieval. Checklist definitions remain reference material.

The rest of the current enrichment comes from NCBI SRA EFetch. Current sequencing SDRF rendering uses accessions, BioSample/GEO links, library values, instrument model, run alias, FASTQ URI/hash, and read length. PubMed ESummary is driven by PubMed IDs already parsed from GEO MINiML, not discovered through ENA.

## Compact GEO precedence and overwrite behavior

No ENA response writes GEO country, latitude/longitude, organism, source characteristics, title, or protocols. GEO sample-level library fields and instrument model win over conflicting SRA-derived values and produce audit warnings. A differing SRA `geo_sample` link also leaves the GEO accession in place. The sole active ENA precedence rule is a non-empty ENA FASTQ report replacing the NCBI-derived file list for the same run.

## Future `ena2json`/MAGE-TAB coverage

Browser XML plus Portal file reports, taxonomy, optional cross-references, and PubMed lookup can supply the structural content for sequencing MAGE-TAB: study/publication IDF rows, sample SDRF characteristics, experiment/library protocols, assay hardware, and run/files. Portal search alone is insufficient because only indexed fields are returned. Even full XML may omit publication, factors, ontology IDs, or detailed wet-lab protocols, so a future converter needs explicit absent-value handling, lossless extension storage, provenance, and conflict rules rather than assuming every MAGE-TAB field is derivable.

## Native converter field contract

The current [native field mapping](../codebase.md#native-archive-contract) and
[fidelity rules](../codebase.md#native-archive-fidelity) govern native imports.
They preserve full structured records in `extensions.insdc` 1.0 with MINiML 3.0.
BioProject accession retrieval uses exact PRJA-to-UID resolution and validated
ArchiveID identity, including legacy PRJDA identifiers. Explicit assembly
versions stay distinct. ENA indexed read/base counts stay in run-level
`indexed_statistics`; they are not inferred read lengths. Literal FTP filename
characters survive usable URI encoding.

Optional GEO/ArrayExpress enrichment preserves matched material/protocol paths,
scoped factors and units, native membership and additive file branches. Database,
contributor and organization references are reconciled with their declarations.
Files export their explicit links as `Comment[File URI]`. Operational diagnostics
belong only in logs and optional reports; absent metadata never justifies an
invented molecule, date, overall design or protocol step.
