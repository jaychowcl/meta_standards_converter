<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# Expected SRA and linked-record fields

This inventory documents the provider contract consumed by `sra2json`; it is not a parallel JSON schema. Cardinalities describe the provider model: `1`, `0..1`, `0..*`, and `1..*`. Archive-generated fields may be present in retrieval responses even when they are absent from submission XSDs.

## Object hierarchy and core contracts

| Object | XPath or field | Cardinality | Primitive/value contract | Field class |
|---|---|---:|---|---|
| package | `/EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE` | `1..*` | Composite retrieval unit; do not treat it as a biological entity | archive-injected |
| study | `STUDY/@accession`, `IDENTIFIERS/PRIMARY_ID` | `1` public record | `(S|E|D)RP` accession; BioProject may also appear as `EXTERNAL_ID[@namespace='BioProject']` | archive-injected identifier |
| study | `DESCRIPTOR/STUDY_TITLE`, `STUDY_ABSTRACT`, `STUDY_TYPE/@existing_study_type` | `0..1` each | Unicode free text except XSD-controlled study type | optional/free-form and controlled |
| study | `STUDY_LINKS/STUDY_LINK/*` | `0..*` | URL or database/ID cross-reference, including GEO and PubMed | optional cross-reference |
| sample | `SAMPLE/@accession`, `IDENTIFIERS/PRIMARY_ID` | `1` public record | `(S|E|D)RS` accession | archive-injected identifier |
| sample | `IDENTIFIERS/EXTERNAL_ID` | `0..*` | Namespace-qualified identifier; commonly BioSample `SAMN...` and GEO `GSM...` | optional cross-reference |
| sample | `TITLE`, `SAMPLE_NAME/TAXON_ID`, `SCIENTIFIC_NAME`, `COMMON_NAME` | title `0..1`; taxon normally `1` | Free text plus integer NCBI Taxonomy identifier | required/optional by submission context |
| sample | `SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/{TAG,VALUE,UNITS}` | `0..*` | Repeating name/value/unit tuples; package/catalogue supplies stronger contracts where recognized | checklist-driven or free-form |
| experiment | `EXPERIMENT/@accession`, `IDENTIFIERS/PRIMARY_ID` | `1` public record | `(S|E|D)RX` accession | archive-injected identifier |
| experiment | `STUDY_REF`, `DESIGN/SAMPLE_DESCRIPTOR` | `1` each in normal read experiments | References to study and sample records | required reference |
| library | `DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_{STRATEGY,SOURCE,SELECTION}` | `1` each in XSD model | XSD enumerations; exact spellings and case are provider values | required/controlled |
| library | `DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/*` | `1` | Choice of `SINGLE` or `PAIRED`; paired records can carry nominal length/deviation attributes | required/controlled |
| experiment | `PLATFORM/*/INSTRUMENT_MODEL` | normally `1` | Platform-specific XSD enumeration or provider text | required/controlled |
| run | `RUN_SET/RUN/@accession` | `1..*` per package | `(S|E|D)RR` accession | archive-injected identifier |
| run | `RUN/@alias`, `EXPERIMENT_REF` | alias `0..1`; ref `1` | Submitter name and experiment reference | optional/required |
| run | `SRAFiles/SRAFile` and `Alternatives` | `0..*` | Filename, URL, MD5, byte size, semantic/super type, provider/access mode | archive-injected |
| run | `Statistics/Read/@average` | `0..*` | Decimal/string read-length statistic per read index | calculated/archive-injected |
| BioSample | `/BioSampleSet/BioSample/@accession` | `1` | `SAMN...` identifier | archive-injected identifier |
| BioSample | `Description/Organism`, `Attributes/Attribute` | organism `1`; attrs `0..*` | Taxonomy plus attribute name, harmonized/display names, optional unit, and text value | package-driven/free-form |
| BioProject | `/RecordSet/DocumentSummary/Project` descendants | `1` | Project accession/ID, title, description, data type/scope, target, links and publications | mixed required/optional |
| PubMed | `/eSummaryResult/DocSum/Item[@Name=...]` | item-dependent | PMID, title, repeating authors, DOI, journal/source, dates, status | archive-injected/calculated summary |

The exhaustive primitive, enumeration, and occurrence rules live in `schemas/` and `schemas/ncbi/`. BioSample requirements must be resolved through both `biosample-packages.xml` and `biosample-attributes.xml`: the same attribute can be required for one package, optional for another, or accepted as submitter-defined free text.

## Accession and cross-reference behavior

Do not infer hierarchy only from prefixes. Resolve explicit `STUDY_REF`, `SAMPLE_DESCRIPTOR`, `EXPERIMENT_REF`, and namespace-qualified identifiers. Preserve primary and secondary accessions rather than replacing one with another. NCBI, ENA, and DDBJ prefixes (`S`, `E`, and `D`) identify the original INSDC broker, not a different semantic object type.

Publication is optional. `SRX017289` demonstrates a PubMed link on the study and a separate ESummary response. `SRX7812918` deliberately has no PubMed fixture; absence is represented by the missing file and must not be converted into an empty publication record or a failed required join.

## Fields MSC consumes today

The study-scoped `sra2json` converter resolves projects, studies, samples,
experiments, and runs through SRA ESearch history and batched EFetch. It joins
unique BioSample, BioProject, and PubMed records, then emits one typed MSC
MINiML 3.0 package per resolved study. The mapping consumes:

| MSC value | Provider source | Use |
|---|---|---|
| series identity and biology | study accession, identifiers, descriptor, type, links, contacts and organization | `series`, accession databases, contributors and publication discovery |
| sample identity and biology | sample identifiers/name/attributes plus linked BioSample | sample/channel taxonomy, descriptions, and ordered duplicate-preserving characteristics |
| explicit design | experiment design descriptions | typed protocols without inferring biology from titles |
| library and platform | library descriptor and platform instrument model | each `sra_run`; promoted to sample fields only under unanimous experiment agreement |
| assay/run hierarchy | study/sample/experiment/run references, aliases and statistics | structural assay paths, scan names, read statistics and accessions |
| files | run `SRAFiles/SRAFile` and alternatives | ordered file metadata with URI, size, checksum and provider semantics |
| publication | study PMID followed by PubMed ESummary | IDF publication identifiers, descriptions, DOI and authors when available |
| audit evidence | every consumed response | `source.documents` URI/media type/SHA-256 and `extensions.insdc` provenance, conflicts, warnings and unmapped values |

The older GEO `MINiMLEnricher` path remains intentionally compact and keeps its
existing subset: GEO-linked SRA accessions, experiment/run/library/instrument
fields, read lengths, and ENA file reports. `sra2json` does not change that
behavior.

## Precedence and optional origin enrichment

Base `sra2json` is source-faithful and does not harmonize or infer fields. Composite SRA XML is primary; linked records fill gaps, while conflicting alternatives remain in `extensions.insdc.conflicts`. Literal missing-value terms and duplicate attributes remain values rather than becoming null or being deduplicated.

`--enrich-geo` and `--enrich-ae` are mutually exclusive and disabled by default. When a compatible broker link is present but enrichment is off, the converter warns visibly and records the warning. Explicit enrichment requires a compatible link and successful retrieval. Origin values win only for submitted study/sample biology, factors, protocols, contacts, and publication descriptions after an explicit one-to-one accession match. INSDC identifiers, status, library/run/file metadata and provider provenance remain additive. Ambiguous or unmatched samples are retained with diagnostics rather than guessed.

## MAGE-TAB coverage and limitations

The emitted package is structurally valid input to `json2ae --no-enrich`: study descriptions/publications become IDF content; sample organism/attributes become SDRF sources and characteristics; explicit experiment descriptions become protocols; library/platform/run/files become assays and data-file columns. Provider metadata does not guarantee experimental factors, ontology identifiers, detailed wet-lab protocols, publications, or even a sequencing hierarchy. The output is therefore a source-faithful representation and can yield structurally valid but scientifically incomplete MAGE-TAB. A valid zero-run study is retained as metadata-only MINiML with a degradation warning.
