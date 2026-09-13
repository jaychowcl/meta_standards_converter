<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# Expected SRA and linked-record fields

This is a converter-design inventory, not a final JSON schema. Cardinalities describe the provider model: `1`, `0..1`, `0..*`, and `1..*`. Archive-generated fields may be present in retrieval responses even when they are absent from submission XSDs.

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

## Compact GEO enrichment subset

The separate compact GEO enrichment route has a narrower contract than native `sra2json`. During GEO conversion, `MINiMLEnricher` extracts SRA accessions only from GEO sample relations of type `SRA`, then calls NCBI EFetch. It consumes the following subset:

| Current MSC value | Provider source | Use |
|---|---|---|
| `study` | `STUDY/@accession`, primary ID, `STUDY_REF`, then study external ID fallback | Stored in each `sample.sra_run`; also collected into the historically named `sample.ena_accession` list |
| `experiment`, `run`, `sample` | object accession or primary ID | `Comment[ENA_EXPERIMENT]`, `Comment[ENA_RUN]`, and `Comment[ENA_SAMPLE]` in sequencing SDRF |
| `biosample`, `geo_sample` | sample `EXTERNAL_ID` namespaces | BioSample export/provenance; `geo_sample` is comparison evidence only |
| `library_layout`, `library_selection`, `library_source`, `library_strategy` | experiment library descriptor | Sequencing library comments and protocol text |
| `instrument_model` | platform-specific `INSTRUMENT_MODEL` | Sequencing assay comment |
| `scan_name` | run alias, otherwise run accession | SDRF scan name |
| `fastq_files` | ENA file report first, otherwise NCBI `SRAFile` records | Read filenames, URI, MD5, and file selection |
| `read_lengths` | run `Statistics/Read/@average` | Preserved in enriched JSON; not a core GEO characteristic |
| PubMed publication tuple | GEO `series.pubmed_id` followed by PubMed ESummary | IDF publication ID, DOI, authors, title and harmonized publication status |

That compact route does **not** fetch the separate BioSample or BioProject responses, parse sample attributes such as `geo_loc_name`/`lat_lon`, consume study abstract/type, or discover PubMed IDs from the SRA study link. Native `sra2json` retrieves and maps these records.

## Compact GEO precedence and overwrite behavior

SRA/ENA enrichment does not overwrite GEO biological, geographic, or characteristic fields. It writes only `sample.sra_accession`, `sample.sra_run`, and `sample.ena_accession`; invoking enrichment on a package already containing those three enrichment slots replaces those slots, but it does not rewrite the GEO sample object.

For MAGE-TAB rendering, existing GEO sample-level library fields and instrument model explicitly win over differing SRA run values, with an audit warning. A differing SRA `geo_sample` accession also produces a warning and the GEO sample accession remains the assay name. ENA can replace the NCBI-derived FASTQ list for the same run when ENA supplies a non-empty file report; that is file provenance precedence, not GEO metadata precedence.

## Future `sra2json`/MAGE-TAB coverage

Together, the SRA composite, BioSample, BioProject, and optional PubMed response contain enough information to construct a useful sequencing MAGE-TAB package: study description/publications for IDF; sample organism and attributes for SDRF source characteristics; library/platform information for protocols and assays; and run/file information for scans and data files. They do not guarantee complete experimental protocols, ontology identifiers, factors, or publication links. A converter must preserve free-form attributes and provenance, apply package/checklist rules, model missing publications, and define conflict resolution before claiming a lossless mapping.

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
