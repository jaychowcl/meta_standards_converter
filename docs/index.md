<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# Documentation index

Use the [README](../README.md) for installation and first use. For onboarding,
read [Architecture](codebase.md#architecture), [OOP design](codebase.md#oop-design),
then one [workflow](codebase.md#principal-workflows). This index routes into the
canonical handoff by stable anchor, purpose and search keywords.

## Canonical routes

- id: architecture
  title: Architecture
  anchor: architecture
  purpose: Start here for MSC 8 purpose, subsystem ownership and the whole-system map.
  keywords: MSC 8, architecture, converters, source, model, output
  link: [Open section](codebase.md#architecture)

- id: system-context-and-boundaries
  title: System context and boundaries
  anchor: system-context-and-boundaries
  purpose: Maps actors, external services, trust/process boundaries, and failure ownership.
  keywords: boundary, GEO FTP, NCBI, ENA, PubMed, OLS, BioStudies, filesystem, Nextflow, Docker, Atlas v1
  link: [Open section](codebase.md#system-context-and-boundaries)

- id: architectural-decisions
  title: Architectural decisions
  anchor: architectural-decisions
  purpose: Records evidence-backed decisions and observed design choices with consequences.
  keywords: decision, rationale, evidence, CLI, MAGE-TAB, projectors, rootless
  link: [Open section](codebase.md#architectural-decisions)

- id: design-invariants-and-expectations
  title: Design invariants and expectations
  anchor: design-invariants-and-expectations
  purpose: Defines compatibility, validation, lifecycle, security, overwrite, and failure rules.
  keywords: invariant, validation, compatibility, precedence, overwrite, logging, partial, rootless
  link: [Open section](codebase.md#design-invariants-and-expectations)

- id: proposed-enriched-miniml-core
  title: Enriched MINiML-compatible core
  anchor: proposed-enriched-miniml-core
  purpose: Defines MSC MINiML 3.0 protocols, assay paths, occurrence-local hz groups, units, provenance, document-scoped SDRF ordering, migration, and semantic MAGE-TAB boundaries.
  keywords: MSC MINiML 3.0, MAGE-TAB, hz fields, harmonized value, canonical IDF labels, units, factors, protocols, assay paths, v3-only boundary, semantic round trip
  link: [Open section](codebase.md#proposed-enriched-miniml-core)

- id: data-contracts
  title: Data contracts and preservation
  anchor: data-contracts
  purpose: Distinguish canonical metadata, export schemas, provenance and round-trip guarantees.
  keywords: MINiML 3, patch 3.1, Atlas 1, H5AD 2, assay 3, transport 1, hz, provenance
  link: [Open section](codebase.md#data-contracts)

- id: component-relationships-and-data-flow
  title: Component relationships and data flow
  anchor: component-relationships-and-data-flow
  purpose: Shows control/data exchange, ownership, lifecycle, and failure propagation.
  keywords: relationship, data flow, converter, fetcher, parser, projector, NFCoreRunner, Atlas
  link: [Open section](codebase.md#component-relationships-and-data-flow)

- id: oop-design
  title: Design and object relationships
  anchor: oop-design
  purpose: Understand class composition, inheritance, protocols, injected collaborators and state ownership.
  keywords: OOP, design, classes, inheritance, composition, lifecycle, dependency injection
  link: [Open section](codebase.md#oop-design)

- id: entrypoints-and-interfaces
  title: Entrypoints and interfaces
  anchor: entrypoints-and-interfaces
  purpose: Inventories supported CLI, Python, Docker, Compose, and operational interfaces.
  keywords: entrypoint, CLI, Python, Docker, Compose, geo2ae, geo2json, json2ae, ae2json, json2h5ad, json2tsv, json2obs
  link: [Open section](codebase.md#entrypoints-and-interfaces)

- id: orchestrators-and-core-types
  title: Orchestrators and core types
  anchor: orchestrators-and-core-types
  purpose: Explains state ownership and responsibilities of converters, constructors, planners, runners, and fetchers.
  keywords: orchestrator, AEConstructor, AEParser, JSON2H5ADConverter, DatasetCombinationPolicy, scientific compatibility, JSON2DelimitedConverter, RateLimitedRequester
  link: [Open section](codebase.md#orchestrators-and-core-types)

- id: public-api-reference
  title: Public API reference
  anchor: public-api-reference
  purpose: Read the major integration contracts before looking up exact source signatures.
  keywords: converter, model, projector, validation, constructor, failure, side effects
  link: [Open section](codebase.md#public-api-reference)

- id: atlas-v1-reader
  title: Atlas v1 reader facade
  anchor: atlas-v1-reader
  purpose: Defines the standalone versioned wire reader, dataset adaptation, warnings, v1 rejection, and dependency boundary.
  keywords: AtlasV1Reader, AtlasV1Error, schema_version, golden fixture, harmonized, status 2.0, metadata.packages, DatasetPackageGroup, standalone, no dependency
  link: [Open section](codebase.md#atlas-v1-reader)

- id: principal-workflows
  title: Principal workflows
  anchor: principal-workflows
  purpose: Choose a conversion or legacy import flow and trace validation through terminal outcomes.
  keywords: seven converters, eight commands, workflow, stages, failure, pseudocode
  link: [Open section](codebase.md#principal-workflows)

- id: workflow-geo2ae
  title: geo2ae: GEO to MAGE-TAB
  anchor: workflow-geo2ae
  purpose: Traces GEO retrieval, parsing, enrichment, construction, writing, and errors.
  keywords: geo2ae, GSE, GEO FTP, MINiML, AEConstructor, IDF, SDRF
  link: [Open section](codebase.md#workflow-geo2ae)

- id: workflow-geo2json
  title: geo2json: GEO to parsed JSON
  anchor: workflow-geo2json
  purpose: Traces GEO-to-JSON parsing, guarded direct-parent publication inheritance, optional enrichment, writing, and errors.
  keywords: geo2json, GSE, GEOParser, contributor order, PYTHONHASHSEED, GSE60450, parent publication, SubSeries, SuperSeries, reciprocal, PubMed, publication_inheritance, RelatedSeriesParseResult, related series, degraded, partial, safe error, enrichment, JSON
  link: [Open section](codebase.md#workflow-geo2json)

- id: workflow-json2ae
  title: json2ae: MINiML/Atlas v1 JSON to MAGE-TAB
  anchor: workflow-json2ae
  purpose: Traces native MINiML or canonical Atlas v1 loading, filtering, validation, enrichment, typed semantic MAGE-TAB reconstruction, and writing.
  keywords: json2ae, JSON, Atlas v1, harmonized, warning, MAGE-TAB, canonical IDF labels, publication status, protocol ontology, round trip, overlay, IDF, SDRF
  link: [Open section](codebase.md#workflow-json2ae)

- id: workflow-ae2json
  title: ae2json: MAGE-TAB to parsed JSON
  anchor: workflow-ae2json
  purpose: Traces bounded local, policy-approved HTTPS, and paginated BioStudies resolution through strict typed JSON output.
  keywords: ae2json, BioStudies, E-MTAB-1, duplicate term sources, unresolved references, E-MTAB-6486, pagination, IDF, SDRF, canonical labels, legacy label aliases, factor normalization, material type, resource profile, egress, AEParser, typed model
  link: [Open section](codebase.md#workflow-ae2json)

- id: workflow-json2h5ad
  title: json2h5ad: MINiML/Atlas v1 JSON to H5AD
  anchor: workflow-json2h5ad
  purpose: Traces grouping, asset planning, raw/processed paths, sample projection, catalogue-only publication, fail-closed compatibility evidence, and partial failures.
  keywords: json2h5ad, Atlas, MINiML, SourcePlanner, DatasetCombinationPolicy, per-sample catalogue, no expression integration, Entrez, gene symbol, all unknown evidence, memory report, admission, force memory, 70 percent, 90 percent, sequential durable checkpoint, scope index, asset performance, typed annotations, organism, reference, modality, feature namespace, deprecated allow unverified combination, msc_miniml, packages_json, ConversionResult, BatchConversionResult, nf-core, AnnData, lazy Scanpy, processed checkpoint, resume
  link: [Open section](codebase.md#workflow-json2h5ad)

- id: workflow-json2tsv
  title: json2tsv: MINiML/Atlas JSON to sample manifest
  anchor: workflow-json2tsv
  purpose: Traces grouping, projection, fail-closed validation, ordering, and TSV output.
  keywords: json2tsv, Atlas, projector, package_source, metadata_service, manifest injection, TabularConversionResult, allow_invalid
  link: [Open section](codebase.md#workflow-json2tsv)

- id: workflow-json2obs
  title: json2obs: MINiML/Atlas JSON to AnnData metadata
  anchor: workflow-json2obs
  purpose: Traces catalogue-backed observation-row aggregation without expression integration, single-sample var, and typed sample-namespaced uns publication.
  keywords: json2obs, obs CSV, row aggregation, no matrix combination, backed read, single-sample var CSV, sample-namespaced uns JSON, cell_id, feature_id, atomic, partial, memory admission, force memory, processed checkpoint, resume
  link: [Open section](codebase.md#workflow-json2obs)

- id: workflow-miniml-migrate
  title: miniml-migrate: explicit legacy source import
  anchor: workflow-miniml-migrate
  purpose: Import supported legacy source JSON directly as v3, including validation and overwrite behavior.
  keywords: miniml-migrate, MINiMLV1Migrator, unversioned, 1.0, v2 rejected, diagnostics
  link: [Open section](codebase.md#workflow-miniml-migrate)

- id: extension-and-change-guidance
  title: Extension and change guidance
  anchor: extension-and-change-guidance
  purpose: Defines supported extension patterns, compatibility concerns, and required tests.
  keywords: extension, projector, platform handler, asset, service, round trip, tests
  link: [Open section](codebase.md#extension-and-change-guidance)

- id: insdc-provider-reference-material
  title: INSDC provider reference material
  anchor: insdc-provider-reference-material
  purpose: Routes the vendored SRA/ENA schemas, endpoint contracts, catalogues, fixture chains, integrity manifests, MSC field-consumption boundary, and checklist availability evidence for future converter design.
  keywords: SRA, ENA, INSDC, XSD, OpenAPI, EInfo, BioSample, BioProject, PubMed, Portal, Browser, Xref, taxonomy, checklist, expected fields, manifest, sra2json, ena2json, MAGE-TAB
  link: [Open section](codebase.md#insdc-provider-reference-material)

- id: project-purpose-and-layout
  title: Project Purpose And Layout
  anchor: project-purpose-and-layout
  purpose: Explain project purpose and layout and its implementation boundaries.
  keywords: project purpose, layout
  link: [Open section](codebase.md#project-purpose-and-layout)

- id: configuration
  title: Configuration
  anchor: configuration
  purpose: Configure resources, handlers, assets, references, profiles, logging and runtime precedence.
  keywords: configuration, environment, platform handlers, manifest, params, Nextflow, replacement_profile
  link: [Open section](codebase.md#configuration)

- id: cli
  title: CLI reference
  anchor: cli
  purpose: Look up every argument, alias, default and mutually exclusive option for all eight commands.
  keywords: CLI, flags, argparse, geo2json, geo2ae, json2ae, ae2json, json2tsv, json2h5ad, json2obs, miniml-migrate
  link: [Open section](codebase.md#cli)

- id: python-api-guide
  title: Python API guide
  anchor: python-api-guide
  purpose: Use the converters through Python with concrete examples and result contracts.
  keywords: Python, convert, convert_source, export_manifest, typed packages, results
  link: [Open section](codebase.md#python-api-guide)

- id: docker-guide
  title: Docker guide
  anchor: docker-guide
  purpose: Build and run the image or provision the supported rootless Compose workflow.
  keywords: Docker, Compose, mount, daemon, nfcore-runner, installation
  link: [Open section](codebase.md#docker-guide)

- id: request-helper
  title: Host-aware request policy
  anchor: request-helper
  purpose: Defines conservative provider defaults, cross-process hostname scheduling and cooldowns, in-flight limits, bounded retries, and safe telemetry.
  keywords: RateLimitedRequester, HostRequestGate, HostRequestCooldownDeferred, slot, lease, model concurrency, RequestSettings, NCBIApplicationIdentity, tool, email, api_key, hostname, flock, Retry-After, max_inline_wait, jitter, provider_attempts, retry_count, rate_wait_seconds, 403, 429
  link: [Open section](codebase.md#request-helper)

- id: runtime-behavior
  title: Runtime Behavior
  anchor: runtime-behavior
  purpose: Find Python requirements, package dependencies, release evidence and runtime boundaries.
  keywords: MSC 8, installation, dependencies, packaging, Nextflow, Docker, release policy
  link: [Open section](codebase.md#runtime-behavior)

- id: end-to-end-geo2ae-flow
  title: GEO to MAGE-TAB
  anchor: end-to-end-geo2ae-flow
  purpose: Explain geo to mage-tab and its implementation boundaries.
  keywords: geo to mage-tab
  link: [Open section](codebase.md#end-to-end-geo2ae-flow)

- id: end-to-end-json2ae-flow
  title: JSON to MAGE-TAB
  anchor: end-to-end-json2ae-flow
  purpose: Explain json to mage-tab and its implementation boundaries.
  keywords: json to mage-tab
  link: [Open section](codebase.md#end-to-end-json2ae-flow)

- id: end-to-end-ae2json-flow
  title: End-To-End ae2json Flow
  anchor: end-to-end-ae2json-flow
  purpose: Explain end-to-end ae2json flow and its implementation boundaries.
  keywords: end-to-end ae2json flow
  link: [Open section](codebase.md#end-to-end-ae2json-flow)

- id: geo2json-vs-ae2json
  title: geo2json Versus ae2json
  anchor: geo2json-vs-ae2json
  purpose: Explain geo2json versus ae2json and its implementation boundaries.
  keywords: geo2json versus ae2json
  link: [Open section](codebase.md#geo2json-vs-ae2json)

- id: json2h5ad-flow
  title: End-To-End json2h5ad Flow
  anchor: json2h5ad-flow
  purpose: Explain end-to-end json2h5ad flow and its implementation boundaries.
  keywords: end-to-end json2h5ad flow
  link: [Open section](codebase.md#json2h5ad-flow)

- id: h5ad-metadata-schema-v1
  title: H5AD metadata schema 2.0
  anchor: h5ad-metadata-schema-v1
  purpose: Defines canonical dotted observation columns, normalized sample values, schema markers, source-column treatment, and globally unique observation identifiers.
  keywords: H5AD, schema 1.0, obs, msc_metadata, msc_assay, typed annotations, assay parameters, sample_values, canonical columns, multivalue, observation ID
  link: [Open section](codebase.md#h5ad-metadata-schema-v1)

- id: json2tabular-flow
  title: End-To-End json2tsv Manifest Flow
  anchor: json2tabular-flow
  purpose: Explain end-to-end json2tsv manifest flow and its implementation boundaries.
  keywords: end-to-end json2tsv manifest flow
  link: [Open section](codebase.md#json2tabular-flow)

- id: reference-annotation-flow
  title: Reference And Annotation Flow
  anchor: reference-annotation-flow
  purpose: Explains FASTA/GTF/GFF validation, conversion, precedence, and provenance.
  keywords: genome, FASTA, GTF, GFF3, gffread, annotation, checksum
  link: [Open section](codebase.md#reference-annotation-flow)

- id: rootless-json2h5ad-runtime
  title: Rootless json2h5ad Runtime
  anchor: rootless-json2h5ad-runtime
  purpose: Documents the hardened Docker/Nextflow process boundary and filesystem contract.
  keywords: rootless, Docker, Compose, socket, ACL, Nextflow, noexec
  link: [Open section](codebase.md#rootless-json2h5ad-runtime)

- id: rootless-acceptance-2026-07-31
  title: Rootless acceptance evidence — 2026-07-31
  anchor: rootless-acceptance-2026-07-31
  purpose: Routes to the successful pinned nf-core rootless acceptance summary and its portable report.
  keywords: rootless, acceptance, rnaseq 3.26.0, scrnaseq 4.2.0, H5AD, return code 0
  link: [Open section](codebase.md#rootless-acceptance-2026-07-31)

- id: parsed-miniml-data-shape
  title: Parsed MINiML data shape
  anchor: parsed-miniml-data-shape
  purpose: Defines the internal JSON package shape shared by converter workflows.
  keywords: MINiML, JSON, package, series, sample, platform, enrichment
  link: [Open section](codebase.md#parsed-miniml-data-shape)

- id: miniml-package-model
  title: MSC MINiML 3.0 package API
  anchor: miniml-package-model
  purpose: Understand immutable v3 entities, codecs, diagnostics, hz groups and retained source-bound patches.
  keywords: MINiMLPackage, MINiMLCodec, strict, hz, patch 3.1, source evidence, extensions
  link: [Open section](codebase.md#miniml-package-model)

- id: workflow-details
  title: Workflow Details
  anchor: workflow-details
  purpose: Explain workflow details and its implementation boundaries.
  keywords: workflow details
  link: [Open section](codebase.md#workflow-details)

- id: public-api-and-callable-reference
  title: Public API and callable reference
  anchor: public-api-and-callable-reference
  purpose: Look up every public-named production definition, constructor, method and property with source links.
  keywords: API, signature, callable, class, protocol, fields, methods, helpers
  link: [Open section](codebase.md#public-api-and-callable-reference)

- id: package-exports
  title: Owning-package exports
  anchor: package-exports
  purpose: Find explicit package exports and the defining modules behind lazy facades.
  keywords: __all__, owning package, converters, sources, expression, miniml, metadata, magetab
  link: [Open section](codebase.md#package-exports)

- id: maintenance-notes
  title: Maintenance Notes
  anchor: maintenance-notes
  purpose: Explain maintenance notes and its implementation boundaries.
  keywords: maintenance notes
  link: [Open section](codebase.md#maintenance-notes)

- id: harmonization-overrides
  title: Converter-owned replacement profiles
  anchor: harmonization-overrides
  purpose: Apply explicit converter replacement profiles while retaining source and harmonized evidence.
  keywords: replacement_profile, hz, fallback, invalid profile, copy, envelope rejection, provenance
  link: [Open section](codebase.md#harmonization-overrides)

- id: test-plan
  title: Test Plan
  anchor: test-plan
  purpose: Choose current offline validation commands and find the tests protecting each subsystem.
  keywords: pytest, offline, e2e, CLI, documentation, import, source, projection, integration
  link: [Open section](codebase.md#test-plan)

- id: live-api-provider-contracts
  title: Live API provider contracts
  anchor: live-api-provider-contracts
  purpose: Owns opt-in PubMed, SRA/ENA, and BioStudies provider checks.
  keywords: live_api, PubMed, SRA, ENA, BioStudies, pagination, no retries
  link: [Open section](codebase.md#live-api-provider-contracts)

- id: durable-artifact-publication
  title: Durable artifact publication
  anchor: durable-artifact-publication
  purpose: Defines immutable generations, the atomic current pointer, compatibility views, digest validation, and recovery failures.
  keywords: artifact, generation, current.json, fsync, rollback, recovery, bundle pointer
  link: [Open section](codebase.md#durable-artifact-publication)

- id: neutral-ae-construction-state
  title: Neutral AE construction state
  anchor: neutral-ae-construction-state
  purpose: Defines the one-way ProtocolRegistry and technology-detection dependency shared by IDF and SDRF construction.
  keywords: ae_common, ProtocolRegistry, technology detection, SDRFConstructor, AEConstructor, import cycle
  link: [Open section](codebase.md#neutral-ae-construction-state)

- id: operational-events-v1
  title: Operational events v1
  anchor: operational-events-v1
  purpose: Shared redacted JSONL/logging envelope and JSON/Prometheus metrics export at provider boundaries.
  keywords: OperationalEventEmitter, schema 1.0, JSONL, logging, redaction, metrics, Prometheus, RateLimitedRequester
  link: [Open section](codebase.md#operational-events-v1)

- id: runtime-contracts-v2
  title: Runtime contracts v2
  anchor: runtime-contracts-v2
  purpose: Defines independent status axes, safe error serialization, standard/large disk and RAM profiles, explicit overrides, memory admission, and disk preflight.
  keywords: status 2.0, execution, completeness, evidence confidence, validation, publication, SafeErrorEnvelope, resource profile, standard, large, 8 GiB, 32 GiB, 70 percent available RAM, 90 percent force memory, disk headroom
  link: [Open section](codebase.md#runtime-contracts-v2)

- id: secure-retrieval-and-xml
  title: Secure retrieval and XML boundaries
  anchor: secure-retrieval-and-xml
  purpose: Defines provider-host/public-address policy, redirects, byte quotas, cache integrity, safe multi-member archives, and non-resolving external-DTD handling.
  keywords: retrieval, SSRF, redirect, private IP, SHA-256, cache lock, capacity reservation, one cache scan, disk preflight, last_used_at, retention report, dry run, active path, minimum retained, quarantine, recovery, NCBI range fallback, Content-Range, compressed response, XML, archive, auxiliary member, traversal, DTD, SYSTEM, PUBLIC, entity
  link: [Open section](codebase.md#secure-retrieval-and-xml)

- id: msc6-source-services
  title: Source services
  anchor: msc6-source-services
  purpose: Explain source services and its implementation boundaries.
  keywords: source services
  link: [Open section](codebase.md#msc6-source-services)

- id: msc6-service-architecture
  title: Service integration and migration
  anchor: msc6-service-architecture
  purpose: Explain service integration and migration and its implementation boundaries.
  keywords: service integration, migration
  link: [Open section](codebase.md#msc6-service-architecture)

- id: converter-test-contracts
  title: Converter end-to-end test contracts
  anchor: converter-test-contracts
  purpose: Explain converter end-to-end test contracts and its implementation boundaries.
  keywords: converter end-to-end test contracts
  link: [Open section](codebase.md#converter-test-contracts)

- id: scoped-library-chemistry
  title: Scoped library chemistry
  anchor: scoped-library-chemistry
  purpose: Explain scoped library chemistry and its implementation boundaries.
  keywords: scoped library chemistry
  link: [Open section](codebase.md#scoped-library-chemistry)

- id: sample-library-routing
  title: Sample and library technology routing
  anchor: sample-library-routing
  purpose: Explains structured chemistry identifiers, sample-specific scRNA/Visium decisions, mixed handler construction and complete output regressions.
  keywords: GEO array categories, spotted DNA/cDNA, GSE100, preset preservation, SC3Pv2, singlecell_type, GSM5388031, GSM9254695, Visium, scRNA, TechnologyDecision, protocol registry, mixed libraries
  link: [Open section](codebase.md#sample-library-routing)

- id: miniml-v3-only-cutover
  title: MINiML compatibility and migration
  anchor: miniml-v3-only-cutover
  purpose: Regenerate saved v2 inputs and distinguish explicit legacy import from canonical decoding.
  keywords: MSC 8, v3 only, v2 rejection, migration, saved input, checkpoint separation
  link: [Open section](codebase.md#miniml-v3-only-cutover)

- id: native-archive-imports
  title: Native archive converters, metadata and enrichment
  anchor: native-archive-imports
  purpose: Native SRA/ENA discovery, parsing, preservation and conversion contracts.
  keywords: sra2json, ena2json, SRASource, ENASource, SRAParser, ENAParser, extensions.insdc, native identity, custom factor names, explicit protocols, native date precision, scoped protocol enrichment, enrichment precedence, CLI reports
  link: [Open section](codebase.md#native-archive-imports)

- id: native-archive-contract
  title: Archive extension and field mappings
  anchor: native-archive-contract
  purpose: Documents archive extension and field mappings.
  keywords: extensions.insdc, files, pool, host, dates, statistics, analysis, assembly
  link: [Open section](codebase.md#native-archive-contract)

- id: native-archive-fidelity
  title: Native retrieval and enrichment fidelity
  anchor: native-archive-fidelity
  purpose: Identity validation, structured workflow merging and file URI export contracts.
  keywords: PRJA, BioProject UID, BioSample UID, Accession search, PRJDA, assembly version, File URI, qualified units, contributor, database, sample_ref, indexed_statistics, harmonization decoding
  link: [Open section](codebase.md#native-archive-fidelity)

- id: native-archive-workflows
  title: Independent native provider workflows
  anchor: native-archive-workflows
  purpose: Documents independent native provider workflows.
  keywords: SRASource, ENASource, History, Portal, Browser, umbrella, sample expansion
  link: [Open section](codebase.md#native-archive-workflows)

- id: native-archive-cli
  title: Native command flags and outcomes
  anchor: native-archive-cli
  purpose: Documents native command flags and outcomes.
  keywords: sra2json, ena2json, include-peer, enrich-from-geo-ae, evidence, accession-resolution evidence, report, overwrite
  link: [Open section](codebase.md#native-archive-cli)

- id: native-archive-api
  title: Native Python interfaces
  anchor: native-archive-api
  purpose: Documents native python interfaces.
  keywords: ArchiveImportResult, StudyImportOutcome, SRA2JSONConverter, ENA2JSONConverter, pure parser
  link: [Open section](codebase.md#native-archive-api)

- id: native-archive-validation
  title: Native import tests
  anchor: native-archive-validation
  purpose: Documents native import tests.
  keywords: TDD, offline, live_api, accession joins, downstream exports
  link: [Open section](codebase.md#native-archive-validation)

## Stable aliases and supporting documents

- [H5AD metadata schema](codebase.md#h5ad-metadata-schema-v1): historical anchor retained; current metadata schema is 2.0.
- [Source services](codebase.md#msc6-source-services) and [service integration](codebase.md#msc6-service-architecture): older release anchors remain valid for the current architecture.
- [Rootless acceptance, 2026-07-31](rootless-acceptance-2026-07-31.md): historical live-run evidence, separate from current offline verification.
- [Source fixture catalogue](../tests/fixtures/README.md) and [test audit](../tests/AUDIT.md): fixture provenance and coverage decisions.
- [INSDC reference material](codebase.md#insdc-provider-reference-material): vendored provider specifications; this is not an installed SRA/ENA converter interface.

## Retrieve a section

From the repository root, print a section and its subsections by anchor:

```bash
python3 - architecture oop-design <<'PYCODE'
import pathlib, re, sys
text = pathlib.Path("docs/codebase.md").read_text()
for anchor in sys.argv[1:]:
    start = text.index(f'<a id="{anchor}"></a>')
    heading = re.search(r"^(#{2,6}) .+$", text[start:], re.M)
    level = len(heading.group(1))
    rest_start = start + heading.end()
    end = re.search(r"^#{1," + str(level) + r"} ", text[rest_start:], re.M)
    print(text[start:rest_start + end.start() if end else None].rstrip())
PYCODE
```

Preserve stable anchors when reorganizing material. Update this index and the
canonical guide together when an interface, workflow, or contract changes.
