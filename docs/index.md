<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# codebase.md Index

Read this routing index first. Retrieve the relevant stable section from
`docs/codebase.md`; the canonical handoff contains implementation evidence,
call chains, branch-complete workflows, and the exhaustive symbol inventory.

```bash
anchor='<anchor-id>'
python3 -c 'import pathlib,sys; p=pathlib.Path("docs/codebase.md").read_text(); a=f"<a id=\"{sys.argv[1]}\"></a>"; s=p.index(a); n=p.find("\\n<a id=",s+len(a)); print(p[s:n if n >= 0 else None])' "$anchor"
```

## Canonical routes

- id: architecture
  title: Architecture
  anchor: architecture
  purpose: Explains repository purpose, topology, ownership, and whole-system navigation.
  keywords: architecture, topology, GEO, MAGE-TAB, JSON, H5AD, TSV, CSV
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

- id: component-relationships-and-data-flow
  title: Component relationships and data flow
  anchor: component-relationships-and-data-flow
  purpose: Shows control/data exchange, ownership, lifecycle, and failure propagation.
  keywords: relationship, data flow, converter, fetcher, parser, projector, NFCoreRunner, Atlas
  link: [Open section](codebase.md#component-relationships-and-data-flow)

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
  purpose: Defines formal exports and classifies all importable production symbols by support evidence.
  keywords: public API, __all__, signature, constructor, property, method, protocol, Asset, AssetDownloader, SourcePlanner, JSON2H5ADConverter, DatasetCombinationPolicy, MINiMLMetadataProvider, MINiMLMetadataService, DatasetBundleRecoveryError, projector, obs_renames, obs_drops, failure, recovery, side effect
  link: [Open section](codebase.md#public-api-reference)

- id: atlas-v1-reader
  title: Atlas v1 reader
  anchor: atlas-v1-reader
  purpose: Defines the standalone versioned wire reader, dataset adaptation, warnings, v1 rejection, and dependency boundary.
  keywords: AtlasV1Reader, AtlasV1Error, schema_version, golden fixture, harmonized, status 2.0, metadata.packages, DatasetPackageGroup, standalone, no dependency
  link: [Open section](codebase.md#atlas-v1-reader)

- id: principal-workflows
  title: Principal workflows
  anchor: principal-workflows
  purpose: Introduces the seven conversion flows and their important terminal outcomes.
  keywords: workflow, branch, stages, pseudocode, terminal, partial result
  link: [Open section](codebase.md#principal-workflows)

- id: extension-and-change-guidance
  title: Extension and change guidance
  anchor: extension-and-change-guidance
  purpose: Defines supported extension patterns, compatibility concerns, and required tests.
  keywords: extension, projector, platform handler, asset, service, round trip, tests
  link: [Open section](codebase.md#extension-and-change-guidance)

- id: runtime-behavior
  title: Runtime behavior and packaging
  anchor: runtime-behavior
  purpose: Defines the standalone 4.0.0 distribution, exact-platform hashed runtime lock, deterministic CycloneDX SBOM, fail-closed release evidence, dependencies, services, commands, and outputs.
  keywords: version 4.0.0, pylock, SHA-256, CycloneDX, SBOM, advisory snapshot, signed artifact, security SLA, standalone, no ThematicAtlases dependency, Python, dependencies, network, Nextflow, console scripts
  link: [Open section](codebase.md#runtime-behavior)

- id: durable-artifact-publication
  title: Durable artifact publication
  anchor: durable-artifact-publication
  purpose: Defines immutable generations, the atomic current pointer, compatibility views, digest validation, and recovery failures.
  keywords: artifact, generation, current.json, fsync, rollback, recovery, bundle pointer
  link: [Open section](codebase.md#durable-artifact-publication)

- id: live-api-provider-contracts
  title: Live API provider contracts
  anchor: live-api-provider-contracts
  purpose: Owns opt-in PubMed, SRA/ENA, and BioStudies provider checks.
  keywords: live_api, PubMed, SRA, ENA, BioStudies, pagination, no retries
  link: [Open section](codebase.md#live-api-provider-contracts)

- id: request-helper
  title: Host-aware request policy
  anchor: request-helper
  purpose: Defines conservative provider defaults, process-wide hostname scheduling, in-flight limits, retries, and safe telemetry.
  keywords: RateLimitedRequester, RequestSettings, NCBIApplicationIdentity, tool, email, contact, hostname, max_in_flight, request_delay, retry, 429
  link: [Open section](codebase.md#request-helper)

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

- id: proposed-enriched-miniml-core
  title: Enriched MINiML-compatible core
  anchor: proposed-enriched-miniml-core
  purpose: Defines MSC MINiML 3.0 protocols, assay paths, occurrence-local hz groups, units, provenance, document-scoped SDRF ordering, migration, and semantic MAGE-TAB boundaries.
  keywords: MSC MINiML 3.0, MAGE-TAB, hz fields, harmonized value, canonical IDF labels, units, factors, protocols, assay paths, v2 migration, semantic round trip
  link: [Open section](codebase.md#proposed-enriched-miniml-core)

- id: miniml-package-model
  title: MSC MINiML 3.0 package API
  anchor: miniml-package-model
  purpose: Defines the immutable Python model, 3.0 codec, v2 migration, strict hz validation, projector iterators, and compatibility diagnostics.
  keywords: MSC 5, MINiML 3.0, MINiMLPackage, MINiMLCodec, MINiMLV2Migrator, HarmonizedValue, hz_disease, collision parentheses, raw preservation, source documents
  link: [Open section](codebase.md#miniml-package-model)

## Workflow routes

- id: workflow-geo2ae
  title: GEO to MAGE-TAB
  anchor: workflow-geo2ae
  purpose: Traces GEO retrieval, parsing, enrichment, construction, writing, and errors.
  keywords: geo2ae, GSE, GEO FTP, MINiML, AEConstructor, IDF, SDRF
  link: [Open section](codebase.md#workflow-geo2ae)

- id: workflow-geo2json
  title: GEO to parsed JSON
  anchor: workflow-geo2json
  purpose: Traces GEO-to-JSON parsing, guarded direct-parent publication inheritance, optional enrichment, writing, and errors.
  keywords: geo2json, GSE, GEOParser, parent publication, SubSeries, SuperSeries, reciprocal, PubMed, publication_inheritance, RelatedSeriesParseResult, related series, degraded, partial, safe error, enrichment, JSON
  link: [Open section](codebase.md#workflow-geo2json)

- id: workflow-json2ae
  title: Parsed JSON to MAGE-TAB
  anchor: workflow-json2ae
  purpose: Traces native MINiML or canonical Atlas v1 loading, filtering, validation, enrichment, typed semantic MAGE-TAB reconstruction, and writing.
  keywords: json2ae, JSON, Atlas v1, harmonized, warning, MAGE-TAB, canonical IDF labels, publication status, protocol ontology, round trip, overlay, IDF, SDRF
  link: [Open section](codebase.md#workflow-json2ae)

- id: workflow-ae2json
  title: MAGE-TAB to parsed JSON
  anchor: workflow-ae2json
  purpose: Traces bounded local, policy-approved HTTPS, and paginated BioStudies resolution through strict typed JSON output.
  keywords: ae2json, BioStudies, E-MTAB-6486, pagination, IDF, SDRF, canonical labels, legacy label aliases, factor normalization, material type, resource profile, egress, AEParser, typed model
  link: [Open section](codebase.md#workflow-ae2json)

- id: workflow-json2h5ad
  title: MINiML or Atlas JSON to H5AD
  anchor: workflow-json2h5ad
  purpose: Traces grouping, asset planning, raw/processed paths, sample projection, catalogue-only publication, fail-closed compatibility evidence, and partial failures.
  keywords: json2h5ad, Atlas, MINiML, SourcePlanner, DatasetCombinationPolicy, per-sample catalogue, no expression integration, Entrez, gene symbol, all unknown evidence, memory report, admission, force memory, 70 percent, 90 percent, sequential durable checkpoint, scope index, asset performance, typed annotations, organism, reference, modality, feature namespace, deprecated allow unverified combination, msc_miniml, packages_json, ConversionResult, BatchConversionResult, nf-core, AnnData, lazy Scanpy, processed checkpoint, resume
  link: [Open section](codebase.md#workflow-json2h5ad)

- id: workflow-json2tsv
  title: MINiML or Atlas JSON to TSV
  anchor: workflow-json2tsv
  purpose: Traces grouping, projection, fail-closed validation, ordering, and TSV output.
  keywords: json2tsv, Atlas, projector, TabularConversionResult, allow_invalid
  link: [Open section](codebase.md#workflow-json2tsv)

- id: workflow-json2obs
  title: MINiML or Atlas JSON to AnnData metadata
  anchor: workflow-json2obs
  purpose: Traces catalogue-backed observation-row aggregation without expression integration, single-sample var, and typed sample-namespaced uns publication.
  keywords: json2obs, obs CSV, row aggregation, no matrix combination, backed read, single-sample var CSV, sample-namespaced uns JSON, cell_id, feature_id, atomic, partial, memory admission, force memory, processed checkpoint, resume
  link: [Open section](codebase.md#workflow-json2obs)

## Detailed evidence routes

- id: public-api-and-callable-reference
  title: Complete callable inventory
  anchor: public-api-and-callable-reference
  purpose: Inventories public-named production callables beyond the formal export boundary.
  keywords: callable, class, method, helper, evidence gap, source
  link: [Open section](codebase.md#public-api-and-callable-reference)

- id: parsed-miniml-data-shape
  title: Parsed MINiML data shape
  anchor: parsed-miniml-data-shape
  purpose: Defines the internal JSON package shape shared by converter workflows.
  keywords: MINiML, JSON, package, series, sample, platform, enrichment
  link: [Open section](codebase.md#parsed-miniml-data-shape)

- id: h5ad-metadata-schema-v1
  title: H5AD metadata schema 1.0
  anchor: h5ad-metadata-schema-v1
  purpose: Defines canonical dotted observation columns, normalized sample values, schema markers, source-column treatment, and globally unique observation identifiers.
  keywords: H5AD, schema 1.0, obs, msc_metadata, msc_assay, typed annotations, assay parameters, sample_values, canonical columns, multivalue, observation ID
  link: [Open section](codebase.md#h5ad-metadata-schema-v1)

- id: harmonization-overrides
  title: Harmonization overrides
  anchor: harmonization-overrides
  purpose: Defines Agentic Curator envelope recognition, opt-in destination replacement, hz retention, fallback, and provenance.
  keywords: harmonization overrides, Agentic Curator, typed annotations, ECTO, PCL, exposure, cell state, profile, destination, provenance
  link: [Open section](codebase.md#harmonization-overrides)

- id: rootless-json2h5ad-runtime
  title: Rootless json2h5ad runtime
  anchor: rootless-json2h5ad-runtime
  purpose: Documents the hardened Docker/Nextflow process boundary and filesystem contract.
  keywords: rootless, Docker, Compose, socket, ACL, Nextflow, noexec
  link: [Open section](codebase.md#rootless-json2h5ad-runtime)

- id: rootless-acceptance-2026-07-31
  title: Rootless acceptance evidence
  anchor: rootless-acceptance-2026-07-31
  purpose: Routes to the successful pinned nf-core rootless acceptance summary and its portable report.
  keywords: rootless, acceptance, rnaseq 3.26.0, scrnaseq 4.2.0, H5AD, return code 0
  link: [Open section](codebase.md#rootless-acceptance-2026-07-31)

- id: reference-annotation-flow
  title: Reference and annotation flow
  anchor: reference-annotation-flow
  purpose: Explains FASTA/GTF/GFF validation, conversion, precedence, and provenance.
  keywords: genome, FASTA, GTF, GFF3, gffread, annotation, checksum
  link: [Open section](codebase.md#reference-annotation-flow)

- id: test-plan
  title: Test plan
  anchor: test-plan
  purpose: Routes maintainers to behavioral and documentation verification coverage, including the 2026-08-10 deterministic result.
  keywords: tests, 587 passed, 3 skipped, 89 subtests, pytest, no network, fake process, CLI, documentation, acceptance, regression
  link: [Open section](codebase.md#test-plan)
