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
  keywords: boundary, GEO FTP, NCBI, ENA, PubMed, OLS, BioStudies, filesystem, Nextflow, Docker, Atlas v2
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
  keywords: orchestrator, AEConstructor, AEParser, JSON2H5ADConverter, JSON2DelimitedConverter, RateLimitedRequester
  link: [Open section](codebase.md#orchestrators-and-core-types)

- id: public-api-reference
  title: Public API reference
  anchor: public-api-reference
  purpose: Defines formal exports and classifies all importable production symbols by support evidence.
  keywords: public API, __all__, signature, constructor, property, method, protocol, Asset, SourcePlanner, JSON2H5ADConverter, projector, failure, side effect
  link: [Open section](codebase.md#public-api-reference)

- id: atlas-v2-reader
  title: Atlas v2 reader
  anchor: atlas-v2-reader
  purpose: Defines the standalone versioned wire reader, dataset adaptation, warnings, v1 rejection, and dependency boundary.
  keywords: AtlasV2Reader, AtlasV2Error, schema_version, golden fixture, harmonized, metadata.packages, DatasetPackageGroup, standalone, no dependency
  link: [Open section](codebase.md#atlas-v2-reader)

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
  purpose: Defines the standalone 4.0.0 distribution, dependency/runtime requirements, external-service boundaries, commands, and outputs.
  keywords: version 4.0.0, standalone, JSONDataOutputOrchestrator, no ThematicAtlases dependency, Python, dependencies, network, Nextflow, console scripts
  link: [Open section](codebase.md#runtime-behavior)

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
  purpose: Traces the distinct GEO-to-JSON path, optional enrichment, writing, and errors.
  keywords: geo2json, GSE, GEOParser, enrichment, JSON
  link: [Open section](codebase.md#workflow-geo2json)

- id: workflow-json2ae
  title: Parsed JSON to MAGE-TAB
  anchor: workflow-json2ae
  purpose: Traces native MINiML or canonical Atlas v2 loading, filtering, validation, enrichment, round-trip restoration, construction, and writing.
  keywords: json2ae, JSON, Atlas v2, harmonized, warning, MAGE-TAB, round trip, overlay, IDF, SDRF
  link: [Open section](codebase.md#workflow-json2ae)

- id: workflow-ae2json
  title: MAGE-TAB to parsed JSON
  anchor: workflow-ae2json
  purpose: Traces local, HTTP, and BioStudies resolution through typed JSON output.
  keywords: ae2json, BioStudies, IDF, SDRF, AEParser, typed model
  link: [Open section](codebase.md#workflow-ae2json)

- id: workflow-json2h5ad
  title: MINiML or Atlas JSON to H5AD
  anchor: workflow-json2h5ad
  purpose: Traces grouping, asset planning, raw/processed paths, projection, aggregation, and partial failures.
  keywords: json2h5ad, Atlas, MINiML, ConversionResult, BatchConversionResult, nf-core, AnnData, lazy Scanpy
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
  purpose: Traces shared AnnData assembly and combined obs, optional var, and typed uns publication.
  keywords: json2obs, obs CSV, var CSV, uns JSON, cell_id, feature_id, atomic, partial
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

- id: h5ad-metadata-schema-v3
  title: H5AD metadata schema 3.0
  anchor: h5ad-metadata-schema-v3
  purpose: Defines canonical dotted observation columns, normalized sample values, schema markers, source-column treatment, and globally unique observation identifiers.
  keywords: H5AD, schema 3.0, obs, msc_metadata, sample_values, canonical columns, multivalue, observation ID, migration
  link: [Open section](codebase.md#h5ad-metadata-schema-v3)

- id: rootless-json2h5ad-runtime
  title: Rootless json2h5ad runtime
  anchor: rootless-json2h5ad-runtime
  purpose: Documents the hardened Docker/Nextflow process boundary and filesystem contract.
  keywords: rootless, Docker, Compose, socket, ACL, Nextflow, noexec
  link: [Open section](codebase.md#rootless-json2h5ad-runtime)

- id: reference-annotation-flow
  title: Reference and annotation flow
  anchor: reference-annotation-flow
  purpose: Explains FASTA/GTF/GFF validation, conversion, precedence, and provenance.
  keywords: genome, FASTA, GTF, GFF3, gffread, annotation, checksum
  link: [Open section](codebase.md#reference-annotation-flow)

- id: test-plan
  title: Test plan
  anchor: test-plan
  purpose: Routes maintainers to behavioral and documentation verification coverage.
  keywords: tests, pytest, no network, fake process, CLI, documentation, acceptance, regression
  link: [Open section](codebase.md#test-plan)
