<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# meta_standards_converter Codebase Handoff

This is the canonical handoff for the live package under
`src/meta_standards_converter`. It covers all seven conversion paths, their
runtime boundaries, extension contracts, and the evidence needed to change them
safely.

<a id="architecture"></a>
## Architecture

`meta_standards_converter` is a Python library and seven-command toolkit for
moving study metadata and expression assets among GEO MINiML, the package's
parsed JSON model, ArrayExpress MAGE-TAB, delimited sample tables, and AnnData
H5AD. CLI modules are thin batch adapters. Converter classes own use-case
orchestration; fetchers and parsers own repository-specific I/O; MAGE-TAB
constructors and H5AD/tabular projectors own output models.

```text
Users: CLI / Python / Docker / rootless Compose
                         |
                  seven converters
       +-----------------+--------------------+
       |                 |                    |
 GEO/BioStudies      local JSON/files    explicit/discovered assets
       |                 |                    |
 fetch + parse       validate/group       plan/download/run nf-core
       |                 |                    |
 enrich (optional)  project/construct     normalize/project AnnData
       +-----------------+--------------------+
                         |
      JSON | IDF/SDRF | TSV/CSV | H5AD + provenance
```

The package owns conversion policy and emitted artifacts. GEO, NCBI, ENA,
BioStudies, OLS, nf-core, Nextflow, Docker, and the Atlas v1 wire contract are
delegated boundaries. MSC implements that wire boundary locally and does not
import its producer package. See
[Project purpose and layout](#project-purpose-and-layout) for the module tree
and [Runtime behavior](#runtime-behavior) for concrete dependencies.

<a id="system-context-and-boundaries"></a>
## System context and boundaries

| Boundary | Data crossing it | Repository responsibility | Failure owner |
| --- | --- | --- | --- |
| CLI/Python callers | Accessions, paths, options, injected projectors | Validate supported combinations and return/write deterministic results | CLI records per-input failures; most library calls raise, while H5AD group aggregation records per-group failures |
| GEO FTP | GSE accession, MINiML tar/XML | Fetch safely, scope packages, and avoid leaking query data to logs | `GEOWebFetcher` and `RateLimitedRequester` |
| NCBI/ENA/PubMed/OLS | Accessions and public metadata | Rate-limit, retry, parse, and attach enrichment without owning upstream availability | Service fetcher or harmonizer; enrichment records partial failures where supported |
| BioStudies/HTTP/local MAGE-TAB | IDF/SDRF references and text | Resolve exactly one IDF plus SDRFs, parse and retain round-trip evidence | `AEWebFetcher`/`AEParser` |
| Filesystem | JSON, MAGE-TAB, matrices, H5AD, manifests, references | Validate paths, protect existing output unless overwrite is explicit, and write provenance | Converter owning the artifact |
| Nextflow/nf-core/runtime | FASTQ samplesheet, references, workflow parameters | Pin default revisions, reserve converter-owned parameters, execute and collect outputs | `NFCoreRunner`; subprocess failures become conversion failures |
| Docker/rootless Compose | Image, socket, output tree, caches | Provide a hardened supported process boundary, not a daemon | Provisioning/helper scripts and connected daemon operator |
| Atlas v1 JSON | Versioned document and harmonized dataset metadata | Validate the delegated wire format with `AtlasV1Reader` and emit organization-neutral package groups | MSC owns reading/adaptation; the producer owns atlas creation and the canonical model |

The base Python path has no persistent database. Files are the durable boundary;
in-memory dictionaries, typed MAGE-TAB models, and AnnData objects are
conversion-scoped state. Remote metadata is public study metadata, but logs
must not include URL queries, request parameters, XML, metadata payloads,
credentials, or tokens.

<a id="architectural-decisions"></a>
## Architectural decisions

<a id="decision-thin-cli-adapters"></a>
### AD-001: Keep CLI adapters thin and converters reusable

- **Status:** Observed
- **Decision:** Each console script parses a batch and delegates one input at a time to a converter class; converters expose the programmatic workflow.
- **Rationale:** Not documented.
- **Consequences:** CLI failures can be isolated per input, while Python callers retain exceptions and injectable collaborators.
- **Affected components:** `cli/*`, `converters/*`, logging configuration, and output status handling.
- **Evidence:** [`pyproject.toml`](../pyproject.toml), [`cli/geo2ae.py`](../src/meta_standards_converter/cli/geo2ae.py), and [`converters/geo2ae.py`](../src/meta_standards_converter/converters/geo2ae.py).

<a id="decision-magetab-round-trips"></a>
### AD-002: Preserve MAGE-TAB round trips beside the mapped core

- **Status:** Observed
- **Decision:** `ae2json` maps semantic MAGE-TAB content into MSC MINiML 3.0; `json2ae` regenerates ordered IDF/SDRF tables from that typed model.
- **Rationale:** Not documented.
- **Consequences:** Metadata semantics are editable and format-independent. Raw row layout is not replayed; unsafe consolidation of heterogeneous SDRF document graphs is rejected.
- **Affected components:** `AEParser`, `ae_model`, `MINiMLV1Migrator`, and `AEConstructor`.
- **Evidence:** [`ae_parser.py`](../src/meta_standards_converter/magetab/parser.py), [`ae_model.py`](../src/meta_standards_converter/magetab/semantics.py), and [`migration.py`](../src/meta_standards_converter/miniml/migration.py).

<a id="decision-expression-source-planning"></a>
### AD-003: Select expression assets before normalizing AnnData

- **Status:** Observed
- **Decision:** `json2h5ad` resolves explicit and discovered sources into a plan, then dispatches processed assets or raw nf-core execution before a common normalization/projection stage.
- **Rationale:** Not documented.
- **Consequences:** Source precedence is centralized; incompatible samples remain as per-sample H5ADs and make the aggregate result partial.
- **Affected components:** `AssetManifest`, `SourcePlanner`, `NFCoreRunner`, and `JSON2H5ADConverter`.
- **Evidence:** [`json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py) and [`test_h5ad_pipeline.py`](../tests/expression/test_h5ad_pipeline.py).

<a id="decision-fail-closed-projectors"></a>
### AD-004: Extend emitted metadata through fail-closed projector protocols

- **Status:** Observed
- **Decision:** AnnData and delimited outputs accept injected projector protocols; explicit tabular projectors replace the default MSC projection.
- **Rationale:** Not documented.
- **Consequences:** Downstream adapters can add organization-specific fields without coupling them into MSC, while collisions, invalid vector lengths, and projector errors stop unsafe output unless tabular invalid-output mode is explicit.
- **Affected components:** `AnnDataMetadataProjector`, `TabularMetadataProjector`, H5AD normalization, and delimited converters.
- **Evidence:** [`json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py), [`json2tabular.py`](../src/meta_standards_converter/converters/json2tsv.py), and their projector tests.

<a id="decision-rootless-runner"></a>
### AD-005: Isolate raw processing behind a restricted rootless runner

- **Status:** Documented
- **Decision:** The supported Compose workflow connects a locked runner to its rootless Docker socket, mounts one output tree, drops capabilities, and uses a read-only container filesystem.
- **Rationale:** Limit host access while still allowing Nextflow's Docker profile to launch containers and produce project-owner-readable artifacts.
- **Consequences:** Provisioning requires Linux root access once; raw inputs, caches, work, and outputs must live below the dedicated output root.
- **Affected components:** `Dockerfile`, `compose.yaml`, `provision-rootless-json2h5ad.sh`, and `json2h5ad-compose.sh`.
- **Evidence:** [`compose.yaml`](../compose.yaml), [`provision-rootless-json2h5ad.sh`](../scripts/provision-rootless-json2h5ad.sh), and [`test_docker_artifacts.py`](../tests/test_docker_artifacts.py).

<a id="design-invariants-and-expectations"></a>
## Design invariants and expectations

- A parsed source contains at least one package group and at least one
  convertible sample; packages remain scoped to their Series.
- JSON consumers accept native parsed MINiML or `schema_version = "1.0"` Atlas
  documents. Atlas v1 `accessions` envelopes fail with explicit pinned-v1
  guidance; MSC has no runtime/build dependency on ThematicAtlases.
- CLI batch commands continue after an input failure and return `1` if any
  input fails. Most programmatic workflows propagate exceptions; H5AD group
  aggregation is the explicit exception and converts per-group exceptions
  into `BatchConversionResult.failures`.
- `geo2ae`/`json2ae` use the same selected platform handler for IDF and SDRF.
- MAGE-TAB regeneration preserves typed-model semantics. Raw table layout is
  not part of the runtime model.
- H5AD source priority remains manifest, explicit specification, then JSON
  discovery; processed H5AD/matrix sources outrank raw FASTQ within a tier.
- Existing H5AD, tabular, and provenance outputs are protected unless the
  relevant overwrite option is explicit. H5AD sample, combined, and manifest
  artifacts publish as one staged bundle with backup/rollback on commit error.
- Metadata projectors cannot silently replace existing AnnData keys or emit
  axis vectors of the wrong length. Tabular collisions and invalid projections
  fail closed unless `allow_invalid=True`; projector-reported H5AD errors use
  the same explicit policy, while structural errors always fail.
- External calls use service-specific requester settings. Safe telemetry never
  logs credentials, URL queries, raw XML, or study payloads.
- Remote expression assets cross `RetrievalPolicy`: provider hosts and public
  addresses are checked before every request/redirect, cross-scheme redirects
  and URL credentials are rejected, and object/run/cache/disk ceilings are
  mandatory. Cache reuse requires a matching SHA-256 sidecar.
- XML provider bodies and GEO archives use the selected typed resource profile.
  XML is byte-bounded; one ordinary external SYSTEM/PUBLIC DTD may be stripped
  without resolution, while entities, internal subsets, malformed declarations,
  and misplaced declarations fail closed. GEO archives are streamed to a
  disk-preflighted temporary file and may contain safe auxiliary files and
  directories alongside exactly one expected bounded XML member.
- Converter-owned nf-core input/output/reference parameters override additional
  params JSON; Nextflow config is for infrastructure and resources.
- Rootless Compose refuses a daemon without rootless security mode and confines
  writable state to the configured output tree.

<a id="proposed-enriched-miniml-core"></a>
## Enriched MINiML-compatible core

MSC MINiML 3.0 folds semantic MAGE-TAB content into the typed package itself:
`series.protocols` owns protocol identity and details, `series.assay_paths`
owns ordered node and protocol-application paths, named values own units and
typed `annotations`, channel-level `annotations` describe harmonized scalar
channel fields such as `source`, declaration lists retain QC/replicate/normalization
semantics, and `source.documents` records source-document provenance. Raw IDF
and SDRF table layouts and the former `mage_tab` replay sidecar are deliberately
outside the runtime representation.

MAGE-TAB construction is therefore semantic and deterministic. IDF ordering
is a renderer responsibility, while assay-path order and repeated
characteristic/parameter occurrences remain data. The SDRF renderer reads the
v2 `name` field (with `tag` only as a migration fallback) and unwraps typed
ontology values such as channel `source` and `molecule` instead of serializing
their JSON object representation. IDF renderers emit the MAGE-TAB 1.1 labels
`Publication Status Term Source REF`, `Publication Status Term Accession
Number`, `Protocol Term Source REF`, and `Protocol Term Accession Number`.
The parser retains the four former MSC labels as input-only aliases and gives
the canonical rows precedence when both spellings are present. Parsing retains
characteristic units and ontology companions, keeps Protocol Contact distinct from per-application
Protocol Performer, and preserves generic IDF/SDRF comments as named comments.
Blank and external Protocol REF values and external sample names are retained
with compatibility diagnostics. Unknown layout remains outside the semantic
model and is reported rather than replayed.

SDRF row order is scoped to each `SourceDocument`. The public
`meta_standards_converter.magetab.semantics.render_miniml_assay_documents`
returns one rendered table per document, preserving repeated Sample Name paths,
protocol-application performer/date/comments, ontology companion columns,
`Unit[type]`, factor qualifiers, and repeated headers. The legacy AE constructor
has a single SDRF-table return shape, so it consolidates compatible document
graphs and raises for heterogeneous layouts instead of silently losing order or
path multiplicity.

Runtime converters accept only explicit schema `2.0` packages. Legacy or
unversioned packages enter through `MINiMLV1Migrator`, which folds supported
sidecar semantics once, converts legacy `hz_*` channel/characteristic fields
and `pre_hz_label` into typed characteristic annotations, and reports dropped
source-layout evidence. No runtime workflow emits or consumes `hz_*` fields;
harmonized values live only in typed annotation records.

**Current-state evidence:** [`ae_parser.py`](../src/meta_standards_converter/magetab/parser.py),
[`ae_model.py`](../src/meta_standards_converter/magetab/semantics.py),
[`migration.py`](../src/meta_standards_converter/miniml/migration.py),
and [`ae_constructor.py`](../src/meta_standards_converter/magetab/constructor.py).

<a id="component-relationships-and-data-flow"></a>
## Component relationships and data flow

| Source | Destination | Interface and direction | Data/lifecycle | Failure behavior |
| --- | --- | --- | --- | --- |
| CLI module | converter | `convert(...)` control call | One converter instance per command; inputs processed in order | Logs per-input exception and records non-zero status |
| GEO converters | `GEOWebFetcher` → `GEOParser` | internal calls around HTTP/XML | GSE becomes Series-scoped package dictionaries | Fetch/parse errors propagate |
| GEO/JSON converters | `MINiMLEnricher` | optional internal mutation | Adds PubMed and SRA/ENA evidence to a package | Enricher records service-specific misses where implemented |
| GEO/JSON converters | `AEConstructor` | `miniml2magetab` then optional `magetab2file` | Package becomes IDF/SDRF row collections and files | Validation/handler/write errors propagate |
| `ae2json` | `AEWebFetcher` → `AEParser` | resolve and parse | IDF/SDRF text becomes package + typed sidecar | Invalid source cardinality or MAGE-TAB fails |
| tabular converters | `JSONPackageSource` → `MINiMLMetadataProvider` → projectors | load/group, obtain format-neutral sample metadata, then `project_sample` | Dataset groups become ordered row maps | Collisions/invalid output fail closed by default |
| `JSON2H5ADConverter` | planner/downloader/runner | plan, localize, or process | Per-conversion state becomes sample AnnData | Per-sample failures retained; aggregate can be partial |
| tabular and H5AD converters | `MINiMLMetadataProvider` | public injected service calls | Study/sample identity, canonical metadata, and modality use one scientific interpretation | Provider failures propagate without partial output |
| H5AD converter | metadata projectors | `project_sample`/`project_combined` | Projection applied before each write | Collision/shape/projector error fails conversion |
| Atlas/MINiML source | `JSON2H5ADConverter` | `AtlasV1Reader` → `JSONPackageSource.load` then group conversion | Harmonized Atlas v1 datasets or ordinary packages become one conversion per dataset | Invalid versions/shapes raise before aggregation; non-harmonized states warn; per-group conversion failures are retained |
| `NFCoreRunner` | Nextflow/nf-core | subprocess/system call | Samplesheet + reference + params produce pipeline results | Exit/output failure becomes a recorded pipeline failure |

The orchestrating converter owns transient conversion state and output
decisions. Fetchers own network protocol details; parsers own input shape;
constructors/projectors own output shape; callers own retrying a failed
top-level conversion.

<a id="entrypoints-and-interfaces"></a>
## Entrypoints and interfaces

The supported public entrypoints are:

<a id="interface-cli"></a>
- seven console scripts registered in `pyproject.toml`: `geo2ae`, `geo2json`,
  `json2ae`, `ae2json`, `json2h5ad`, `json2tsv`, and `json2obs`;
<a id="interface-python"></a>
- direct Python converter classes, the twenty-four formal exports from
  `meta_standards_converter.converters`, and the four-name
  `meta_standards_converter.atlas_v1` facade;
<a id="interface-docker"></a>
- a Docker image that accepts any installed console command;
<a id="interface-rootless-compose"></a>
- rootless Compose and its two operational scripts for raw `json2h5ad`.

Every CLI parser contract, option, default, and batch status rule is listed
under [CLI](#cli). Docker and rootless process contracts are described under
[Rootless json2h5ad runtime](#rootless-json2h5ad-runtime). No HTTP server,
database service, or plugin discovery mechanism is exposed.

<a id="orchestrators-and-core-types"></a>
## Orchestrators and core types

<a id="orchestrator-atlas-v1-reader"></a>
- `AtlasV1Reader` owns version detection, structural/cross-reference/summary
  validation, harmonized-dataset selection, skipped-state diagnostics, and
  adaptation to locally owned immutable reader results. `JSONPackageSource`
  owns the subsequent MINiML package grouping/deduplication contract.

<a id="orchestrator-metadata-converters"></a>
- `geo2ae`, `geo2json`, `json2ae`, and `ae2json` coordinate metadata-only
  conversions. `AEConstructor` owns MAGE-TAB handler selection and writing;
  `AEParser` owns reverse mapping and round-trip extensions.
<a id="orchestrator-json2h5ad-converter"></a>
- `JSON2TSVConverter`, `JSON2H5ADConverter`, and `JSON2OBSConverter` own the respective JSON-origin workflows. Its manifest,
  H5AD, and AnnData-metadata methods back the three thin CLI wrappers.
  `JSON2H5ADConverter` owns the expression conversion lifecycle beneath it.
  `SourcePlanner`, `AssetManifest`, and `AssetDownloader` resolve inputs;
  `ReferenceResolver`, `AnnotationConverter`, and `NFCoreRunner` own raw-data
  execution; an injected `DatasetCombinationPolicy` owns fail-closed
  organism/reference/modality/feature-namespace evidence inspection but refuses
  implicit matrix combination. Per-sample H5ADs plus a truthful catalogue
  manifest are the expression output; result dataclasses expose complete and
  partial outcomes.
<a id="orchestrator-json2delimited-converter"></a>
- `JSON2DelimitedConverter` owns JSON grouping, sample iteration, projection,
  column ordering, validation policy, and file output. `JSON2TSVConverter`
  selects TSV or CSV through `output_format`. Both tabular and H5AD conversion
  depend on the public `MINiMLMetadataProvider` contract instead of calling
  another converter's private helpers; `MINiMLMetadataService` is the default.
<a id="core-rate-limited-requester"></a>
- `RateLimitedRequester` is the shared external-call boundary.
  `GEOWebFetcher`, `AEWebFetcher`, `PubmedWebFetcher`, and `INSDCWebfetcher`
  apply repository-specific URL and response semantics.
- `HostRequestGate` persists per-host request-start timestamps and provider
  cooldowns under an owner-only per-user runtime directory. Its `flock`
  boundary coordinates unrelated local processes as well as MSC instances;
  an unavailable implicit XDG runtime path falls back to an owner-only `/tmp`
  directory, while an explicit unsafe path still fails closed.
- `NCBIApplicationIdentity` supplies validated application tool/contact
  parameters and an optional secret API key to every PubMed/SRA E-utilities
  request without exposing their values in logs. A missing email warns once
  and retains conservative pacing.
- `OperationStatusV2`, `SafeErrorEnvelope`, and `ResourceProfile` are the
  shared status, persistence-safe error, and resource-policy vocabulary.
  `RetrievalService` consumes the resource profile behind the supported
  `AssetDownloader` facade.
<a id="core-magetab-construction"></a>
- Neutral `ae_common.ProtocolRegistry` and technology/file detection feed
  `IDFConstructor` and `SDRFConstructor`; typed MAGE-TAB
  records, and technology handlers form the MAGE-TAB construction subsystem.
  Constructor and SDRF modules now depend one-way on that neutral module; they
  contain no mutual or late imports. `ae_constructor.ProtocolRegistry` remains
  a v1 compatibility re-export of the same class.

Definitions, signatures, state, internal calls, external operations, and
failure behavior are detailed in
[Public API and callable reference](#public-api-and-callable-reference).

<a id="public-api-reference"></a>
## Public API reference

The formal support boundary is the twenty-four names in
`meta_standards_converter.converters.__all__` (including the public
`AssetDownloader`) plus the four names in
`meta_standards_converter.atlas_v1.__all__`. CLI converter classes are also
supported through their registered commands. Other non-underscored
module-level symbols are inventoried later because Python makes them
importable, but the repository contains no export declaration or compatibility
statement for them; treat those as **evidence-gap**, not stable API.

<a id="atlas-v1-reader"></a>
### Atlas v1 reader facade

- `AtlasV1Reader.load(path) -> AtlasV1ReadResult` reads UTF-8 JSON;
  `from_mapping(value) -> AtlasV1ReadResult` accepts an already decoded object.
- `AtlasV1ReadResult.datasets` contains immutable `AtlasV1Dataset` records for
  datasets whose status is `harmonized`; each record carries `dataset_id`,
  `source_repository`, `source_ordinal`, and copied MINiML-compatible metadata.
  `warnings` describes every skipped state and its document diagnostics.
- `JSONPackageSource` turns each retained dataset into one
  `DatasetPackageGroup`. Ordinary metadata is one package; metadata whose exact
  shape is `{"packages": [...]}` expands the validated object list within that
  same group. Sample deduplication and conflicting-metadata rejection then run
  across all contained packages.
- `AtlasV1Error` is the fail-closed `ValueError` subclass for unsupported
  versions, legacy unversioned envelopes, malformed collections, duplicate IDs, broken
  publication references, invalid dataset metadata, and inconsistent summary
  counts.
- A dataset's optional `harmonization.status` is parsed as the shared
  `OperationStatusV2` contract. Invalid or legacy status mappings fail closed;
  documents that omit the field remain compatible with Atlas schema 1.0.
- `SCHEMA_VERSION` is currently `"1.0"`. The producer-owned golden fixture is
  copied verbatim to `tests/fixtures/contracts/atlas-document-v1.json`; tests
  consume it without importing ThematicAtlases.
- Qualified production symbols are
  `meta_standards_converter.atlas_v1.reader.AtlasV1Dataset`,
  `meta_standards_converter.atlas_v1.reader.AtlasV1Error`,
  `meta_standards_converter.atlas_v1.reader.AtlasV1ReadResult`, and
  `meta_standards_converter.atlas_v1.reader.AtlasV1Reader`.
- Side effects are limited to reading the supplied path. The reader performs no
  network, subprocess, database, or output writes.

<a id="api-asset"></a>
### `Asset`

- **Signature:** `Asset(scope_id, path, kind, role="primary", source="json", members=(), features_path=None, barcodes_path=None, orientation="auto", md5=None, study_scope=None, reference=None, annotation_source=None, annotation_format=None, annotation_sha256=None, effective_annotation=None)`.
- **Inputs:** sample/study scope, source path/kind, role and optional matrix,
  checksum, reference, and annotation metadata.
- **Outputs:** frozen description of one processed, matrix, 10x, or raw
  expression source.
- **Failures:** construction performs no custom validation; consuming planners
  and converters validate supported combinations.
- **Side effects:** none.
- **Support:** formal export for source-planner and adapter integrations.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-json2h5ad-converter"></a>
### `JSON2H5ADConverter`

- **Signature:** `JSON2H5ADConverter(planner=None, pipeline_runner=None, downloader=None, metadata_projectors=None, package_source=None, retrieval_policy=None, resource_profile="standard", resource_overrides=None, metadata_service=None, combination_policy=None)`; `convert(...)` and `convert_source(...)` own the documented expression workflow.
- Dataset, study, sample, checkpoint, and output identities are validated as safe single path components before publication. `series.iid` is the canonical native-package study identity, so an ArrayExpress IID is not displaced by an earlier GEO secondary accession.
- The injected/default combination policy owns multi-sample scientific
  compatibility evidence. All-unknown dimensions are missing, Entrez IDs and
  gene symbols are distinct namespaces, and `combine()` fails with guidance
  because catalogue conversion never treats an outer join as integration. The
  converter facade retains source, normalization, and publication orchestration.
- `DatasetBundleRecoveryError` preserves the original publication error, rollback errors, and surviving recovery paths when an overwrite cannot be fully restored.
- **Inputs:** native MINiML or Atlas v1 JSON, source/reference/runtime options,
  and optional public collaborators.
- **Outputs:** single or batch conversion results plus a transactional
  per-sample H5AD catalogue bundle and manifest declaring no expression
  integration.
- **Failures:** structural, validation, filesystem, external process, and
  publication failures follow the strict/permissive and batch contracts in
  [JSON to H5AD](#workflow-json2h5ad).
- **Side effects:** may read/download assets, run nf-core, and publish staged
  H5AD/provenance artifacts.
- **Support:** formal export and composition boundary.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

The non-exported injection implementation and its typed scientific rejection
are importable as
`meta_standards_converter.converters.dataset_combination.DatasetCombinationPolicy`
and
`meta_standards_converter.converters.dataset_combination.DatasetCompatibilityError`.
They are documented implementation seams rather than additions to
`converters.__all__`.

<a id="api-source-planner"></a>
### `SourcePlanner`

- **Signature:** `SourcePlanner()`; `discover(packages) -> list[Asset]` and
  `plan(packages, explicit_assets=None, force_reprocess=False) -> dict[str, Asset]`.
- **Inputs:** MINiML packages and optional explicit assets.
- **Outputs:** one ranked `Asset` per sample; subclasses may specialize
  discovery while retaining converter-owned lifecycle policy.
- **Failures:** unsupported or missing sources raise `ValueError`.
- **Side effects:** planning does not download, execute, or publish.
- **Support:** formal export and source-planning extension point.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-anndata-metadata-projection"></a>
### `AnnDataMetadataProjection`

- **Signature:** `AnnDataMetadataProjection(obs={}, var={}, uns={}, obs_renames={}, obs_drops=(), warnings=(), errors=())`.
- **Inputs:** mappings of additions for AnnData axes/unstructured metadata, atomic observation rename/drop requests, and warning and validation-error strings.
- **Outputs:** frozen projector-result dataclass.
- **Failures:** construction performs no validation; application rejects missing transform sources, duplicate/colliding targets, rename/drop overlap, addition collisions, and wrong-length axis values.
- **Side effects:** none.
- **Support:** formal export.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-anndata-projection-error"></a>
### `AnnDataProjectionError`

- **Signature:** `AnnDataProjectionError(errors)`; subclass of `ValueError`.
- **Inputs:** projector-reported validation error strings.
- **Outputs:** fail-closed exception exposing the normalized `errors` tuple.
- **Failures:** raised before final bundle publication unless `allow_invalid=True`.
- **Side effects:** none; staged files are discarded.
- **Support:** formal export.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-anndata-metadata-projector"></a>
### `AnnDataMetadataProjector`

- **Signature:** structural `Protocol` with `project_sample(*, adata, context) -> AnnDataMetadataProjection`.
- **Inputs:** current sample AnnData plus one context. Former
  `project_combined` methods on third-party objects are ignored because no
  combined expression object exists.
- **Outputs:** metadata additions and warnings.
- **Failures:** projector exceptions, collisions, invalid vectors, and wrong result types fail conversion.
- **Side effects:** the converter, not the projector contract, owns applying returned additions.
- **Support:** formal export and injection extension point; it has no `@runtime_checkable`, so runtime protocol checks are unsupported.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-metadata-projection-context"></a>
### `MetadataProjectionContext`

- **Signature:** `MetadataProjectionContext(sample, package, study_accession, sample_accession, asset, base_metadata)`.
- **Inputs:** read-only mappings, identifiers, selected `Asset`, and base metadata.
- **Outputs:** frozen H5AD-projector context.
- **Failures:** no custom validation.
- **Side effects:** none.
- **Support:** formal export; these are all current fields.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

<a id="api-json-data-output-orchestrator"></a>
### `JSON2OBSConverter`

- **Signature:** `JSON2OBSConverter(h5ad_converter=None, components=None).convert(...)`.
- **Inputs:** JSON source, output destination, optional per-sample sidecars, and H5AD conversion options.
- **Outputs:** `AnnDataMetadataBatchResult` with aggregated observation tables and optional sample components.
- **Failures:** Preserves group errors, overwrite refusal, and atomic publication/recovery errors.
- **Side effects:** Converts selected assets, then publishes tables through the artifact bundle service.
- **Support:** Public MSC 6 API. TSV manifest publication is owned by `JSON2TSVConverter.export_manifest`; H5AD publication by `JSON2H5ADConverter.convert`.
- **Source:** [`converters/json2obs.py`](../src/meta_standards_converter/converters/json2obs.py).

<a id="api-anndata-metadata-export-result"></a>
### `AnnDataMetadataExportResult`

- **Signature:** `AnnDataMetadataExportResult(dataset_id, obs, var, uns, obs_path, var_path, uns_path, manifest_path, warnings=(), errors=(), partial=False)`.
- **Inputs:** one dataset's row-aggregated sample observation metadata, optional published paths, and diagnostics.
- **Outputs:** in-memory `obs`, optional `var`/`uns`, artifact locations, and a compact `to_dict()` summary.
- **Failures:** construction performs no custom validation; the orchestrator validates and serializes its data.
- **Side effects:** none.
- **Support:** formal export and successful per-dataset result for `json2obs`.
- **Source:** [`converters/json2obs.py`](../src/meta_standards_converter/converters/json2obs.py).

<a id="api-anndata-metadata-batch-result"></a>
### `AnnDataMetadataBatchResult`

- **Signature:** `AnnDataMetadataBatchResult(source, conversions=(), warnings=(), failures=())`.
- **Inputs:** source path, completed dataset exports, cross-dataset warnings, and keyed failures.
- **Outputs:** immutable batch status, `partial` state, and compact `to_dict()` summary.
- **Failures:** construction performs no custom validation; conversion failures are retained in `failures`.
- **Side effects:** none.
- **Support:** formal export and multi-dataset result for `json2obs`.
- **Source:** [`converters/json2obs.py`](../src/meta_standards_converter/converters/json2obs.py).

<a id="api-json2tsv-converter"></a>
### `JSON2TSVConverter`

- **Signature:** `JSON2TSVConverter(metadata_projectors=None, package_source=None, metadata_service=None, *, output_format="tsv")`; inherited `convert_source(source, destination, *, allow_invalid=False, overwrite=False) -> TabularConversionResult`.
- **Inputs:** parsed MINiML JSON or canonical Atlas v1, TSV/CSV format, and destination.
- **Outputs:** selected delimited file and result metadata.
- **Failures:** invalid formats, source, projector, collision, fail-closed diagnostic, and protected-output errors propagate.
- **Side effects:** creates the destination parent and writes TSV or CSV.
- **Support:** formal export.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

<a id="api-miniml-metadata-service"></a>
### `MINiMLMetadataProvider` and `MINiMLMetadataService`

- **Signature:** `MINiMLMetadataProvider` is the structural service contract;
  `MINiMLMetadataService()` is its default stateless implementation.
- **Inputs:** parsed MINiML package/sample mappings.
- **Outputs:** canonical study/sample identity, tuple-valued or rendered sample
  metadata, platform resolution, normalized values, and modality.
- **Failures:** malformed collaborator inputs and ontology protocol mapping
  failures propagate to the owning converter.
- **Side effects:** none; the service does not read files, access the network,
  execute processes, or publish artifacts.
- **Support:** both names are formal exports and the supported injection seam
  shared by tabular and AnnData output adapters.
- **Source:** [`metadata/interpretation.py`](../src/meta_standards_converter/metadata/interpretation.py).

<a id="api-msc-metadata-projector"></a>
### `MSCMetadataProjector`

- **Signature:** `MSCMetadataProjector()` and `project_sample(*, context: TabularMetadataContext) -> TabularMetadataProjection`.
- **Inputs:** one sample context.
- **Outputs:** canonical dotted `msc.*` values and ordered base columns.
- **Failures:** converter validation applies to the returned projection.
- **Side effects:** none.
- **Support:** formal export and default projector when no explicit projectors are supplied.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

<a id="api-tabular-conversion-result"></a>
### `TabularConversionResult`

- **Signature:** `TabularConversionResult(row_count, columns, dataset_ids, warnings=(), errors=(), output_path=None)`.
- **Inputs:** immutable output summary values.
- **Outputs:** frozen result; `partial` is `True` exactly when `errors` is non-empty.
- **Failures:** no custom validation.
- **Side effects:** none.
- **Support:** formal export; `partial` is its public property.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

<a id="api-tabular-metadata-context"></a>
### `TabularMetadataContext`

- **Signature:** `TabularMetadataContext(package, sample, dataset_id, study_accession, sample_accession, base_metadata)`.
- **Inputs:** group/package/sample mappings, identifiers, and base metadata.
- **Outputs:** frozen context passed to every tabular projector.
- **Failures:** no custom validation.
- **Side effects:** none.
- **Support:** formal export; these are all current fields.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

<a id="api-tabular-metadata-projection"></a>
### `TabularMetadataProjection`

- **Signature:** `TabularMetadataProjection(values, columns=(), warnings=(), errors=())`.
- **Inputs:** projected values, preferred column order, and diagnostics.
- **Outputs:** frozen per-sample projection.
- **Failures:** converter rejects wrong types/collisions; errors raise `TabularProjectionError` unless `allow_invalid=True`.
- **Side effects:** none.
- **Support:** formal export.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

<a id="api-tabular-metadata-projector"></a>
### `TabularMetadataProjector`

- **Signature:** structural `Protocol` with `project_sample(*, context: TabularMetadataContext) -> TabularMetadataProjection`.
- **Inputs:** one immutable sample context.
- **Outputs:** projected values, order, warnings, and errors.
- **Failures:** projector exceptions propagate; converter validates type and collisions.
- **Side effects:** none required.
- **Support:** formal export/injection extension point; not runtime-checkable.
- **Source:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py).

Supported production symbols are inventoried below by qualified name. Their
signatures, inputs, outputs, exceptions, side effects, and decisive internal or
external calls are documented in the linked legacy callable sections that
follow this canonical overview.

- MAGE-TAB:
  `meta_standards_converter.magetab.constructor.validate_platform_handler`,
  `meta_standards_converter.magetab.constructor.ProtocolRegistry`,
  `meta_standards_converter.magetab.constructor.AEConstructor`,
  `meta_standards_converter.magetab.idf.IDFConstructor`,
  `meta_standards_converter.magetab.semantics.MAGETabModelError`,
  `meta_standards_converter.magetab.semantics.build_model`,
  `meta_standards_converter.magetab.semantics.validate_model`,
  `meta_standards_converter.magetab.semantics.render_model`,
  `meta_standards_converter.magetab.semantics.overlay_core`,
  `meta_standards_converter.magetab.parser.normalized_label`,
  `meta_standards_converter.magetab.parser.AEParser`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFAttr`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFNode`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFEdge`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFPath`,
  `meta_standards_converter.magetab.sdrf.constructor.ColumnGroup`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFAudit`,
  `meta_standards_converter.magetab.sdrf.constructor.SDRFConstructor`,
  `meta_standards_converter.magetab.sdrf.constructor.normalized_extension`,
  `meta_standards_converter.magetab.sdrf.constructor.classify_file`,
  `meta_standards_converter.sources.magetab.TextResource`,
  `meta_standards_converter.sources.magetab.MAGETabInput`, and
  `meta_standards_converter.sources.magetab.AEWebFetcher`.
- CLI:
  `meta_standards_converter.cli.ae2json.main`,
  `meta_standards_converter.cli.geo2ae.main`,
  `meta_standards_converter.cli.geo2json.main`,
  `meta_standards_converter.cli.json2ae.main`,
  `meta_standards_converter.cli.json2h5ad.main`,
  `meta_standards_converter.cli.json2tsv.main`,
  `meta_standards_converter.cli.json2obs.main`,
  `meta_standards_converter.cli.common.add_platform_handler_arguments`,
  `meta_standards_converter.cli.common.print_platform_handlers`,
  `meta_standards_converter.cli.common.add_logging_arguments`,
  `meta_standards_converter.cli.common.log_level`, and
  `meta_standards_converter.cli.common.configure_logging`.
- Conversion:
  `meta_standards_converter.converters.ae2json.AE2JSONConverter`,
  `meta_standards_converter.converters.geo2ae.GEO2AEConverter`,
  `meta_standards_converter.converters.geo2json.GEO2JSONConverter`,
  `meta_standards_converter.converters.json2ae.JSON2AEConverter`,
  `meta_standards_converter.expression.assets.Asset`,
  `meta_standards_converter.metadata.projection.anndata.MetadataProjectionContext`,
  `meta_standards_converter.metadata.projection.anndata.AnnDataMetadataProjection`,
  `meta_standards_converter.metadata.projection.anndata.AnnDataProjectionError`,
  `meta_standards_converter.metadata.projection.anndata.AnnDataMetadataProjector`,
  `meta_standards_converter.expression.assets.AssetManifest`,
  `meta_standards_converter.expression.assets.AssetDownloader`,
  `meta_standards_converter.expression.catalogue.PipelineRun`,
  `meta_standards_converter.expression.catalogue.ConversionResult`,
  `meta_standards_converter.expression.catalogue.BatchConversionResult`,
  `meta_standards_converter.expression.catalogue.RawProcessingResult`,
  `meta_standards_converter.expression.references.ReferenceResolver`,
  `meta_standards_converter.expression.references.AnnotationConverter`,
  `meta_standards_converter.expression.nfcore.NFCoreRunner`,
  `meta_standards_converter.expression.planning.SourcePlanner`,
  `meta_standards_converter.converters.json2h5ad.JSON2H5ADConverter`, and
  `meta_standards_converter.converters.json2h5ad.JSON2H5ADConverter`.
- Tabular and JSON source:
  `meta_standards_converter.metadata.projection.tabular.TabularMetadataContext`,
  `meta_standards_converter.metadata.projection.tabular.TabularMetadataProjection`,
  `meta_standards_converter.metadata.projection.tabular.TabularMetadataProjector`,
  `meta_standards_converter.metadata.projection.tabular.TabularConversionResult`,
  `meta_standards_converter.metadata.projection.tabular.TabularProjectionError`,
  `meta_standards_converter.metadata.projection.tabular.MSCMetadataProjector`,
  `meta_standards_converter.converters.json2tsv.JSON2DelimitedConverter`,
  `meta_standards_converter.converters.json2tsv.JSON2TSVConverter`,
  `meta_standards_converter.converters.json2tsv.JSON2TSVConverter`,
  `meta_standards_converter.converters.json2obs.JSON2OBSConverter`,
  `meta_standards_converter.expression.components.AnnDataMetadataExportResult`,
  `meta_standards_converter.expression.components.AnnDataMetadataBatchResult`,
  `meta_standards_converter.sources.json.DatasetPackageGroup`,
  `meta_standards_converter.sources.json.SourceLoadResult`, and
  `meta_standards_converter.sources.json.JSONPackageSource`,
  `meta_standards_converter.metadata.interpretation.MINiMLMetadataProvider`,
  and `meta_standards_converter.metadata.interpretation.MINiMLMetadataService`.
- Fetch, parse, enrich, harmonize, and helpers:
  `meta_standards_converter.metadata.enrichment.MINiMLEnricher`,
  `meta_standards_converter.miniml.geo_parser.GEOParser`,
  `meta_standards_converter.miniml.geo_parser.RelatedSeriesParseResult`,
  `meta_standards_converter.sources.geo.GEOWebFetcher`,
  `meta_standards_converter.metadata.ontology_mappings.GEO2OLS`,
  `meta_standards_converter.metadata.ontology_mappings.Harmonizer`,
  `meta_standards_converter.metadata.ontology_mappings.Pubmed2OLS`,
  `meta_standards_converter.helpers.json_helper.JSONHandler`,
  `meta_standards_converter.helpers.request_helper.HostRequestCooldownDeferred`,
  `meta_standards_converter.helpers.request_helper.HostRequestGate`,
  `meta_standards_converter.helpers.request_helper.NCBIApplicationIdentity`,
  `meta_standards_converter.helpers.request_helper.RequestSettings`,
  `meta_standards_converter.helpers.request_helper.RateLimitedRequester`,
  `meta_standards_converter.sources.insdc.INSDCWebfetcher`, and
  `meta_standards_converter.sources.pubmed.PubmedWebFetcher`.

<a id="principal-workflows"></a>
## Principal workflows

<a id="workflow-geo2ae"></a>
### `geo2ae`: GEO to MAGE-TAB

```text
GSE -> fetch MINiML --failure--> exception
    -> parse packages --failure--> exception
    -> enrich each --failure----> exception
    -> build IDF/SDRF --failure-> exception
    -> [out?] write files / else return in memory
```

1. The CLI validates accessions/options and calls `geo2ae.convert`.
2. `GEOWebFetcher.fetch_gse_miniml` performs the GEO FTP HTTP retrieval.
3. `GEOParser.parse` scopes packages and optionally traverses related Series.
4. Each package is enriched and passed to `AEConstructor.miniml2magetab`.
5. `out` writes IDF/SDRF; otherwise only the in-memory list is returned.

Pseudocode: `fetch -> parse -> for package: enrich -> construct -> [write] -> list`.

**Evidence:** [`converters/geo2ae.py`](../src/meta_standards_converter/converters/geo2ae.py), [`cli/geo2ae.py`](../src/meta_standards_converter/cli/geo2ae.py), and [`geo_webfetcher.py`](../src/meta_standards_converter/sources/geo.py).

<a id="workflow-geo2json"></a>
### `geo2json`: GEO to parsed JSON

```text
GSE -> fetch --failure--> exception
    -> parse --failure--> exception
    -> [enrich?] yes -> [one direct parent publication?]
                         yes -> bounded parent fetch + reciprocal/unique-PMID guard
                         no  -> retain child publication state
                      -> enrich each --failure--> exception
                  no --------------------------> retain parsed packages
    -> [out?] JSON file / else return only
```

1. The CLI calls `geo2json.convert` once per accession and continues after failures.
2. GEO fetch and parsing are shared with `geo2ae`.
3. With enrichment enabled, a child that has no direct publication and exactly one
   `SubSeries of` parent may inherit exactly one parent PubMed ID. The parent must
   reciprocally identify the child as a `SuperSeries of` member. Retrieval failure,
   non-reciprocal linkage, multiple parents, or multiple PubMed IDs leaves the child
   unchanged. This lookup does not recursively traverse related studies or add the
   parent as another output package.
4. Inherited publication provenance is retained under
   package `extensions.publication_inheritance`; the ordinary PubMed enricher then
   resolves title, DOI, authors, and status from the inherited identifier.
5. `enrich=False` bypasses parent-publication, PubMed, and SRA/ENA enrichment.
6. `json2file` creates the output directory and writes `{GSE}.json` when requested.

Pseudocode: `packages = parse(fetch(gse)); [inherit guarded parent PMID]; [enrich packages]; [write]; return`.

**Evidence:** [`converters/geo2json.py`](../src/meta_standards_converter/converters/geo2json.py) and [`cli/geo2json.py`](../src/meta_standards_converter/cli/geo2json.py).

<a id="workflow-json2ae"></a>
### `json2ae`: MINiML/Atlas v1 JSON to MAGE-TAB

```text
path -> AtlasV1Reader/JSONPackageSource -> invalid/version/v1 -> exception
                                      \-> non-harmonized dataset -> warning + skip
     -> [enrich?] -> construct each package --failure--> exception
     -> [out?] IDF/SDRF files / else in-memory MAGE-TAB list
```

1. `JSONPackageSource` accepts one parsed MINiML object, a non-empty package
   list, or a canonical Atlas v1 document.
2. `AtlasV1Reader` validates version, IDs, references and summary; it retains
   `harmonized` datasets and emits warnings for other states and their
   diagnostics. V1 envelopes are rejected. No remaining groups raises
   `ValueError`.
3. Exact Atlas metadata wrappers `{"packages": [...]}` expand within the
   dataset's single group; non-object entries fail before conversion.
4. Every retained package must be an object with a usable study accession,
   and all packages are validated before collaborator calls.
5. Optional enrichment precedes `AEConstructor.miniml2magetab`.
6. `AEConstructor` renders deterministic IDF/SDRF tables from native protocols,
   assay paths, named characteristics, units, and typed annotations. It does
   not replay raw source tables or consult a `mage_tab` sidecar, and it emits
   canonical MAGE-TAB 1.1 publication-status and protocol ontology companion
   labels rather than the former MSC aliases.
7. `out` controls writing; construction errors propagate.

Pseudocode: `validate(flatten(source.load(path))); warn(skipped); for package:
[enrich] -> construct -> [write]; return`.

**Evidence:** [`converters/json2ae.py`](../src/meta_standards_converter/converters/json2ae.py), [`ae_constructor.py`](../src/meta_standards_converter/magetab/constructor.py), and [`ae_model.py`](../src/meta_standards_converter/magetab/semantics.py).

<a id="workflow-ae2json"></a>
### `ae2json`: MAGE-TAB to parsed JSON

```text
local/HTTP/accession -> resolve IDF + SDRF(s) --failure--> exception
                     -> parse/model/round-trip --failure--> exception
                     -> [out?] accession JSON / else return package list
```

1. `AEWebFetcher.resolve` accepts a bounded local IDF, a policy-approved HTTPS IDF, or a BioStudies accession and optional SDRF overrides.
2. Every API, IDF, and SDRF response is host/public-address checked, manually redirected under the selected profile, streamed under the per-file ceiling, and charged to the aggregate run limit. Local reads stop at the same ceiling.
3. Resolution requires exactly one IDF and at least one SDRF.
4. `AEParser.parse` maps core fields and retains typed/source provenance. It registers referenced accession databases, maps non-MINiML source material types to lossless characteristics, keeps extract molecules in the typed molecule field, and maps MAGE-TAB factors such as `compound` to the nearest XSD factor while retaining the original type value. Canonical publication/protocol ontology companion rows and the four former MSC aliases map to the same typed values; regenerated output is always canonical.
5. The frozen E-MTAB-6486 IDF/SDRF contract exercises the ENA secondary accession, repeated material columns, and `compound` factor through strict `MINiMLCodec` validation.
6. `out` writes a sanitized accession filename; otherwise no file is created.

Pseudocode: `resolved = fetcher.resolve(source); package = parser.parse(resolved); [write]; return [package]`.

**Evidence:** [`converters/ae2json.py`](../src/meta_standards_converter/converters/ae2json.py), [`ae_webfetcher.py`](../src/meta_standards_converter/sources/magetab.py), and [`ae_parser.py`](../src/meta_standards_converter/magetab/parser.py).

<a id="workflow-json2h5ad"></a>
### `json2h5ad`: MINiML/Atlas v1 JSON to H5AD

The source may be ordinary parsed MINiML JSON (one object or list) or a
canonical Atlas v1 document. `AtlasV1Reader` validates and selects harmonized
datasets; `JSONPackageSource` adapts each selected dataset to one package group
and carries skipped-state warnings. An exact `metadata.packages` wrapper
expands multiple MINiML packages inside that group without changing the Atlas
dataset ID or output containment scope.

```text
path -> missing/invalid/no groups --------------------------> exception
     -> one group via convert -> package conversion --------> ConversionResult
     -> multiple groups via convert/any via convert_source
          -> each group in child directory when multiple
          -> success ---------------------------------------> conversions[id]
          -> exception -------------------------------------> failures[]; continue
          -> aggregate -------------------------------------> BatchConversionResult
package conversion -> processed normalize / raw reference + nf-core
          -> estimate memory; oversized sample retained as partial failure
          -> admitted samples checkpoint sequentially; release each matrix
          -> one H5AD per sample; no expression-matrix join
          -> catalogue manifest declaring expression_integration=none
```

1. `convert(json_path, ...)` returns `ConversionResult` for one group or
   `BatchConversionResult` for multiple groups.
2. `convert_source(json_path, ...)` always returns `BatchConversionResult`,
   including for one group.
3. With multiple groups, `_convert_groups` places each group in an output-root
   child directory named for `dataset_id`; one group uses the root directly.
4. Planning applies manifest, explicit, then discovered asset precedence. It
   indexes candidates once by sample/study scope, avoiding a full asset scan
   for every sample while preserving source order and rank behavior.
5. Processed assets normalize directly; ordinary H5AD and delimited paths do
   not import Scanpy, while 10x HDF5/MTX branches import it lazily. Raw assets
   call reference resolution and Nextflow/nf-core through `NFCoreRunner`.
6. Per-sample failures and allowed projector errors can make
   `ConversionResult.partial`; per-group exceptions are caught in
   `BatchConversionResult.failures`. Organism, reference, modality, and feature
   namespace heterogeneity is catalogued rather than joined and is not itself a
   conversion failure. `allow_unverified_combination=True` /
   `--allow-unverified-combination` is retained as a deprecated compatibility
   option and ignored with a warning.
7. Invalid path/source/no-group/no-sample and unsafe dataset-ID conditions raise before aggregation; successful
   groups and diagnostics survive later failures.
8. Bound MAGE-TAB Parameter Values are projected to dotted `obs` columns and
   every typed attribute occurrence is retained in `uns["msc_mage_tab"]`.
   Raw named characteristics use `msc.characteristics.<name>`; their typed
   ontology annotations use separate
   `msc.characteristics.harmonized_<field>` value/ID/ontology columns. Typed
   organism annotations take precedence over raw organism labels in each
   sample artifact.
9. `uns["msc_miniml"]` retains the flattened query table and provenance plus
   `packages_json`, a deterministic JSON encoding of the complete MSC MINiML
   package list. JSON encoding is intentional because HDF5 cannot represent
   heterogeneous lists of nested MINiML objects natively; `json.loads`
   reconstructs the data-model shape losslessly.
10. Before loading an expression matrix, conversion estimates peak resident
   memory. Normal admission is `min(profile ceiling, 70% available RAM)`, where
   `standard` is 8 GiB and `large` is 32 GiB. Oversized samples are skipped and
   recorded in `memory_report`, making the result partial. Only
   `--resume --force-memory` bypasses the fixed ceiling; its non-bypassable
   limit is 90% of current host/cgroup availability.
11. Every admitted sample is normalized and atomically checkpointed before its
   AnnData matrix is released. The default checkpoint root is
   `{out}/.processed/{study}`; `--processed-checkpoint-dir` overrides it. The
   fingerprint covers source JSON, sample, asset, orientation, and converter
   version. With `--resume`, valid matching checkpoints are copied into the
   catalogue without materialising their expression matrices; incomplete,
   corrupt, colliding, or stale checkpoint pairs are ignored.

Pseudocode: `load -> if one and convert: convert_packages; else for group: try convert_packages into child; except record; return batch`.

**Evidence:** [`JSON2H5ADConverter`](../src/meta_standards_converter/converters/json2h5ad.py), [`JSONPackageSource`](../src/meta_standards_converter/sources/json.py), and [`cli/json2h5ad.py`](../src/meta_standards_converter/cli/json2h5ad.py).

<a id="workflow-json2tsv"></a>
### `json2tsv`: MINiML/Atlas JSON to sample manifest

```text
source -> load/group --failure--> exception
       -> project each sample
          -> bad type/collision ---------------------------> exception
          -> errors + !allow_invalid ----------------------> TabularProjectionError
          -> errors + allow_invalid -----------------------> partial result
       -> destination exists + !overwrite ----------------> FileExistsError
       -> write selected TSV/CSV + JSON manifest ----------> result
```

1. `JSONPackageSource.load` accepts native MINiML or canonical Atlas v1 data.
2. The converter builds base metadata and invokes every projector per sample.
3. Preferred columns precede sorted extras; diagnostics are deduplicated.
4. Validation is fail-closed unless `allow_invalid=True`.
5. The selected table and result JSON are staged and published as one bundle.

Pseudocode: `load -> project -> validate -> order -> protect -> write selected delimiter + result JSON`.

**Evidence:** [`converters/json2tsv.py`](../src/meta_standards_converter/converters/json2tsv.py), [`sources/json.py`](../src/meta_standards_converter/sources/json.py), and [`cli/json2tsv.py`](../src/meta_standards_converter/cli/json2tsv.py).

<a id="workflow-json2obs"></a>
### `json2obs`: MINiML/Atlas JSON to AnnData metadata

```text
JSON + expression assets -> per-sample H5AD catalogue
       -> no sample H5AD ----------------------------------> failure
       -> concatenate obs metadata rows only --------------> .obs.csv
       -> include_var + one sample ------------------------> .var.csv
       -> include_var + multiple samples ------------------> explicit failure
       -> include_uns -------------------------------------> typed sample-namespaced .uns.json
       -> atomic component bundle + JSON result manifest
```

1. Asset resolution, raw processing, normalization, projectors, and smart IDs
   exactly match `json2h5ad`; expression matrices are not combined.
2. Sample `obs` tables are read in backed mode and row-aggregated with
   `cell_id`; source and canonical names remain unchanged.
3. Optional `var` uses `feature_id` only for a single-sample catalogue because
   a union feature table would imply an unperformed integration. Optional
   `uns` uses tagged JSON for nested mappings, arrays, and DataFrames and is
   namespaced by sample for multi-sample catalogues.
4. Every batch conversion publishes beneath its dataset-ID directory, even when other groups fail and only one dataset succeeds.
5. The CLI prints only a compact JSON summary to stdout, sends logs to stderr, and returns status `1` for partial/failure outcomes.

Pseudocode: `catalogue -> backed-read each sample metadata -> concatenate obs rows -> serialize selected components -> atomic publish -> result`.

**Evidence:** [`converters/json2obs.py`](../src/meta_standards_converter/converters/json2obs.py) and [`cli/json2obs.py`](../src/meta_standards_converter/cli/json2obs.py).

<a id="extension-and-change-guidance"></a>
## Extension and change guidance

- Add a converter workflow behind a reusable converter first, then register a
  thin CLI and update `pyproject.toml`, parser-contract tests, README interface
  tables, this handoff, and the index.
- Add organization-specific H5AD or table metadata through the public projector
  protocols. Preserve collision, length, ordering, warning, and fail-closed
  contracts; do not import downstream organization packages here.
- Add a MAGE-TAB technology through the shared platform registry and matched
  IDF/SDRF handler selection. Update handler-list tests and document the new
  metadata-based detection rule.
- Add an expression source by extending `Asset` classification, precedence,
  localization, and normalization together; test explicit, manifest, and
  discovery paths plus partial failure and overwrite behavior.
- Add an external service behind `RateLimitedRequester` with explicit timeout,
  delay, retry, safe-logging, and response validation tests.
- Any round-trip model change requires unchanged restoration, edited overlay,
  duplicate-header alignment, and fingerprint compatibility tests.
- Any raw-runner change requires Docker artifact tests, rootless-daemon
  rejection, mount/ACL review, and corresponding README/codebase updates.

<a id="insdc-provider-reference-material"></a>
## INSDC provider reference material

The repository vendors offline, source-faithful provider evidence for designing
future SRA- and ENA-native conversions. This reference material does not implement `sra2json` or `ena2json`; no new CLI, converter, or runtime endpoint is registered.

| Route | Contents and contract |
|---|---|
| [SRA snapshot](sra/README.md) | All eight `SRA.*.xsd` files from INSDC SRA 1.5; NCBI BioSample/BioProject schemas; SRA EInfo and BioSample catalogues; NCBI composite SRA, BioSample, BioProject, and available PubMed fixtures; integrity manifest |
| [SRA expected fields](sra/expected-fields.md) | Study/sample/experiment/run hierarchy, XPath/cardinality/value classes, linked-record joins, strict-validation limits, current MSC subset, precedence, and future MAGE-TAB coverage |
| [ENA snapshot](ena/README.md) | All eight in-scope `ENA.*.xsd` files, four OpenAPI documents, nine result-type field catalogues, 18 controlled vocabularies, checklist evidence, Browser/Portal/file-report/Xref/Taxonomy fixtures, integrity manifest |
| [ENA expected fields](ena/expected-fields.md) | Endpoint/column and XML contracts, accession classes, semicolon-aligned files, missing-value semantics, current MSC subset, precedence, and future MAGE-TAB coverage |
| [Checklist availability](ena/checklist-availability.md) | Snapshot reconciliation: Portal declared 47 checklist IDs, while the documented Browser route supplied 31 XML records and returned HTTP 404 for 16 |

Two fixture chains cover complementary conditions: `SRX017289` links
`SRP002056` → `SRS011830` → `SRX017289` → `SRR037073` with BioProject,
BioSample, GEO, and PMID `20133686`; `SRX7812918` links `SRP250911` →
`SRS6225446` → `SRX7812918` → `SRR11192680` with BioProject/BioSample,
paired files, geographic sample attributes, and no publication. Provider files
remain byte-for-byte snapshots. Each provider manifest records URL, retrieval
timestamp, content type, available HTTP validators, SHA-256, local path, and
fixture accession; offline tests recalculate every digest.

The active implementation boundary remains deliberately narrow.
`MINiMLEnricher` discovers SRA accessions from GEO sample relations and uses
NCBI SRA EFetch for accessions, library values, instrument model, runs, file
fallbacks, and read lengths. It uses ENA only for a four-column read-run file
report, preferring a non-empty ENA FASTQ list over NCBI file entries. PubMed
ESummary is driven by PubMed IDs already present in GEO MINiML. It does not
currently fetch the separate BioSample/BioProject fixtures or consume ENA
Browser, search, Xref, Taxonomy, analysis, assembly, or checklist records.

SRA/ENA enrichment does not overwrite GEO geographic or biological fields.
The sequencing SDRF handler explicitly chooses GEO sample-level library values
and instrument model over conflicting SRA values, logs the disagreement, and
keeps the GEO accession when the SRA `geo_sample` differs. ENA precedence is
limited to file metadata for the same run. These observed rules inform future
converter design but do not define precedence for an SRA- or ENA-native input.

The ENA schemas retain their provider `schemaLocation` values. Resolve local
SRA imports through `docs/sra/schemas/` rather than rewriting XSD bytes;
complete `ENA.webin.xsd` compilation additionally needs excluded `EGA.*.xsd`
dependencies. Likewise, NCBI's archive-oriented `EXPERIMENT_PACKAGE_SET`
contains extensions beyond submission-oriented SRA 1.5 schemas, so fixture
verification asserts well-formedness and field contracts rather than claiming
strict whole-document XSD validity.

**Evidence:** [`tests/test_provider_reference_material.py`](../tests/test_provider_reference_material.py),
[`sources/insdc.py`](../src/meta_standards_converter/sources/insdc.py),
[`metadata/enrichment.py`](../src/meta_standards_converter/metadata/enrichment.py),
and [`magetab/sdrf/constructor.py`](../src/meta_standards_converter/magetab/sdrf/constructor.py).

<a id="project-purpose-and-layout"></a>
## Project Purpose And Layout

`meta_standards_converter` converts biological study metadata between repository standards. The current package focuses on GEO MINiML to ArrayExpress/MAGE-TAB-style output.

```text
src/meta_standards_converter/
├── cli/                       # unchanged command names and option contracts
├── converters/                # GEO2JSON, GEO2AE, AE2JSON, JSON2AE/TSV/H5AD/OBS
├── sources/                   # GEO, MAGE-TAB/BioStudies, PubMed, INSDC, JSON groups
├── miniml/                    # typed models, codec, patches, pure GEO XML parser
├── magetab/
│   ├── parser.py              # IDF/SDRF decoding
│   ├── constructor.py         # ordered evidence and build orchestration
│   ├── writer.py              # IDF/SDRF file publication
│   ├── semantics.py           # semantic overlay and round-trip contracts
│   ├── protocols.py           # operation-local registry
│   ├── technology.py          # technology selection
│   ├── idf.py                 # network-free IDF row construction
│   └── sdrf/                  # model, constructor, renderer, technology handlers
├── metadata/                  # enrichment, interpretation, overrides, provenance
│   └── projection/            # tabular, AnnData, assay semantics
├── expression/                # assets, planning, readers, normalization, memory
│                              # references, nfcore, checkpoints, catalogue, components
├── atlas_v1/                  # unchanged standalone Atlas input contract
├── helpers/                   # shared JSON, requests and host gates
├── retrieval.py               # shared source and asset policy
├── runtime_contracts.py       # resource and failure contracts
└── artifact_bundle.py         # durable publication and recovery

```

Tests cover parser packaging, converter orchestration, CLI flags, AE constructor composition, IDF behavior, and SDRF rendering:

```text
tests/miniml/test_geo_parser.py
tests/converters/test_geo2ae.py
tests/converters/test_geo2json.py
tests/converters/test_json2ae.py
tests/converters/test_ae2json.py
tests/sources/test_ae_webfetcher.py
tests/expression/test_json2h5ad.py
tests/cli/test_cli_geo2ae.py
tests/cli/test_cli_geo2json.py
tests/cli/test_cli_json2ae.py
tests/cli/test_cli_ae2json.py
tests/cli/test_cli_json2h5ad.py
tests/test_project_scripts.py
tests/magetab/test_ae_constructor.py
tests/magetab/test_ae_sdrf_handlers.py
tests/metadata/test_miniml_enricher.py
tests/test_request_helper.py
tests/sources/test_geo_webfetcher.py
tests/sources/test_insdc_webfetcher.py
tests/sources/test_pubmed_webfetcher.py
tests/test_retrieval.py
tests/test_runtime_contracts.py
tests/test_xml_safety.py
tests/GSE328265_family.xml
```

<a id="runtime-behavior"></a>
## Runtime Behavior

- Distribution version `6.0.0` keeps typed immutable MINiML packages the Python conversion boundary. It uses
  H5AD metadata schema 1.0 and
  consumes Atlas document schema 1.0 and MINiML ledger schema 1.0;
  neither build metadata nor production imports depend on ThematicAtlases.
- The package requires Python `>=3.10`.
- The 2026-08-02 local compatibility point passed the complete deterministic
  suite on Python 3.12 with python-dateutil 2.9.0.post0, requests 2.34.2,
  AnnData 0.13.2, h5py 3.16.0, NumPy 2.4.6, pandas 3.0.5, Scanpy 1.12.3, and
  SciPy 1.18.0. Declared next-major ceilings contain that tested point; they are
  compatibility bounds, not claims that every intervening version was tested.
- Base runtime dependencies are `requests>=2.31,<3` and
  `python-dateutil>=2.8.2,<3`; the `h5ad` extra bounds AnnData `<1`, Scanpy
  `<2`, NumPy `<3`, pandas `<4`, SciPy `<2`, and h5py `<4` while retaining the
  documented minimum versions.
- `dependency-provenance/pylock.python312-linux-x86_64.toml` locks the complete
  external base/H5AD Python 3.12/Linux x86_64 resolution by exact version,
  artifact URL, and SHA-256. `runtime.python312-linux-x86_64.cdx.json` is the
  deterministic CycloneDX 1.6 inventory derived from that lock.
  `release-policy.json` requires a verified signed artifact manifest and an
  approved offline advisory snapshot no older than seven days, with
  critical/high/medium/low remediation SLAs of 2/7/30/90 days. Those trusted
  operator artifacts are deliberately absent, leaving the composing release
  gate blocked instead of fabricating security evidence.
- The `geo2ae`, `geo2json`, `json2ae`, `ae2json`, `json2h5ad`, `json2tsv`, and `json2obs` console scripts point to their matching modules under `meta_standards_converter.cli`.
- Network calls are owned by platform fetchers and routed through `RateLimitedRequester`: `GEOWebFetcher` handles GEO FTP MINiML tarballs and related-series traversal, `AEWebFetcher` handles BioStudies discovery and HTTP(S) MAGE-TAB text, `INSDCWebfetcher` handles NCBI SRA EFetch plus ENA Portal file reports, and `PubmedWebFetcher` handles NCBI PubMed ESummary publication metadata.
- Default request settings are derived from the standard resource profile and
  enforced across the process by normalized hostname: 10-second connect and
  60-second read timeouts with at most four network workers. Service-specific
  request delays and three bounded retries remain in force. The opt-in large
  profile uses 30/300-second timeouts and eight network workers.
- Library logging propagates safe structured telemetry to caller handlers.
  DEBUG records service/host, attempt, status, timeout, and duration without URL
  queries or request parameters. INFO records retries, GEO fetch sizes/duration,
  MINiML structural counts, related-series progress, and enrichment hit/failure
  totals. XML, parsed metadata, publication content, tokens, and credentials are
  never logged.
- `GEO2AEConverter.convert()` keeps parsed and enriched GEO metadata in memory for MAGE-TAB construction.
- `GEO2JSONConverter.convert()` returns parsed GEO package JSON, enriched by default, and can write `{accession}.json`.
- `JSON2AEConverter.convert()` loads one parsed package object or a non-empty package list, enriches it by default, and returns or writes MAGE-TAB outputs.
- `AE2JSONConverter.convert()` resolves one IDF and one or more SDRFs, returns one MINiML-compatible package in a list, and can write `{accession}.json`.
- `JSON2H5ADConverter.convert()` selects per-sample H5AD, matrix, or raw FASTQ sources; normalizes them into AnnData; and writes a per-sample H5AD catalogue without matrix integration.
- `json2tsv --format {tsv,csv}` emits the neutral MSC sample projection; `json2obs` row-aggregates only sample observation metadata. Both are orchestrator methods with CLI wrappers, and both share the public `MINiMLMetadataProvider` interpretation boundary with `json2h5ad`.
- When `out` is supplied, `GEO2AEConverter.convert()` writes `{accession}.idf.txt` and `{accession}.sdrf.txt`.
- `geo2ae` `out` controls MAGE-TAB output only; use `geo2json` for parsed JSON snapshots.
- Processed `json2h5ad` conversion requires the `h5ad` extra. Raw processing directly on the host additionally requires Nextflow, Java, and a supported execution profile/runtime. The project image includes Java 21, pinned Nextflow, the Docker CLI, and `.[h5ad]`.

<a id="end-to-end-geo2ae-flow"></a>
## End-To-End geo2ae Flow

```text
main(argv)
  -> parse CLI args
  -> instantiate GEO2AEConverter()
  -> for each GSE accession:
       GEO2AEConverter.convert(gse, related_series, remove_empty, out, platform_handler)
       continue to the next accession if a conversion fails
  -> return 1 if any accession failed, else 0

GEO2AEConverter.convert(gse, related_series, remove_empty, out, platform_handler=None)
  -> GEOWebFetcher.fetch_gse_miniml(gse)
  -> GEOParser.parse(miniml, remove_empty=remove_empty, related_series=related_series)
  -> MINiMLEnricher.enrich(data) for each parsed package
  -> instantiate one shared AEConstructor()
  -> for each enriched metadata package:
       AEConstructor.miniml2magetab(data, platform_handler=platform_handler)
  -> if out:
       AEConstructor.magetab2file(magetab, out) for each MAGE-TAB payload
  -> return list of MAGE-TAB payloads
```

`GEO2JSONConverter.convert(gse, related_series, remove_empty, enrich, out)` follows the same GEO fetch and parse stages, optionally enriches each parsed package through `MINiMLEnricher`, writes `{gse}.json` when `out` is truthy, and returns the list of JSON packages without invoking AE/MAGE-TAB construction.

The persisted JSON and H5AD workflows are documented separately under End-To-End json2ae Flow and End-To-End json2h5ad Flow.

External calls in the live conversion path are isolated behind fetchers:

- `GEOWebFetcher.fetch_gse_miniml()` calls GEO FTP through `RateLimitedRequester(service="geo_ftp")`.
- `PubmedWebFetcher.pubmed_summary()` calls NCBI PubMed ESummary through `RateLimitedRequester(service="ncbi_eutils")`.
- `INSDCWebfetcher.fetch_sra_runs()` calls NCBI SRA EFetch and ENA Portal file reports through service-specific `RateLimitedRequester` instances.

CLI behavior:

- Positional `gse` accepts one or more GEO Series accessions.
- `--related`, `--related-series`, and `--get-related-series` are aliases that enable related super/subseries traversal.
- `--remove-empty` is the default and removes empty parsed MINiML fields before conversion.
- `--keep-empty` preserves empty parsed MINiML fields.
- `--out DIR` defaults to the current directory.
- `--platform-handler KEY` forces the same validated technology key through IDF and SDRF generation.
- `--list-platform-handlers` requires no GSE input, prints the stable keys one per line, and exits without conversion.
- Failed accessions log an error and traceback to stdout through the configured logger, then later accessions still run.

<a id="end-to-end-json2ae-flow"></a>
## End-To-End json2ae Flow

```text
main(argv)
  -> parse JSON paths, --out, --no-enrich, platform-handler, and logging flags
  -> instantiate JSON2AEConverter()
  -> for each JSON path:
       JSON2AEConverter.convert(json_path, enrich, out, platform_handler)
       continue to the next path if conversion fails
  -> return 1 if any path failed, else 0

JSON2AEConverter.convert(json_path, out, enrich=True, platform_handler=None)
  -> fail if the path does not exist or JSON decoding fails
  -> normalize one package object to a one-element list
  -> require a non-empty list of package objects
  -> validate each package contains a usable study accession; GSE accessions must be numeric
  -> for each package in order:
       MINiMLEnricher.enrich(data) when enrich=True
       AEConstructor.miniml2magetab(data, platform_handler=platform_handler)
  -> if out:
       AEConstructor.magetab2file(magetab, out) for each payload
  -> return the ordered list of MAGE-TAB payloads
```

Enrichment is enabled by default to match the live `geo2ae` path and may call PubMed, NCBI SRA, and ENA. `--no-enrich` or `enrich=False` makes conversion operate on the supplied JSON without those enrichment calls. The input is otherwise not rewritten. All packages are validated before enrichment or construction begins, and related Series require no separate traversal flag because their packages are already represented in the input list.

The converter uses one injected or default `MINiMLEnricher` and `AEConstructor` per instance. It delegates file naming, IDF/SDRF validation, TSV rendering, and overwrite behavior to `AEConstructor.magetab2file()`. Logs contain paths, package indexes, counts, and stages rather than metadata payloads.

Both MAGE-TAB-producing CLIs expose `--platform-handler KEY` and standalone `--list-platform-handlers`. A forced handler bypasses unchanged source-table reuse and typed-model-only rendering so the chosen IDF/SDRF handlers always run; unsupported MAGE-TAB extension fields are restored afterward where possible. Omitting the option preserves automatic detection and existing round-trip behavior.

<a id="end-to-end-ae2json-flow"></a>
## End-To-End ae2json Flow

```text
main(argv)
  -> parse one or more sources, optional repeated --sdrf, --out,
     --resource-profile/--resource-override, --source-host, and logging flags
  -> require exactly one source when --sdrf is present
  -> for each source:
       AE2JSONConverter.convert(source, out, sdrf_sources)
       continue to the next source if conversion fails
  -> return 1 if any source failed, else 0

ae2json(resource_profile, resource_overrides, source_hosts)
  .convert(source, out=None, sdrf_sources=None)
  -> AEWebFetcher.resolve(source, sdrf_sources)
       existing path: bounded-read the IDF and relative/local/HTTPS SDRFs
       HTTPS URL: validate egress, stream IDF text, and resolve approved SDRFs
       accession: paginate bounded BioStudies JSON, query study info, then
                  validate and stream only the IDF/SDRF metadata files
  -> AEParser.parse(resolved_input)
       parse the IDF and rectangular SDRF tables
       map known investigation, publication, contributor, protocol, sample,
       platform, factor, characteristic, SRA/FASTQ, and array-file metadata
       merge repeated sample/platform records in first-seen order
       keep the first conflicting scalar and record a warning
       set series.iid from the primary ArrayExpress/investigation accession
       normalize factor/material values into the strict MINiML 3.0 vocabulary
       preserve original values in typed fields/characteristics and source digests
  -> if out, write [{package}] to {study_accession}.json
  -> return [package]
```

IDF labels are matched case- and whitespace-insensitively. The parser accepts general MAGE-TAB inputs rather than only files emitted by this project, including MSC's former `Status Term ...` and `Protocol Type Term ...` ontology-companion aliases. Canonical spellings win if both are present. Its public output is immutable MSC MINiML 3.0: repository source format and document name/URI/media type/SHA-256 live under `source`, while typed protocols, variables, assay paths, samples, platforms, and accessions live in their canonical model fields. Raw IDF/SDRF bodies and the former `mage_tab` runtime sidecar are not retained. `series.iid` prefers the explicit ArrayExpress accession and cannot be displaced by a GEO or ENA secondary accession. Values outside the XSD vocabulary are normalized only where required for strict validation, with the original scientific value retained in the adjacent typed value or characteristic rather than discarded.

Accession resolution calls `GET /api/v1/files/{accession}` to discover exactly one IDF and at least one SDRF, calls `GET /api/v1/studies/{accession}/info` for the HTTPS base, and downloads only those metadata files beneath `Files/`. API JSON and MAGE-TAB text are UTF-8/BOM decoded only after streamed declared/actual byte checks. Every URL and redirect is restricted to HTTPS, approved provider or explicit exact hosts, and public DNS answers. Referenced assay data is not downloaded.

<a id="geo2json-vs-ae2json"></a>
## geo2json Versus ae2json

Both converters return a list of package dictionaries using the same MINiML-compatible core vocabulary, but they do not promise identical keys or values. `geo2json` is a generic projection of actual GEO MINiML XML; `ae2json` is an explicit MAGE-TAB projection with a separate typed extension for information the core cannot represent.

| Behavior | `geo2json` | `ae2json` |
|---|---|---|
| Input | One GEO `GSE...` accession | IDF path/URL or BioStudies/ArrayExpress accession plus one or more SDRFs |
| Parsing | Generic XML element/attribute mapping to snake_case | Explicit IDF-row and SDRF-column mappings |
| Packages | One package per selected GEO Series, including related Series when requested | One consolidated package per resolved IDF and its SDRFs |
| Root `version` | MINiML document version such as `0.5.0` | Normalized `magetabv1.1` |
| Root `schema_location` | MINiML XSD location from the XML root | BioStudies MAGE-TAB v1.1 specification URL |
| `series.iid` | GEO Series IID, normally the `GSE...` accession | Primary ArrayExpress accession, with investigation-accession fallback |
| Series/sample fields | Every field present in the selected MINiML elements | Only fields with an explicit exact core mapping |
| Factors | Native MINiML `Variable` records when supplied | IDF experimental factors plus SDRF factor values |
| Protocols | Native GEO protocol elements | Complete typed records in `series.protocols` |
| Platform/contributors | Referenced GEO records with their available MINiML fields | IDF/SDRF projection, usually a smaller record |
| PubMed/SRA enrichment | Enabled by default and optional with `enrich=False`/`--no-enrich` | No remote enrichment stage; publication, ENA, run, and file metadata come from IDF/SDRF fields |
| Empty fields | Removed by default or retained with `remove_empty=False`/`--keep-empty` | Omitted unless a mapped source value exists; required extension structure remains present |
| Assay-row multiplicity | MINiML samples plus optionally enriched `sra_run` lists | Core samples/runs may consolidate rows; every modeled SDRF row remains an independent `series.assay_paths` record |
| Unsupported metadata | Remains available when it exists as an XML element/attribute | Supported protocol/assay extensions become typed steps; source document names, URIs, media types, and digests remain as provenance, while unmodeled raw layout is deliberately dropped |
| Lossless MAGE-TAB round trip | Not applicable | Semantic typed content is reconstructable; exact raw row order/layout requires retaining the original IDF/SDRF documents outside MINiML 3.0 |

Representative GEO output:

```json
{
  "version": "0.5.0",
  "schema_location": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML http://www.ncbi.nlm.nih.gov/geo/info/MINiML.xsd",
  "series": {"iid": "GSE123", "accession": [{"value": "GSE123", "database": "GEO"}]},
  "sample": [],
  "platform": [],
  "contributor": []
}
```

Representative AE output:

```json
{
  "miniml_schema_version": "2.0",
  "source": {
    "format": "MAGE-TAB",
    "documents": [{"kind": "idf", "name": "E-MTAB-1.idf.txt", "sha256": "..."}]
  },
  "series": {
    "iid": "E-MTAB-1",
    "accession": [{"value": "E-MTAB-1", "database": "ArrayExpress"}],
    "protocols": [],
    "assay_paths": []
  },
  "sample": [],
  "platform": [],
  "contributor": []
}
```

The shared core makes downstream processing reusable; it does not imply
field-for-field parity between repositories. Consumers that need exact source
layout must retain the original MAGE-TAB documents; MINiML 3.0 preserves their
identity/digests and the supported typed scientific semantics, not raw tables.

<a id="json2h5ad-flow"></a>
## End-To-End json2h5ad Flow

```text
JSON2H5ADConverter.convert(json_path, out, asset_manifest, asset_specs, force_reprocess, ...)
  -> load ordinary parsed MINiML JSON or validate a canonical Atlas v1 document
  -> group packages by dataset; fail if no convertible groups
  -> one group: convert directly
  -> multiple groups: convert each below out/{dataset_id}, recording group exceptions
  -> AssetManifest loads explicit CSV/TSV and CLI mappings
  -> SourcePlanner invokes injected AssetDiscovery and selects per sample:
       explicit H5AD > explicit matrix > JSON H5AD > JSON matrix > raw FASTQ
  -> if force_reprocess: require raw FASTQ for every sample
  -> NFCoreRunner groups raw samples by detected modality
       -> ReferenceResolver validates catalogue/custom reference and annotation combinations
       -> AnnotationConverter validates local files and converts GFF3 to shared GTF
       -> write nf-core samplesheet and params.json
       -> subprocess.run(nextflow run nf-core/{scrnaseq|rnaseq}, shell=False)
       -> discover scrnaseq H5AD or rnaseq count/TPM matrices
  -> load H5AD (including .h5ad.gz), 10x HDF5, 10x MTX, or delimited matrices
  -> normalize sparse AnnData with canonical dotted msc.* obs fields
  -> retain reversible sample values in uns["msc_metadata"] and provenance in uns
  -> invoke ordered metadata projectors for additive sample obs/var/uns metadata
  -> flatten the permitted MINiML metadata into uns["msc_miniml"]
  -> estimate peak resident memory against fixed and live-availability bounds
       -> oversized: record skip in memory_report and continue
       -> admitted: normalize, checkpoint immediately, release matrix
  -> copy durable checkpoints into one normalized H5AD per sample
  -> do not join expression matrices or invoke combined-study projector callbacks
  -> write JSON catalogue manifest with expression_integration=none
  -> return ConversionResult for one group or BatchConversionResult for multiple
```

`AssetDownloader` streams HTTP(S)/FTP processed assets into an output-local cache and verifies an MD5 when supplied. Source files are never modified. Gzip-compressed H5AD assets are expanded into a temporary `.h5ad` only while AnnData reads them; the cached download remains compressed. Study-level H5AD splitting recognizes canonical `msc.sample.accession` and the external generic columns `geo_accession`, `sample_id`, `sample`, and `gsm_accession`.

<a id="h5ad-metadata-schema-v3"></a>
<a id="h5ad-metadata-schema-v1"></a>
### H5AD metadata schema 1.0

Converter-owned observation metadata uses only dotted names grouped under `msc.sample`, `msc.series`, `msc.platform`, `msc.archive`, `msc.library`, `msc.instrument`, `msc.protocol`, `msc.database`, `msc.asset`, `msc.expression`, `msc.characteristics`, and `msc.observation`. The converter does not generate underscore aliases. Existing underscore-style columns from an input H5AD remain opaque source columns: normalization preserves but neither interprets nor validates them. Custom projectors retain ownership of their injected names.

Stable observation fields cover sample/study accessions, title/description, organism and taxid, organism part, developmental stage, disease, genotype, biological source, material/provider/molecule, platform, SRA/ENA/BioSample/run accessions, library fields, instrument, modality, asset provenance, and database identity. Organism resolution evaluates each channel independently. Database identity uses `public_id`, then `iid`, then `name`. Every raw characteristic becomes `msc.characteristics.<normalized_name>`; typed annotations become `msc.characteristics.harmonized_<field>` plus identifier and ontology companions. Sample-bound assay-node material types take precedence over the typed `material_type` characteristic, legacy channel material, normalized molecule, and organism-part fallbacks. Native assay parameters become `msc.assay.parameter.<name>.*`, with their occurrence ledger in `uns["msc_assay"]`. Missing values are empty in `obs`; repeated values are case-insensitively de-duplicated in source order and displayed with `; ` separators.

`uns["msc_metadata"]` declares schema version `1.0` and contains the authoritative normalized `sample_values` DataFrame with `sample_accession`, `field`, `ordinal`, `value`, and `value_type`. It stores one row per non-empty canonical value, so embedded semicolons and list cardinality remain recoverable without parsing the display string. Each sample H5AD contains only its sample rows. `uns["msc_miniml"]` remains the complete typed source ledger at schema 1.0. H5AD provenance and manifests separately declare the H5AD metadata schema version; the catalogue manifest additionally declares `artifact_kind = per_sample_h5ad_catalogue`, `expression_integration = none`, and a non-verified combination state.

Normalization copies each incoming index into `msc.observation.original_id`. An identifier is already sample-qualified when its accession occurs case-insensitively as a token bounded by the start/end or `-`, `_`, `.`, or `:`. Qualified identifiers are preserved; other identifiers receive `-{sample_accession}`. Repeated candidates receive source-order numeric suffixes. A sequential used-ID ledger qualifies any remaining collision before each sample checkpoint is written, so global uniqueness does not require retaining earlier expression matrices.

`uns["msc_miniml"]` contains schema/policy metadata, the source JSON path and SHA-256, and a typed long-form `fields` DataFrame (`package_index`, `entity_type`, `entity_id`, `path`, `value`, `value_type`). GSM files contain the sample plus its series and transitively referenced platform, contributor, and database records without following `sample_ref`; both scalar references and real MINiML `{"ref": "..."}` objects are resolved. GSE files contain all package entities. Protocol descriptions remain in this table. `msc.protocol.types`, source refs, and accessions preserve exact typed ontology values in sample-bound assay-path order; when no path binds the sample, all declared protocols retain study order. `Harmonizer.geoprotocols2efo()` is only the compatibility fallback when no applicable typed protocol exists. Publication records are whitelisted to PubMed ID, DOI, title, authors, status, and status ontology fields; abstracts, full text, article bodies, sections, and other publication content are not embedded. GEO series summary and overall design remain experiment metadata.

Persisted local provenance paths in the manifest and H5AD metadata are relative to the containing artifact and declare `path_base = artifact_parent`; parallel scope fields distinguish internal, external (`../...`), and remote locations. Remote URLs remain unchanged. Recorded Nextflow command path arguments are also relative, while generated runtime configs retain the absolute paths required by Nextflow resume. `ConversionResult` continues to return absolute paths in memory. ANSI-stripped, de-duplicated Nextflow warnings are stored on each pipeline run and promoted to the manifest's top-level warnings without changing a successful return code.

Before writing nf-core samplesheets, `NFCoreRunner` upgrades raw `ftp://` URLs from the known ENA and NCBI archive hosts to their equivalent `https://` endpoints. This avoids truncated Java FTP transfers through rootless container networking while leaving unknown FTP servers unchanged.

Raw processing pins `nf-core/scrnaseq` 4.2.0 and `nf-core/rnaseq` 3.26.0 by default. The scrnaseq 4.2.0 floor includes the upstream strict-syntax fixes required by the pinned Nextflow 26 runtime; 4.1.0 contains a reference to a missing `conf/test_multiome.config` and fails during config parsing. The runner requires Nextflow, Java, and the selected Docker/Podman/Apptainer/Singularity runtime. `scrnaseq` discovers H5ADs across the results tree, associates outputs when the sample accession appears in the filename or containing directories, and prefers CellBender-filtered, then filtered, then raw output. This includes QCATCH names such as `GSM1_filtered_quants.h5ad`, not only `*_matrix.h5ad`. For each sample, `rnaseq` selects its named numeric column from the merged count and TPM tables, tolerates the standard text `gene_name` column, writes counts to sparse `X`, writes aligned TPM values to `layers["tpm"]`, and preserves `gene_name` in `var`.

When `META_STANDARDS_REQUIRE_ROOTLESS_DOCKER` is truthy and the Docker profile is selected, `NFCoreRunner._preflight()` queries `docker info` before creating workflow files. An unreachable daemon or security options without `rootless` abort the conversion before Nextflow starts. Other deployments retain the existing runtime-presence checks.

Catalogue publication preserves successful per-sample outputs regardless of
expression modality, organism, declared reference build, or feature namespace;
those differences are descriptive heterogeneity, not an integration claim.
The compatibility helper treats every all-unknown dimension as missing,
classifies all-numeric gene IDs as Entrez rather than symbols, and returns
unknown for mixed/ambiguous namespaces. Its `combine()` method fails with
guidance to use an explicit scientific integration workflow.

`JSONPackageSource` delegates canonical documents to `AtlasV1Reader`. The
reader retains only `harmonized` datasets, preserves canonical dataset IDs,
and reports other states with diagnostics; v1 envelopes are rejected. Exact
`metadata.packages` wrappers expand within the dataset group before shared
sample deduplication.
Native MINiML grouping still deduplicates identical samples and rejects
conflicting duplicates. `JSON2H5ADConverter.convert_source()` runs every group separately;
multi-study output uses one child directory per study and returns
`BatchConversionResult`. The existing `convert()` contract remains a
`ConversionResult` for one study and returns the batch result only for a
multi-study source. Source warnings and per-study failures remain auditable
without discarding successful studies.

GEO 10x `.matrix.mtx`, `.barcodes.tsv`, and `.genes.tsv`/`.features.tsv`
companions are grouped by the standard `SourcePlanner`. The reader localizes
all members into a temporary Scanpy-compatible directory, including legacy
gzip-compressed genes trios.

<a id="json2tabular-flow"></a>
## End-To-End json2tsv Manifest Flow

```text
JSON2TSVConverter.export_manifest(source, outdir, output_format)
  -> AtlasV1Reader/JSONPackageSource loads MINiML or harmonized Atlas v1 metadata
  -> group packages by study and visit every sample in source order
  -> build TabularMetadataContext with normalized MINiML sample metadata
  -> invoke ordered TabularMetadataProjector objects
  -> reject projector column collisions
  -> fail closed on projection errors unless allow_invalid=True
  -> write preferred columns first and remaining columns sorted
  -> return TabularConversionResult with rows, datasets, warnings, and errors
```

With no explicit projectors, `MSCMetadataProjector` emits stable dotted
`msc.sample`, `msc.series`, `msc.platform`, `msc.archive`, `msc.library`,
`msc.instrument`, `msc.protocol`, `msc.database`, and `msc.expression`
columns followed by sorted `msc.characteristics.*` columns. Supplying an
explicit projector list replaces that default contract, allowing private
schemas to own the complete table without leaking into this repository.

<a id="reference-annotation-flow"></a>
### Reference And Annotation Flow

Raw nf-core processing accepts these reference combinations:

```text
--genome KEY
--genome KEY --gtf annotation.gtf[.gz]
--genome KEY --gff annotation.gff|gff3[.gz]
--fasta genome.fa[.gz] --gtf annotation.gtf[.gz]
--fasta genome.fa[.gz] --gff annotation.gff|gff3[.gz]
```

`ReferenceResolver.resolve()` rejects simultaneous GTF/GFF, annotation without a genome or FASTA, and FASTA without an annotation. Explicit catalogue keys bypass organism inference. Otherwise, confirmed inference supports human/taxid 9606 as `GRCh38` and mouse/taxid 10090 as `GRCm39`.

`AnnotationConverter.prepare()` requires local readable reference files, resolves them to absolute paths, computes the source annotation SHA-256, and passes GTF through unchanged. GFF3 is converted with `gffread SOURCE -T -o OUTPUT` into `nfcore/{study}/reference/{source_sha256}.gtf`. Existing non-empty checksum-addressed output is reused across mixed `scrnaseq`/`rnaseq` runs and resume attempts. Conversion failures remove the temporary output and stop before Nextflow.

Generated nf-core parameters include `genome` plus the explicit/effective `gtf`, or fully custom `fasta` plus `gtf`. Generated `input`, `outdir`, and reference fields override the same keys from a user params file. Annotation source path, input format, SHA-256, and effective GTF path are stored on nf-core assets and pipeline runs, then written to per-sample H5AD `uns["meta_standards_converter"]` and the JSON provenance manifest. The converter validates files and formats but cannot prove FASTA/annotation assembly or chromosome-name compatibility.

<a id="rootless-json2h5ad-runtime"></a>
## Rootless json2h5ad Runtime

The deterministic suite was refreshed on 2026-08-10 and reported
`587 passed, 3 skipped` (plus 89 unittest subtests). The public wire contract is Atlas document schema 1.0
and converter output uses H5AD metadata schema 1.0.

`Dockerfile` builds the application image with Python 3.12, Java 21, Nextflow 26.04.2 verified by SHA-256, Docker CLI 29.6.2, `gffread`, and the H5AD extra. It contains no Docker daemon.

`scripts/provision-rootless-json2h5ad.sh` is the administrative boundary. It installs rootless prerequisites, creates the locked `nfcore-runner` account, allocates a non-overlapping 65,536-ID subordinate range, enables its user service, and configures ACLs. The build context is read-only to the runner; `.out/json2h5ad` is the only writable project path.

`scripts/json2h5ad-compose.sh` must run as `nfcore-runner`. It resolves that user's socket, with `ROOTLESS_DOCKER_SOCKET` as the only override seam, refuses a daemon without the `rootless` security option, prepares output-local home and Nextflow caches, and invokes `compose.yaml` with absolute paths. Compose passes the same absolute `JSON2H5AD_OUT` into the container that it uses for the working directory and bind mount. The application process uses container UID/GID 0:0, which the rootless daemon maps to the unprivileged host `nfcore-runner` identity; using UID 1001 inside the container would instead map to a subordinate host UID without output or socket access.

```text
administrator provisions nfcore-runner once
  -> rootless dockerd listens at /run/user/<uid>/docker.sock
  -> runner helper validates docker info SecurityOptions
  -> rootless Compose creates the converter container
       -> only .out/json2h5ad and the rootless socket are mounted
       -> json2h5ad preflight independently verifies rootless mode
       -> Nextflow asks the same unprivileged daemon for per-process containers
  -> nf-core tasks and converter outputs remain owned by the runner namespace
  -> H5AD mode 0660 keeps the provisioned project-user ACL effective
```

The Compose container drops all capabilities, enables `no-new-privileges`, uses a read-only root filesystem, and receives a `noexec` tmpfs `/tmp`. Nextflow alone is pointed through `NXF_OPTS` at a separate executable 2 GiB `/nextflow-tmp`; this is required because its AWS/S3 client extracts a native library before staging iGenomes references. Host and container output paths are identical because the sibling nf-core task containers must bind the Nextflow work files by host-visible absolute path. The system rootful socket is never mounted. Final H5AD files use mode `0660`: this preserves the output directory's named project-user ACL while denying access to other users. Provisioning and Compose verify effective runner/project-owner ACL access; filesystems that map the rootless writer to `nobody:nogroup` are accepted when those checks pass, and no privileged ownership repair is attempted.

<a id="rootless-acceptance-2026-07-31"></a>
### Rootless acceptance evidence — 2026-07-31

The dedicated acceptance run on 2026-07-31 completed `nf-core/rnaseq` 3.26.0
and `nf-core/scrnaseq` 4.2.0 with return code 0 and non-partial H5AD results.
See [`rootless-acceptance-2026-07-31.md`](rootless-acceptance-2026-07-31.md).

<a id="parsed-miniml-data-shape"></a>
## Parsed MINiML Data Shape

`GEOParser.parse()` returns `list[MINiMLPackage]`, with one canonical,
self-contained package per top-level MINiML `Series`. `AEParser.parse()`
returns one `MINiMLPackage` in the same canonical representation. Package
objects implement the read-only mapping interface used by legacy callers.

```python
[
    {
        "miniml_schema_version": "2.0",
        "source": {"format": str, "version": str | None, "documents": list[dict]},
        "database": list[dict],
        "organization": list[dict],
        "contributor": list[dict],
        "platform": list[dict],
        "sample": list[dict],
        "series": dict,
    }
]
```

`AEParser.parse()` uses the same core package vocabulary but identifies its
source dialect and documents under `source`. Its `series.iid` is the explicit
ArrayExpress accession, then an ArrayExpress-form investigation accession,
then another ArrayExpress-classified accession, with the investigation
accession as fallback. GEO secondary accessions remain in `series.accession`
and do not displace an available ArrayExpress IID.

`miniml_schema_version` versions MSC's JSON representation. It is independent
of `source.version`, which records the source MINiML or MAGE-TAB dialect.

Top-level package keys are singular. Parser keys inside each parsed XML element are original XML names converted to snake_case. Repeated XML elements also keep the singular snake_case key and point to a list.

Examples:

```text
Sample-Ref          -> sample_ref
Pubmed-ID           -> pubmed_id
Data-Table          -> data_table
Supplementary-Data  -> supplementary_data
Raw-Data            -> raw_data
```

Element text mapping is generic:

- Plain leaf elements with no attributes or child elements become strings.
- Elements with attributes become dictionaries containing those attributes.
- When an attribute-bearing or child-bearing element also has text, the text is stored as `value`.
- Values remain strings; the parser does not coerce dates, numbers, booleans, or ontology identifiers.
- Namespaces are stripped to local names.
- Non-`version` root attributes are attached to each package under snake_case keys.

Example:

```xml
<Characteristics>whole larval tissue</Characteristics>
<Characteristics tag="time">30 Days</Characteristics>
```

parses as:

```python
{"characteristics": ["whole larval tissue", {"tag": "time", "value": "30 Days"}]}
```

`MINiMLEnricher` adds remote lookup results without changing the raw parsed GEO fields:

- `series.pubmed_publication`: one dict per `series.pubmed_id`, with `pubmed_id`, `doi`, `author_list`, `title`, `status`, `status_term_source_ref`, and `status_term_accession_number`.
- `sample.*.sra_accession`: SRA/ENA/DDBJ accessions extracted from SRA sample relations.
- `sample.*.ena_accession`: deduplicated study/project accessions such as `ERP137216` collected from SRA run enrichment.
- `sample.*.sra_run`: run dicts returned by `INSDCWebfetcher.fetch_sra_runs()`, including study accession, library metadata, run/sample IDs, read lengths, instrument model, and per-FASTQ `filename`/`uri`/`md5`.

<a id="miniml-package-model"></a>
## MSC MINiML 3.0 package API

MSC owns the unified metadata representation in
`meta_standards_converter.miniml`. The dependency-free Python model and its
codec are derived from the repository's MINiML XSD. The Python model is the
sole structural authority; MSC does not publish or package a parallel JSON
Schema. The model covers the package entities (`Database`, `Organization`, `Contributor`,
`Platform`, `Sample`, and `Series`) and reusable accession, reference, status,
person, channel, table/data, variable, repeat, organism, relation, and link
structures. Protocols, protocol applications, assay nodes and paths, ontology
values, named values, occurrence-local harmonized values, comments, and source documents
fold MAGE-TAB semantics into the same immutable representation.
PubMed publications, SRA/ENA accessions, SRA runs, and FASTQ file records are
first-class typed enrichments.

The wire discriminator is `miniml_schema_version: "3.0"`. Runtime decoding
rejects unversioned, 1.x, and unknown versions. MINiML 2.0 remains readable and
is deterministically migrated by `MINiMLV2Migrator`; `MINiMLV1Migrator` and the
`miniml-migrate` command provide explicit one-way migration to 3.0. Canonical
collections are always lists, while `series` remains a single object.
`MINiMLCodec.decode`/`decode_many` return immutable packages plus structured
compatibility diagnostics. Strict mode promotes blocking diagnostics to
`MINiMLCompatibilityError`; policy `miniml-3.0-source-compat-v1` accepts only
sample-title `xsd_uniqueness` warnings when every sample sharing that title has
a non-empty unique iid. Accepted warnings remain in the decode result, titles
remain byte-for-byte source faithful, and duplicate/missing identities plus all
other diagnostics remain blocking. The public
`MINIML_STRICT_COMPATIBILITY_POLICY_VERSION` constant identifies this policy
for cache and runtime provenance. `encode`/`encode_many` are the serializer boundary,
while `load`/`dump` provide deterministic, atomic UTF-8 JSON publication.
`MINiMLPackage.from_mapping()` and `load()` remain direct model conveniences.
Construction and codec encode/decode canonicalize and validate nested typed
objects; direct construction recursively freezes mapping and collection inputs,
and duplicate top-level identifiers or malformed entity shapes raise
`MINiMLModelError`.

`MINiMLPackage.validate()` is compatibility-first. XSD vocabulary deviations,
checksum formats, unresolved references, and channel-count mismatches are
reported as structured `MINiMLValidationIssue` diagnostics instead of rejecting
historically accepted data. External protocol/sample references are therefore
preserved as warnings. Wire validation is performed by `MINiMLCodec` and
`MINiMLPackage.from_mapping()`; there is no secondary schema contract for
consumers to reconcile.

Harmonized evidence is stored beside the raw occurrence, never in an
`annotations` array. A group uses `hz_<field>`, optional `hz_<field>_id`,
`hz_<field>_onto`, and `hz_<field>_hierarchy_depth`; collisions use aligned
`(1)`, `(2)`, ... suffixes. Named characteristic lists store adjacent named
rows, while ontology/value objects store the same keys as members. Every
companion requires a value. `HarmonizedValue`, `iter_harmonized_values`, and
the mapping helpers are the shared validated projector interface.
`append_harmonized_value(destination, value, *, name_key=None)` is the shared
writer: it validates the existing destination, reuses an exact semantic
identity, allocates the next free collision index for a distinct value, and
writes aligned mapping members or `name`/`tag` rows without replacing raw
evidence.

MSC 5.2 adds the generic `MINiMLHarmonizationPatch` 3.1 contract without
changing the MINiML 3.0 discriminator. An operation contains an occurrence
path, a typed harmonized value and optional bounded source evidence. The
`exact_value` and `exact_span` claims are Unicode-normalized and revalidated at
their source value pointer; `interpreted` evidence remains auditable but does
not claim literal authorship. `apply_miniml_harmonization_patch` writes or
deduplicates the occurrence-local `hz_*` group, partitions document patches
into package-local fragments, rebases their pointers, and attaches them in
application order under `extensions.msc_harmonization`. Fingerprints omit only
that retained extension, so unrelated extensions remain identity-bearing and
reapplication is idempotent. Package-list paths always carry a leading package index, including lists with
exactly one package; partitioning therefore follows the document shape rather
than treating a one-package list as a package mapping.
`iter_harmonization_patches`,
`iter_harmonization_operations`, and `harmonization_provenance_index` provide
one validated downstream view after one canonical package decode; internal
iteration reuses that mapping instead of recursively decoding it. Canonical
package encoding also hoists every
legacy `series.extensions` entry to package `extensions`; identical duplicates
deduplicate, conflicts fail closed, and the reserved `msc_harmonization` key
cannot enter through the legacy series surface.

Generic converters expose this evidence rather than using it as an implicit
replacement policy. TSV/CSV emits deterministic indexed
`msc.harmonization.<field>.*` columns, MAGE-TAB adds machine-readable adjacent
`Comment[msc_harmonization_*]` columns without synthesizing a new
`Characteristics[...]`, and H5AD/obs carries the fragment ledger in
`uns["msc_harmonization"]` while the self-contained package remains in
`uns["msc_miniml"]`. The opt-in `harmonization_overrides` policy remains a
separate contract.

`NamedValue` carries a
typed ontology value, optional unit ontology, `unit_type`, and qualifier;
`Variable.type` is an ontology value so factor type source/accession companions
round-trip without string flattening. Source documents retain document role,
source URI, media type, and SHA-256 of the consumed UTF-8 content (never the raw
body), plus document-scoped order. Series fields retain experiment design ontology,
experiment date, contacts and roles, and generic IDF comments.

Public symbols are exported from `meta_standards_converter.miniml`; neither a
schema-path helper nor JSON Schema package data is public. Contract coverage
lives in `tests/test_msc_miniml_v2.py`,
`tests/test_miniml_migration_cli.py`, `tests/magetab/test_magetab_miniml_v2.py`, and
`tests/miniml/test_geo_parser.py`. Cross-boundary stabilization coverage lives in
`tests/test_miniml_stabilization.py`.

The complete qualified model API is
`meta_standards_converter.miniml.model.Accession`,
`meta_standards_converter.miniml.model.Address`,
`meta_standards_converter.miniml.model.Channel`,
`meta_standards_converter.miniml.model.Characteristics`,
`meta_standards_converter.miniml.model.Contributor`,
`meta_standards_converter.miniml.model.DataColumn`,
`meta_standards_converter.miniml.model.Database`,
`meta_standards_converter.miniml.model.DataTable`,
`meta_standards_converter.miniml.model.FASTQFile`,
`meta_standards_converter.miniml.model.InstrumentModel`,
`meta_standards_converter.miniml.model.MINiMLModelError`,
`meta_standards_converter.miniml.model.MINiMLPackage`,
`meta_standards_converter.miniml.model.MINiMLValidationIssue`,
`meta_standards_converter.miniml.model.Organization`,
`meta_standards_converter.miniml.model.Organism`,
`meta_standards_converter.miniml.model.Person`,
`meta_standards_converter.miniml.model.Platform`,
`meta_standards_converter.miniml.model.PubMedPublication`,
`meta_standards_converter.miniml.model.Reference`,
`meta_standards_converter.miniml.model.Relation`,
`meta_standards_converter.miniml.model.Repeat`,
`meta_standards_converter.miniml.model.Sample`,
`meta_standards_converter.miniml.model.Series`,
`meta_standards_converter.miniml.model.SRARun`,
`meta_standards_converter.miniml.model.Status`,
`meta_standards_converter.miniml.model.SupplementLink`,
`meta_standards_converter.miniml.model.TableData`,
`meta_standards_converter.miniml.model.Variable`,
`meta_standards_converter.miniml.codec.MINiMLBatchDecodeResult`,
`meta_standards_converter.miniml.codec.MINiMLCodec`,
`meta_standards_converter.miniml.codec.MINiMLCompatibilityError`, and
`meta_standards_converter.miniml.codec.MINiMLDecodeResult`,
`meta_standards_converter.miniml.harmonization.HarmonizedValue`,
`meta_standards_converter.miniml.harmonization.append_harmonized_value`,
`meta_standards_converter.miniml.harmonization.harmonized_mapping`,
`meta_standards_converter.miniml.harmonization.harmonized_value_mappings`,
`meta_standards_converter.miniml.harmonization.is_harmonized_key`,
`meta_standards_converter.miniml.harmonization.iter_harmonized_values`,
`meta_standards_converter.miniml.harmonization.named_harmonized_rows`,
`meta_standards_converter.miniml.harmonization.next_harmonized_index`,
`meta_standards_converter.miniml.harmonization.parse_harmonized_key`, and
`meta_standards_converter.miniml.harmonization.parse_harmonized_mapping`,
`meta_standards_converter.miniml.patches.MINiMLHarmonizationPatch`,
`meta_standards_converter.miniml.patches.canonical_miniml_document`,
`meta_standards_converter.miniml.patches.miniml_source_fingerprint`,
`meta_standards_converter.miniml.patches.apply_miniml_harmonization_patch`,
`meta_standards_converter.miniml.patches.validate_harmonization_extension_mapping`,
`meta_standards_converter.miniml.patches.iter_harmonization_patches`,
`meta_standards_converter.miniml.patches.iter_harmonization_operations`,
`meta_standards_converter.miniml.patches.harmonization_provenance_index`, and
`meta_standards_converter.metadata.provenance.patch_provenance_columns`.

The typed semantic additions are
`meta_standards_converter.miniml.model.NamedComment`,
`meta_standards_converter.miniml.model.NamedValue`,
`meta_standards_converter.miniml.model.OntologyValue`,
`meta_standards_converter.miniml.model.SourceDocument`,
`meta_standards_converter.miniml.model.SourceInfo`,
`meta_standards_converter.miniml.model.Protocol`,
`meta_standards_converter.miniml.model.ProtocolApplication`,
`meta_standards_converter.miniml.model.AssayNode`,
`meta_standards_converter.miniml.model.AssayPath`,
`meta_standards_converter.miniml.migration.MINiMLMigrationResult`,
`meta_standards_converter.miniml.migration.MINiMLV1Migrator`,
`meta_standards_converter.miniml.migration.MINiMLV2Migrator`,
`meta_standards_converter.magetab.semantics.overlay_miniml_semantics`, and
`meta_standards_converter.magetab.semantics.render_miniml_assay_documents`, and
`meta_standards_converter.cli.miniml_migrate.main`.

MAGE-TAB parsing folds its parser state immediately into native protocols,
declarations, assay paths, attributes, units, comments, and document
provenance. No `mage_tab`, raw table, layout, row/column index, ordinal, or
synthetic object id survives. Construction regenerates ordered IDF rows and
repeated SDRF columns from the model, so the supported round trip is semantic.

**Evidence:** [`model.py`](../src/meta_standards_converter/miniml/model.py),
[`codec.py`](../src/meta_standards_converter/miniml/codec.py),
[`geo_parser.py`](../src/meta_standards_converter/miniml/geo_parser.py),
and [`ae_parser.py`](../src/meta_standards_converter/magetab/parser.py).

<a id="workflow-details"></a>
## Workflow Details

<a id="geo-parse-flow"></a>
### GEO Parse Flow

```text
GEOParser.parse(miniml, remove_empty, related_series)
  -> _parse(miniml)
       -> ET.fromstring(miniml)
       -> _top_level_nodes(root)
       -> _parse_element(each top-level node)
       -> _build_indexes(parsed_top_level)
       -> _series_package(root, each series, indexes)
  -> if related_series:
       _parse_with_related_series(parsed)
  -> if remove_empty:
       remove_empty_fields(each package)
  -> return parsed package list
```

Per-Series packages include only records relevant to that series:

- Samples referenced by `series.sample_ref[*].ref`.
- Platforms referenced by included sample `platform_ref.ref`.
- Contributors referenced by series, sample, or platform `contributor_ref` and `contact_ref`.
- Organizations referenced by included contributors or databases through `organization_ref.ref`.
- Databases referenced by included `accession[*].database` or `status[*].database`.

Missing references are tolerated. The original reference remains in place, and the unresolved target record is omitted from package lists.

The parser uses an XSD-inspired `repeated_children` map for known repeated MINiML fields. Unknown repeated sibling tags still become lists if they occur more than once.

<a id="related-series-flow"></a>
### Related-Series Flow

```text
parse(miniml, related_series=True)
  -> _parse input MINiML into root packages
  -> _parse_with_related_series(root packages)
       -> seed seen_gses from series.accession[*].value
       -> seed queue from superseries/subseries relation entries
       -> fetch unseen related GSE MINiML
       -> _parse related MINiML
       -> append related packages
       -> enqueue newly discovered related GSEs
       -> stop when queue is empty
  -> return root packages plus related packages
```

`GEOParser.parse_related_series(miniml, remove_empty=False, strict=True)` uses
the same traversal logic but returns only related packages, excluding the input
packages, in a list-compatible `RelatedSeriesParseResult`. The result always
contains status 2.0 plus attempted/failed accessions. When `strict=True`, fetch
or parse failures raise. When `strict=False`, successful packages are retained
with degraded/partial/review-required status and persistence-safe error
envelopes; raw provider exception messages are neither returned nor logged.

Related GSE accessions are discovered from `series.relation` entries only when relation `type`, `target`, or `comment` mentions superseries/subseries and contains `GSE` followed by digits.

<a id="idf-and-mage-tab-construction-flow"></a>
### IDF And MAGE-TAB Construction Flow

`AEConstructor` is a coordinator, not an SDRF subclass. It owns an `IDFConstructor` and an `SDRFConstructor`, supplied as optional dependencies or created by default.

```text
AEConstructor.miniml2magetab(data)
  -> create one ProtocolRegistry for the series accession
  -> detect one shared AE technology key
  -> SDRFConstructor._miniml2sdrf(data, protocol_registry, technology_type)
       -> generate SDRF table
       -> register actual non-empty Protocol REF values and required placeholder refs
  -> IDFConstructor.miniml2idf(data, protocol_registry, technology_type)
       -> build IDF rows
       -> emit Protocol rows from the same registry
       -> include empty ["SDRF File"] placeholder after protocol rows
       -> infer term source rows
       -> move Experiment Description after Investigation Title
       -> move top-level Comment[...] rows to the bottom
  -> replace existing SDRF File placeholder with ["SDRF File", sdrf_table]
  -> return MAGE-TAB payload
```

`ProtocolRegistry` normalizes protocol text, reuses refs for identical `(kind, text)` pairs, and names refs as `P-{series_accession}-{n}`. Known protocol kinds map to MAGE-TAB labels such as `Extract-Protocol`, `Hybridization-Protocol`, `Scan-Protocol`, and `Data-Processing`. Required placeholder refs can be created with empty protocol text for protocols that must be present in IDF and SDRF.

`AEConstructor.magetab2file()` normalizes legacy mixed payloads, finds the SDRF row, validates the SDRF payload is a non-empty row table, chooses filenames from `Comment[ArrayExpressAccession]`, `Investigation Accession`, or `Comment[SecondaryAccession]`, writes IDF/SDRF TSV files, and returns the IDF path.

<a id="sdrf-graph-and-rendering-flow"></a>
### SDRF Graph And Rendering Flow

The SDRF code models each rendered row as an `SDRFPath`, an ordered list of nodes and protocol edges.

- `SDRFAttr`: companion column label, value, nested attrs, and required flag.
- `SDRFNode`: primary SDRF column such as `Source Name`, `Extract Name`, `Assay Name`, `Scan Name`, or file columns.
- `SDRFEdge`: `Protocol REF` column with optional attrs.
- `SDRFPath`: one rendered row path.
- `ColumnGroup`: planned primary column plus companion columns.
- `SDRFAudit`: warnings, dropped values, and validation errors.

```text
SDRFConstructor._miniml2sdrf(data, protocol_registry, technology_type=None)
  -> use supplied technology_type, or fall back to _detect_sdrf_technology(data)
  -> select handler class
  -> handler.build()
       -> build_paths()
       -> plan_columns(paths)
       -> render_paths(columns, paths)
  -> store handler.audit as self.last_sdrf_audit
  -> return SDRF table rows
```

Column planning scans every path before rendering rows. Repeated visible labels are disambiguated internally with occurrence keys such as `Array Data File#1`; visible labels can repeat when MAGE-TAB expects repeated columns.

<a id="technology-handler-selection"></a>
### Technology Handler Selection

`AEConstructor._detect_ae_technology()` reads platform technology/title, library source/strategy/selection, sample type, series text, sample text, channel text, file extensions, and SRA relations. `SDRFConstructor._detect_sdrf_technology()` remains available for direct SDRF callers, but delegates to the AE constructor detector.

Selectable handler keys are defined once by `PLATFORM_HANDLER_KEYS`:

```text
plate_single_cell_sequencing
droplet_single_cell_sequencing
tenx_v2_droplet_single_cell_sequencing
tenx_v3_droplet_single_cell_sequencing
single_cell_sequencing
spatial_sequencing
bulk_sequencing
sequencing
array
generic
```

`single_cell_sequencing` and `sequencing` are supported manual dispatch keys even though the current detector normally returns a more specific single-cell variant or `bulk_sequencing`. Unknown keys are rejected by the CLIs and direct `AEConstructor.miniml2magetab()` calls. `--list-platform-handlers` prints this catalog in the displayed order.

The README Guide contains the user-facing configuration reference and a Mermaid graph of this catalog's conceptual specialization. That graph intentionally describes selectable platform relationships rather than the literal private IDF and SDRF class trees, which differ in their 10x implementations.

Detection rules in broad order:

- High-throughput sequencing platform text, SRA relations, library strategy, or sample type `SRA` choose bulk sequencing.
- Single-cell text chooses single-cell sequencing variants.
- `visium` or `spatial` text takes spatial precedence.
- `10x`, `droplet`, or `chromium` text chooses droplet single-cell.
- Within 10x/Chromium droplet context, standalone `v2` and `v3` markers choose the 10x v2/v3 droplet keys; generic 10x text remains `droplet_single_cell_sequencing`.
- Single-cell text without droplet hints chooses plate single-cell.
- Array platform technology or array-like files choose array.
- Everything else uses the generic handler.

Assay terms such as ChIP-seq, ATAC-seq, multiome, methylation, or array assay names do not by themselves create special technology handlers beyond sequencing or array.

<a id="sequencing-sdrf-flow"></a>
### Sequencing SDRF Flow

Sequencing paths use:

```text
Source Name
Protocol REFs for treatment/growth/extraction
Extract Name
Protocol REF for library construction
Assay Name
optional Protocol REF for `sample.scan_protocol`
Scan Name
Factor Value[...] columns
```

Sequencing behavior:

- SRA accessions are extracted from sample relations and cached per handler.
- SRA records contribute library layout/source/strategy/selection, run/experiment/sample accessions, submitted file names, MD5, instrument model, read lengths, and FASTQ metadata.
- GEO values take precedence over conflicting SRA values, and warnings are recorded in `SDRFAudit`.
- `Comment[LIBRARY_SOURCE]` is rendered uppercase in the SDRF only; library source text used for protocol descriptions keeps its original cleaned casing.
- Library construction is distinct from extraction: its synthesized description uses only available layout/strategy/source/selection metadata and its IDF type maps to `nucleic acid library construction protocol` (`EFO_0004184`).
- A populated `sample.scan_protocol` is registered and referenced before the sequencing scan node.
- Shared sequencing paths render FASTQs as `Comment[readN file]`, `Comment[FASTQ_URI]`, and `Comment[MD5]` companion columns.
- Bulk and plate sequencing paths emit one SDRF row per FASTQ URI, with duplicate sample/source metadata allowed across rows.
- More than two FASTQs are preserved; shared sequencing uses additional read comments, while bulk and plate sequencing use additional rows.
- When no SRA FASTQs exist, sequencing raw files from GEO supplementary/raw data can become read comments.
- Sequencing handlers emit supplementary processed assets as `Derived Array Data File` columns while raw reads remain FASTQ comments.
- Single-cell handlers add library construction, technical replicate, read geometry, isolation, or spatial read-index comments where their subclass supports it.

<a id="array-sdrf-flow"></a>
### Array SDRF Flow

Array paths use:

```text
Source Name
Protocol REFs for treatment/growth/extraction
Extract Name
optional labeling Protocol REF
optional Labeled Extract Name
Hybridization Protocol REF
Assay Name
optional Scan Protocol REF
Scan Name
file nodes
Factor Value[...] columns
```

Array behavior:

- Two-channel samples emit one path per channel and warn through the audit object.
- Channel `label` creates a `Labeled Extract Name` with a `Label` companion.
- `Array Design REF` is taken from the sample platform accession and has nested `Term Source REF = ArrayExpress`.
- `.tif` and `.tiff` files render as `Image File`.
- `.cel`, `.gpr`, `.idat`, `.exp`, `.rpt`, `.cab`, and similar raw array files render as `Array Data File`.
- Matrix-like supplementary files render as `Derived Array Data Matrix File`; other processed assets render as `Derived Array Data File`.
- Repeated raw and derived files are preserved as repeated columns and recorded as warnings.

<a id="base-sdrf-behavior"></a>
### Base SDRF Behavior

Base handler behavior shared by generic, sequencing, and array handlers:

- Sample ordering follows `series.sample_ref` first, then any remaining samples.
- A sample with multiple channels emits multiple channel paths and records a warning.
- `Source Name` uses the GSM accession when available, falling back to sample `iid`.
- GEO channel source is preserved as `Comment[Sample_source_name]`.
- Sample title and description render as mapped sample comment columns when present.
- `Characteristics[organism]`, `Characteristics[organism part]`, `Characteristics[developmental stage]`, `Characteristics[disease]`, and `Characteristics[genotype]` are always planned for source nodes and may be blank.
- `Characteristics[organism]` reads each channel organism's legacy `name` field first, then falls back to the parsed MINiML `value` field produced from `<Organism taxid="...">text</Organism>`.
- Explicit `organism part` wins; otherwise `tissue` wins; otherwise channel `source` is used as fallback.
- Repeated characteristics with the same tag are preserved as repeated columns and recorded as warnings.
- Provider and material type are mapped from channel biomaterial provider and molecule/source context.
- Factor values come from series variables when present, otherwise from characteristic tags with more than one value.
- Required blank protocol refs are emitted as blank cells with audit warnings but do not create IDF protocol rows.
- `normalized_extension()` strips archive suffixes such as `.gz`, `.zip`, `.bz2`, and `.xz` before classifying the underlying extension.

Legacy greedy GEO and SRA fallback comment classes are kept only as commented reference code at the bottom of `ae_sdrf_handlers.py`. They have no runtime effect.

<a id="sra-pubmed-and-ontology-enrichment"></a>
### SRA, PubMed, And Ontology Enrichment

- Normal `GEO2AEConverter.convert()` enrichment happens after `GEOParser.parse()` and before `AEConstructor.miniml2magetab()`.
- `MINiMLEnricher` writes PubMed metadata to `series.pubmed_publication` and SRA metadata to `sample.sra_accession`/`sample.sra_run`.
- IDF construction prefers `series.pubmed_publication`; `_lookup_pubmed_id()` remains as a compatibility fallback.
- SDRF construction prefers `sample.sra_run`; relation-based `_lookup_sra()` remains as a compatibility fallback.
- GEO protocol labels map through `Harmonizer().geoprotocols2efo()`.
- PubMed status maps through `Harmonizer().pubstatus2efo()`.
- Term source names combine non-empty `source ref` cells with every declared `database` record. A matching declared record supplies its URL and exact version, including an intentionally missing version; `Harmonizer` is used only when no matching database exists.

<a id="public-api-and-callable-reference"></a>
## Public API And Callable Reference

This section lists public and semi-public callables used by tests or by package orchestration. Many helper methods are intentionally private but documented here because this project currently relies on direct helper behavior in tests and internal composition.

<a id="cli"></a>
### `cli/geo2ae.py`, `cli/geo2json.py`, `cli/json2ae.py`, `cli/ae2json.py`, `cli/json2h5ad.py`, `cli/json2tsv.py`, and `cli/json2obs.py`

`_parser() -> argparse.ArgumentParser`

- `geo2ae` and `geo2json` build command-line parsers for one or more GSE accessions.
- Adds `--related`, `--related-series`, and `--get-related-series` aliases.
- Adds mutually exclusive `--remove-empty` and `--keep-empty` options.
- Adds `--out`, defaulting to `"."`.
- Adds logging controls: repeatable `-v`/`--verbose`, `-q`/`--quiet`, and `--log-file`.
- `geo2ae` and `json2ae` add mutually exclusive `--platform-handler` and `--list-platform-handlers`; list mode runs without positional inputs or converter construction.
- `geo2json` also adds `--no-enrich`, which skips PubMed/SRA enrichment and writes parsed-only JSON.
- `json2ae` accepts one or more parsed MINiML or canonical Atlas v1 JSON paths, adds `--no-enrich`, and writes IDF/SDRF files under `--out`.
- `ae2json` accepts one or more IDF paths, HTTP(S) IDF URLs, or BioStudies accessions. Repeatable `--sdrf` overrides are allowed with exactly one source.
- `json2h5ad` accepts parsed JSON plus `--asset`/`--asset-manifest`, source and matrix controls, catalogue or user FASTA references, `--gtf`/`--gff` annotation overrides, pinned nf-core execution controls, `--resume`, `--processed-checkpoint-dir`, `--overwrite`, and `--allow-invalid`. Resume covers both Nextflow work and fingerprint-valid processed-sample checkpoints.
- `json2tsv` accepts parsed MINiML or Atlas JSON and writes the sample manifest as TSV by default or CSV with `--format csv`; `--out` selects an exact file and `--outdir` derives a filename.
- `json2obs` accepts the same metadata and raw/processed data inputs and processed-checkpoint controls as `json2h5ad`, writes a required row-aggregated `obs.csv` with an explicit `cell_id` column without combining expression matrices, and can add typed `uns.json` or single-sample `var.csv` sidecars.

`main(argv=None) -> int`

- Parses arguments, configures `meta_standards_converter` logging to stdout and optional file output, creates the command's converter, and converts each accession or JSON file in order.
- Default logging emits `WARNING+`; `-v` emits `INFO+`, `-vv` emits `DEBUG+`, and `--quiet` emits `ERROR+`.
- On conversion failure, logs the exception traceback, marks the run failed, and continues.
- Success/progress messages are logged rather than printed; normal success output appears with `-v`.
- Returns `1` if any accession failed, otherwise `0`.

<a id="converter"></a>
### `converters/geo2ae.py`, `converters/geo2json.py`, `converters/json2ae.py`, `converters/ae2json.py`, `converters/json2h5ad.py`, and `converters/json2tsv.py`

`class geo2ae(JSONHandler)`

- Main programmatic converter.
- The class inherits `JSONHandler`, though the converter path does not currently rely on inherited helper methods.
- `__init__(enricher=None, geo_fetcher=None, parser=None, ae_constructor=None)` accepts enrichment, GEO fetcher, parser, and MAGE-TAB constructor dependencies. Defaults are `MINiMLEnricher()`, `GEOWebFetcher()`, `GEOParser(geo_fetcher=self.geo_fetcher)`, and `AEConstructor()`.

`convert(gse, related_series=False, remove_empty=True, out=None, platform_handler=None)`

- Fetches MINiML with `self.geo_fetcher.fetch_gse_miniml(gse=gse)`.
- Parses with `self.parser.parse(miniml, remove_empty=remove_empty, related_series=related_series)`.
- Enriches each parsed package with `self.enricher.enrich(data=meta_json)`.
- Reuses the injected or default `AEConstructor`.
- Converts each enriched package to a MAGE-TAB payload.
- Passes a non-`None` `platform_handler` through to `AEConstructor.miniml2magetab()`.
- Writes each payload when `out` is truthy.
- Returns the list of MAGE-TAB payloads.

`class geo2json(JSONHandler)`

- Main programmatic GEO-to-JSON converter.
- `__init__(enricher=None, geo_fetcher=None, parser=None)` accepts enrichment, GEO fetcher, and parser dependencies with the same defaults as `geo2ae`.

`convert(gse, related_series=False, remove_empty=True, enrich=True, out=None)`

- Fetches and parses MINiML using the same GEO fetcher/parser path as `geo2ae`.
- Enriches each parsed package by default; `enrich=False` returns parsed-only GEO JSON.
- Writes one `{gse}.json` file containing the full package list when `out` is truthy.
- Returns the list of JSON packages.

`class json2ae(JSONHandler)`

- `__init__(enricher=None, ae_constructor=None, package_source=None)` accepts injectable enrichment, MAGE-TAB construction, and JSON package-source collaborators.
- `convert(json_path, out=None, enrich=True, platform_handler=None) -> list[list]` accepts the package form written by `geo2json`, one package object, or a canonical Atlas v1 document.
- Atlas loading retains harmonized dataset metadata, logs shared-reader warnings for skipped states, rejects v1 envelopes, and raises `ValueError` when no convertible groups remain.
- Validates the entire top-level shape, package types, and study accessions before invoking collaborators. Non-GEO accessions are accepted; `GSE...` values retain numeric validation.
- Enriches packages by default; `enrich=False` preserves the supplied metadata and avoids enrichment calls.
- Passes a non-`None` `platform_handler` through to `AEConstructor.miniml2magetab()`.
- Builds all MAGE-TAB payloads in input order, then writes each through the shared `AEConstructor` when `out` is truthy.
- Returns the ordered in-memory MAGE-TAB payload list.

`class ae2json`

- `__init__(fetcher=None, parser=None)` accepts injectable `AEWebFetcher` and `AEParser` collaborators.
- `convert(source, out=None, sdrf_sources=None) -> list[dict]` resolves and parses one MAGE-TAB study and returns a one-package list.
- When `out` is supplied, writes the package list to `{ArrayExpress-or-first-accession}.json`, with path separators made filename-safe.

`class JSON2H5ADConverter`; compatibility alias `class json2h5ad`

- Accepts injectable `SourcePlanner`, `NFCoreRunner`, `AssetDownloader`,
  `DatasetCombinationPolicy`, and ordered `AnnDataMetadataProjector`
  collaborators, plus available-memory and peak-memory-estimator test seams.
- `convert(..., allow_invalid=False, resume=False, force_memory=False) -> ConversionResult | BatchConversionResult` accepts ordinary
  parsed MINiML JSON or a canonical Atlas v1 document. It returns the
  single-group result directly and aggregates multiple groups.
- `convert_source(json_path, out=None, allow_invalid=False, **options) -> BatchConversionResult`
  always aggregates groups. Per-group exceptions populate `failures` and do
  not discard successful conversions.
- `ConversionResult` exposes the compatibility field `combined_h5ad` (always
  `None` for catalogue conversion), `sample_h5ads`, retained pipeline files,
  pipeline commands, warnings/errors/failures, first-sample `primary_h5ad`, and
  `partial`; `memory_report` records every newly assessed sample admission or
  skip and is also persisted in the catalogue manifest.
- `AssetManifest` loads CSV/TSV mappings or `ACCESSION=PATH` CLI specifications. Manifest entries outrank CLI entries, which outrank discovered JSON assets.
- `AssetManifest.load(path: str) -> list[Asset]` reads CSV/TSV, requires
  `scope_id`/`path`, groups raw members, and raises `ValueError` for blank or
  unsupported records. `parse_spec(spec: str) -> Asset` parses the compact
  CLI form and raises on malformed specifications; neither writes files.
- `AssetDownloader.localize(value: str, md5: str | None = None) -> str`
  returns local paths unchanged or streams HTTP(S)/FTP into its cache, verifies
  an optional digest, and raises on transport/checksum failure; downloading is
  its filesystem/network side effect.
- `AssetDownloader.retention_report(*, max_age_seconds,
  min_retained_assets=1, active_paths=(), now=None, apply=False) -> dict`
  returns a dry-run-first integrity/age/source report. Applying a plan moves
  only verified, old, non-active, non-minimum assets and their sidecars into a
  timestamped recoverable quarantine; it never deletes cache data.
- `SourcePlanner.plan(packages, explicit_assets=None, force_reprocess=False)
  -> dict[str, Asset]` selects one asset per sample; `discover(packages) ->
  list[Asset]`, `samples(packages) -> list[str]`,
  `sample_accession(sample) -> str | None`, and
  `classify(path) -> str | None` expose discovery/classification without
  filesystem writes. Planning raises when required raw coverage is absent.
- `NFCoreRunner.process(assets, packages, out, study_accession,
  pipeline="auto", genome=None, fasta=None, gtf=None, gff=None,
  accept_inferred_reference=False, profile="docker", revision=None,
  params_file=None, nextflow_config=None, work_dir=None, resume=False)
  -> RawProcessingResult` validates runtime/reference inputs, writes workflow
  inputs/logs, invokes Nextflow without a shell, discovers outputs, and raises
  on invalid configuration, runtime preflight, subprocess, or output failure.
- `ReferenceResolver` accepts a catalogue `genome` with an optional GTF/GFF override or `fasta` paired with exactly one GTF/GFF; supported organism inference must be explicitly accepted before Nextflow starts.
- `AnnotationConverter` validates local FASTA/annotation paths, records annotation SHA-256, passes GTF through, and converts GFF3 to a shared checksum-addressed GTF through `gffread`.
- Generic delimited matrices require an explicit orientation when it cannot be represented by a study-scoped sample column.
- `_scientific_modules()` loads AnnData, NumPy, pandas, and SciPy only;
  `_scanpy_module()` is called exclusively by 10x HDF5/MTX readers, keeping
  ordinary processed-H5AD conversion independent of Scanpy import side effects.

`MetadataProjectionContext`, `AnnDataMetadataProjection`, and the
`AnnDataMetadataProjector` protocol form the metadata extension
contract. Sample projectors run after standard `_normalize()` processing and
before MINiML attachment/writing. Former `project_combined` callbacks are not
part of the protocol and are not invoked. Scalar
`obs`/`var` values broadcast, vector values must match their axis, and existing
axis or top-level `uns` keys cannot be overwritten. Observation drops and
renames are checked as one operation and applied before additions; invalid
sources, duplicate targets, collisions, and rename/drop overlap fail closed.
Projector warnings are
deduplicated into `ConversionResult.warnings` and the manifest. Reported
`errors` raise `AnnDataProjectionError` before publication unless
`allow_invalid=True`, which records them and returns a partial result.
Structural return-type, collision, and vector-length errors remain
unconditional. All H5ADs and the manifest are staged before a backup/swap
commit; commit failure restores the prior complete bundle. With no projectors,
output is unchanged.

`JSON2TSVConverter`

- Accept injectable ordered `TabularMetadataProjector` collaborators.
- Use `MSCMetadataProjector` only when no explicit projector list is supplied.
- Accept both MINiML package JSON and canonical Atlas v1 JSON through
  `JSONPackageSource`.
- Select TSV or CSV serialization with the validated `output_format` argument.
- Return `TabularConversionResult`; projection errors fail closed unless
  `allow_invalid=True`.

<a id="miniml-enricher"></a>
### `metadata/enrichment.py`

`class MINiMLEnricher`

`__init__(pubmed_fetcher=None, insdc_fetcher=None)`

- Accepts PubMed and INSDC fetcher dependencies for tests or custom network behavior.
- Defaults to `PubmedWebFetcher()` and `INSDCWebfetcher()`.

`enrich(data: dict) -> dict`

- Mutates and returns one parsed MINiML package.
- Calls `enrich_pubmed()` and `enrich_sra()`.

`enrich_pubmed(data: dict) -> dict`

- Deduplicates `series.pubmed_id` values while preserving first-seen order.
- Adds `series.pubmed_publication` records with the fields consumed by `IDFConstructor`.
- On request or XML parse errors, records the PubMed ID with `None` metadata values so late IDF rendering does not retry the lookup.

`enrich_sra(data: dict) -> dict`

- Extracts SRA accessions from `sample.relation` entries whose type is `SRA`.
- Adds `sample.sra_accession`, `sample.sra_run`, and `sample.ena_accession` when fetched runs contain study accessions.
- On request or XML parse errors, keeps the accession and leaves that accession's run contribution empty.

<a id="geo-web-fetcher"></a>
### `sources/geo.py`

`class GEOWebFetcher`

`__init__(requester=None, request_settings=None)`

- Defaults to `RateLimitedRequester(service="geo_ftp")`.
- Accepts a custom requester or GEO request settings for tests and advanced callers.

`url_gse_miniml(gse: str) -> str`

- Requires an accession starting with `GSE`, case-insensitive.
- Converts accessions to the GEO FTP bucket pattern. For example, `GSE234602` becomes bucket `GSE234nnn`.
- Returns the GEO FTP HTTPS URL ending in `{gse}_family.xml.tgz`.

`fetch_gse_miniml(gse) -> str`

- Builds the URL with `url_gse_miniml()`.
- Downloads the `.tgz` archive through the `geo_ftp` requester and calls `raise_for_status()`.
- Streams every tar member without filesystem extraction, accepting safe
  auxiliary regular files/directories but requiring exactly one root
  `{gse}_family.xml`. Absolute/traversal paths, duplicate names, links and
  special files, unexpected XML, excessive member counts, and excessive total
  expanded size fail closed. The bounded expected XML is returned as UTF-8.

<a id="ae-web-fetcher"></a>
### `sources/magetab.py`

`TextResource(name, text, origin)` and `MAGETabInput(idf, sdrfs, source, source_kind)` are immutable transport records used between resolution and parsing.

`class AEWebFetcher`

- `__init__(requester=None, request_settings=None)` defaults to `RateLimitedRequester(service="biostudies")`.
- `resolve(source, sdrf_sources=None) -> MAGETabInput` dispatches existing paths, HTTP(S) URLs, and accession tokens. Missing path-like inputs raise `FileNotFoundError` instead of becoming accession lookups.
- Local and HTTP IDFs use explicit SDRF overrides when supplied; otherwise every `SDRF File` value is resolved relative to the IDF.
- Accession lookup requires exactly one discovered IDF and at least one SDRF. Explicit overrides are rejected for accession sources.
- Remote metadata is fetched as text and never persisted by the fetcher.

<a id="ae-parser"></a>
### `magetab/parser.py`

`class AEParser`

- `parse(source: MAGETabInput) -> dict` parses one IDF plus all SDRFs into the existing MINiML-compatible package shape.
- Root metadata uses MINiML 3.0 with normalized MAGE-TAB format/version and specification provenance under `source`. `series.iid` prefers `Comment[ArrayExpressAccession]`, then an ArrayExpress-form investigation/classified accession, then the investigation accession.
- IDF rows are normalized by case and whitespace. Repeated row values remain ordered and feed investigation, accessions, design/factor, status, publication, contributor, database, and protocol records.
- SDRF headers map source/sample identities, characteristics, factors, protocol refs, platforms, technology, SRA/ENA runs, FASTQ metadata, and array raw/derived files. Repeated sample rows merge without duplicating list values.
- Conflicting scalar values keep the first value and append a warning. Unknown IDF rows are diagnosed but not retained as raw layout; supported generic SDRF values become typed assay-path steps.
- `build_model(...)` creates an internal semantic bridge containing complete protocol columns, QC/replicate/normalization declarations, every SDRF assay path, ordered nodes and protocol references, comments/files, and per-value unit/ontology companions. The migrator folds it into canonical `series.protocols`, variables, and assay paths before the package leaves `parse()`.
- Malformed non-rectangular SDRF rows fail with a filename and column-count error.

<a id="ae-roundtrip"></a>
### Retired raw MAGE-TAB round-trip sidecar

Raw-table round-trip helpers are retired in MSC 4. The disconnected
`ae_handlers/ae_roundtrip.py` module represented the pre-MINiML-2.0
`mage_tab.roundtrip` sidecar and was neither emitted nor consumed by the active
conversion path. Keeping it public would falsely imply exact IDF/SDRF layout
survives the canonical boundary.

`AEParser` now maps supported scientific content into typed MINiML 3.0
protocols, assay paths, values, units, annotations, and source-document
provenance. `MINiMLV1Migrator` explicitly reports `source_layout_dropped` when
it encounters an old raw-table sidecar, and `AEConstructor` reconstructs
semantic MAGE-TAB from the typed model. Operators requiring byte/layout-exact
round trips must retain the original IDF/SDRF source files identified by the
package's source-document records.

<a id="typed-mage-tab-model"></a>
### `magetab/semantics.py`

- `MAGETabModelError` is the public validation failure and
  `validate_model(model)` enforces schema version 1 collections, unique SDRF
  and assay identities, references, step shapes, and scalar harmonization
  annotations.
- `build_model(idf_rows, sdrfs)` creates the version-1 internal semantic bridge consumed immediately by `MINiMLV1Migrator`; it is not a public wire extension.
- `protocols` contains one position-stable record per IDF protocol, including name, arbitrary type, ontology, description, hardware, software, parameters, contact, and performer.
- `declarations` independently stores aligned quality-control, replicate, and normalization terms with source/accession annotations.
- `assay_paths` contains one record per original SDRF data row. Ordered steps distinguish material/assay nodes, protocol references, annotated characteristics/factors/parameters, comments, files, and generic fields. This preserves array assay multiplicity and many-to-one sample relationships.
- Attribute steps keep `Unit`, `Term Source REF`, and `Term Accession Number` as independent fields; barcode/read geometry remains independent comment steps rather than being folded into protocol prose.
- `render_model(model)` regenerates one SDRF directly or consolidates multiple SDRFs by header plus occurrence. `overlay_core(model_rows, core_rows)` unions eligible core fields into that rendering while retaining model-only protocols, identities, annotations, rows, and structural graph columns.
- IDF matching uses normalized row labels and inserts only rows in the mapped allowlist. Legacy MSC publication/protocol companion labels normalize to the four canonical MAGE-TAB 1.1 rows before rendering. SDRF matching uses `(normalized header, occurrence)` keys, so repeated characteristics remain position-stable. Missing core columns are inserted relative to the nearest core-order neighbor; independent curator fields such as `Characteristics[hz_cell_type]`, `Characteristics[hz_cell_type_id]`, and `Characteristics[hz_cell_type_onto]` remain separate rather than being reinterpreted as native ontology companions.
- SDRF values align through the available `Sample Name`, `Source Name`, and `Comment[ENA_RUN]` identities. A core value replaces or populates a model cell only when all matching core rows agree on exactly one value. An unmatched model row keeps its existing value; an ambiguous newly inserted cell remains blank. Core-only rows are not added or broadcast as new assay paths.
- `MINiMLV1Migrator` folds this bridge into the canonical v2 package and drops the internal container. `AEConstructor` renders from those canonical protocol and assay-path fields; no raw-table fingerprint or replay sidecar participates.

<a id="geo-parser"></a>
### `miniml/geo_parser.py`

<a id="geoparser-class-and-parse-methods"></a>
#### GEOParser class and parse methods

`class GEOParser`

- Owns `repeated_children`, an XSD-inspired map of repeated MINiML children by parent tag.
- Parses only recognized top-level package categories: organization, contributor, database, platform, sample, and series.
- `__init__(geo_fetcher=None)` accepts a GEO fetcher dependency for related-series traversal and defaults to `GEOWebFetcher()`.

`parse(miniml, remove_empty=False, related_series=False) -> list[dict]`

- Parses the input MINiML into per-series packages.
- Optionally traverses related super/subseries.
- Optionally removes empty fields after all parsing/traversal.

`_parse(miniml) -> list[dict]`

- Parses one MINiML XML string without related-series fetching or cleanup.
- Builds top-level parsed records, indexes them by `iid`, and creates one package per series.

`parse_related_series(miniml, remove_empty=False, strict=True) -> RelatedSeriesParseResult`

- Parses the input MINiML, seeds a queue from related-series relations, and returns only fetched related packages.
- Deduplicates GSE accessions.
- Raises on fetch/parse failures in strict mode; non-strict mode retains
  successes with a degraded status 2.0 envelope, attempted/failed accessions,
  and safe errors.
- Applies empty cleanup to related packages when requested.

`RelatedSeriesParseResult` subclasses `list[dict]` for compatibility and adds
`status`, `attempted_accessions`, `failed_accessions`, and `summary_dict()`.

`remove_empty_fields(data)`

- Public wrapper around `_remove_empty_fields()`.
- Removes `None`, empty strings, empty lists, and empty dicts recursively.

<a id="parser-reference-resolution"></a>
#### Reference resolution

- `_top_level_nodes(root)` collects known top-level MINiML elements.
- `_build_indexes(parsed_top_level)` creates `iid` lookup maps.
- `_series_package(root, series, indexes)` assembles a scoped package and attaches root attributes.
- `_resolve_samples()`, `_resolve_platforms()`, `_resolve_contributors()`, `_resolve_databases()`, and `_resolve_organizations()` resolve package records from references.
- `_items_for_refs()` preserves first-seen order and deduplicates refs.
- `_reference_values()` walks nested dicts for ref-bearing keys.

<a id="parser-generic-xml-mapping"></a>
#### Generic XML mapping

- `_parse_element(node)` converts XML recursively to strings, dicts, and lists.
- `_child_key(parent_name, child_name)` currently returns singular snake_case child names.
- `_normalized_text(text)` collapses whitespace.
- `_local_name(tag)` strips XML namespaces.
- `_to_snake_case(value)` normalizes tag/attribute names.

<a id="parser-related-series-helpers"></a>
#### Related-series helpers

- `_extract_series_accessions()` returns normalized `GSE` accessions from package series accessions.
- `_extract_related_gse_accessions()` extracts related `GSE` accessions from relation type/target/comment text.
- `_is_related_series_relation()` recognizes superseries/subseries relation text.

<a id="parser-cleanup-and-helpers"></a>
#### Cleanup and helpers

- `_remove_empty_fields(value)` recursively removes empty values.
- `_is_empty_value(value)` defines empty values as `None`, `""`, `[]`, or `{}`.
- `_walk_dicts(value)` recursively yields nested dicts.
- `_as_list(value)` normalizes scalars and `None` to list handling.
- `_attach_namespaced_root_attributes(root, package)` copies non-version root attributes into packages.

<a id="ae-idf-handlers"></a>
### `magetab/idf.py`

`class IDFConstructor`

`__init__(pubmed_fetcher=None)`

- Accepts a PubMed fetcher dependency for tests or custom network behavior.
- Defaults to `PubmedWebFetcher()`.

`miniml2idf(data, protocol_registry=None, technology_type=None) -> list`

- Builds IDF rows in this order before final normalization: MAGE-TAB version, investigation rows, experimental design/factor rows, person rows, date rows, publication rows, experiment description, protocol rows, `SDRF File` placeholder, term source rows, then platform-specific rows. The `_idf_qc_rep_norm()` extension call is intentionally commented out, so QC/replicate/normalization placeholder rows are not included in final IDF output.
- Final normalization moves the first `Experiment Description` row immediately after the first `Investigation Title` row when both are present.
- Final normalization moves every top-level row whose first cell starts with `Comment[` to the bottom while preserving relative order among non-comment rows and among comment rows. SDRF companion columns such as `Comment[FASTQ_URI]` are unaffected because they live inside the SDRF table, not top-level IDF rows.
- `_move_experiment_description_after_title(rows)` performs the experiment-description relocation, `_is_comment_row(row)` identifies top-level comment rows, and `_move_comment_rows_to_bottom(rows)` performs the stable comment partition.

Investigation and experimental rows:

- `_idf_investigations()` uses series title, series accessions, enriched sample study accessions, converted ArrayExpress-style accessions, and related super/subseries GSE accessions.
- `_secondary_accession_pairs()` emits exactly one `Comment[SecondaryAccession]` row and one positionally aligned `Comment[SecondaryAccessionTermSourceRef]` row. Series accessions come first, followed by first-seen `sample.*.ena_accession` values; duplicates are removed case-insensitively while preserving the first rendered value.
- Secondary-accession term sources are inferred by prefix as `GSE -> GEO`, `ERP -> ENA`, `SRP -> SRA`, and `DRP -> DRA`. Unknown Series prefixes retain their declared database, while unknown enriched sample prefixes receive a blank source-ref cell.
- `_to_arrayexpress_accessions()` replaces `GSE` with `E-GEOD-`.
- `Comment[RelatedExperiment]` is emitted when `series.relation` contains superseries/subseries relation text and related `GSE...` accessions. This row records parsed relationships and does not depend on fetching related packages with `--related`.
- `_idf_experimental()` derives experimental factor names from sample channel characteristics whose normalized tag has more than one distinct normalized value; `Experimental Factor Type` currently mirrors the factor names while factor term source/accession rows remain blank.
- `_idf_platform_specific(data, technology_type)` dispatches to a private platform IDF handler.
- `_idf_term_source(magetab, data)` emits all declared databases plus referenced term sources. Declared URL/version values, including blanks, are source-authoritative; only undeclared sources fall back to `Harmonizer` metadata.
- `_idf_persons(data)` prefixes a real structured address with its organization, preserves string addresses verbatim, and leaves a missing address blank rather than copying the affiliation.

Platform IDF handler inheritance mirrors the SDRF platform tree:

```text
_BasePlatformIDFHandler
├── _SequencingPlatformIDFHandler
│   ├── _BulkSequencingPlatformIDFHandler
│   │   └── _PlateSingleCellSequencingPlatformIDFHandler
│   └── _SingleCellSequencingPlatformIDFHandler
│       ├── _DropletSingleCellSequencingPlatformIDFHandler
│       └── _SpatialSequencingPlatformIDFHandler
├── _ArrayPlatformIDFHandler
└── _GenericPlatformIDFHandler
```

The dispatch keys mostly match SDRF: `plate_single_cell_sequencing`, `droplet_single_cell_sequencing`, `tenx_v2_droplet_single_cell_sequencing`, `tenx_v3_droplet_single_cell_sequencing`, `single_cell_sequencing`, `spatial_sequencing`, `bulk_sequencing`, `sequencing`, and `array`; unknown or missing keys use `_GenericPlatformIDFHandler`. The 10x v2/v3 keys intentionally route to `_DropletSingleCellSequencingPlatformIDFHandler`.

Sequencing platform IDF handlers emit empty label-only `Comment[AEExperiment]`, `Comment[AEExperimentType]`, and `Comment[AECurator]` rows. The `single_cell_sequencing`, `droplet_single_cell_sequencing`, 10x v2/v3 droplet, and `plate_single_cell_sequencing` handlers replace the `Comment[AEExperimentType]` row with `RNA-seq of coding RNA from single cells`; other sequencing handlers, including `spatial_sequencing`, leave it blank. Secondary accessions are already consolidated by `_idf_investigations()`, so platform handlers do not emit another secondary-accession row. Sequencing handlers emit `Comment[SequenceDataURI]` from enriched `sample.*.sra_run[*].run` accessions when valid runs exist. Runs are deduplicated, grouped by prefix such as `ERR` or `SRR`, sorted numerically, and collapsed to one ENA data/view URL per prefix group using min-max ranges, for example `http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041`. Missing and malformed run accessions are logged as warnings and skipped. Droplet single-cell IDF handlers append empty label-only `Comment[AEExpectedClusters]`, `Comment[AEAdditionalAttributes]`, and `Comment[AEBatchEffect]` rows. Array, generic, and unknown handlers do not emit these rows.

Person, QC, and date rows:

- `_idf_persons()` uses contributors for names, email, phone, fax, organization-based affiliation, and flattened addresses prefixed with organization when available.
- `_idf_qc_rep_norm()` still returns label-only placeholder rows for quality control, replicate, and normalization, but `miniml2idf()` does not currently include them because its extension call is commented out.
- `_idf_dates()` normalizes parseable status dates to `YYYY-MM-DD`.
- `_normalized_idf_date()` preserves unparseable values and empty strings.
- `Public Release Date` uses the earliest parseable normalized GEO release date; `Comment[GEOReleaseDate]` preserves all GEO release date values.
- `Comment[ArrayExpressSubmissionDate]` uses the current conversion date as one `YYYY-MM-DD` value while GEO update dates are comments.

Publication, experiment, protocol, and term source rows:

- `_idf_publications()` prefers `series.pubmed_publication`; if absent, it reads `series.pubmed_id`, enriches each ID through PubMed ESummary, and maps publication status.
- `_idf_publications()` logs a warning when neither enriched publication records nor usable `series.pubmed_id` values are present, because PubMed-backed publication rows cannot be populated.
- `_lookup_pubmed_id()` delegates to `self.pubmed_fetcher.pubmed_summary()` and returns DOI, author string, title, status term, source ref, and accession.
- `_idf_experiments()` combines series summary and overall design into `Experiment Description`.
- `_idf_protocols()` uses a supplied `ProtocolRegistry` when present; otherwise it falls back to scanning known GEO protocol paths.
- `_idf_protocols_from_registry()` emits protocols registered by the SDRF build and appends missing required protocol definitions.
- `Protocol Parameters` and `Protocol Contact` row definitions are intentionally commented out; `Protocol Hardware` and `Protocol Software` remain and are extended when required protocol placeholders are appended.
- Required protocol definitions are deduped by harmonized protocol type: `sample collection protocol` is required for all IDFs, while `nucleic acid sequencing protocol` is required for sequencing IDFs.
- `_idf_term_source()` scans rows containing `source ref` and emits source name/file/version rows from `Harmonizer().ontologies`.

Current caveats:

- PubMed lookup is fetcher-owned and normally invoked by `MINiMLEnricher`; direct IDF construction can still invoke it as a fallback when enriched records are absent.
- Unknown PubMed statuses are preserved as literal publication-status labels with blank ontology refs.
- Some date rows are intentionally blank placeholders for internal curation.
- `_idf_term_source()` runs before comment rows are moved, so current term-source inference is based on the pre-normalized IDF rows.

<a id="ae-constructor"></a>
### `magetab/constructor.py`

`class ProtocolRegistry`

- Maps protocol kind/text pairs to stable `P-{series_accession}-{n}` refs.
- Reuses the same ref for identical cleaned text under the same kind.
- Tracks kind, MAGE-TAB label, and cleaned description.
- `get_ref(kind: str, text: str | None, label: str | None = None) -> str |
  None` cleans text, returns `None` when it is empty, reuses existing identity,
  or mutates registry state by allocating the next reference.
- `ensure_required(kind, label)` creates or reuses a required placeholder ref even when protocol text is empty.
- `records()` returns records in insertion order.

`class AEConstructor`

`__init__(idf_constructor=None, sdrf_constructor=None)`

- Accepts dependency injection for tests or custom constructors.
- Defaults to `IDFConstructor()` and `SDRFConstructor()`.

`miniml2magetab(data, platform_handler=None) -> list`

- Creates a `ProtocolRegistry` from `_series_accession(data)`.
- Validates a supplied platform-handler key against `PLATFORM_HANDLER_KEYS`, or detects the shared MAGE-TAB technology key with `_detect_ae_technology(data)` when omitted.
- Forced mode skips unchanged round-trip and typed-model-only shortcuts so both selected handlers execute.
- Builds SDRF first so protocol refs are registered.
- Builds IDF with the same registry and technology key.
- Replaces the first `SDRF File` row with the in-memory SDRF table.
- Preserves valid quote and apostrophe characters in in-memory IDF/SDRF values.
- Raises `ValueError` if the IDF lacks an SDRF row.

`magetab2file(magetab, out=None) -> str`

- Creates the output directory.
- Normalizes row shapes with `_normalize_magetab_rows()`.
- Uses `csv.writer` TSV escaping so quote characters survive written IDF/SDRF files.
- Validates and extracts the embedded SDRF table.
- Chooses IDF and SDRF filenames from the MAGE-TAB accession rows.
- Replaces the embedded SDRF table with the SDRF filename in the IDF.
- Writes both files as tab-delimited UTF-8 text and returns the IDF path.

Other helpers:

- `_detect_ae_technology()` chooses `bulk_sequencing`, `plate_single_cell_sequencing`, `droplet_single_cell_sequencing`, `spatial_sequencing`, `array`, or `generic`. Version-specific option names remain accepted explicitly, but automatic selection does not interpret bare versions.
- `_has_array_files()` detects array-like files from platform/sample/series supplementary data and raw data.
- `_normalize_magetab_rows()` accepts row lists, comma-delimited legacy strings, and legacy `"SDRF file", sdrf` pairs.
- Legacy `_strip_quotes()` remains private but is no longer used by construction or file writing.
- `_sdrf_row_index()` finds the SDRF row case-insensitively.
- `_magetab_accession()` searches ArrayExpress, investigation, then secondary accession rows.
- `_safe_filename_token()` removes path separators from accession-derived filenames.
- `_is_table()` validates row-table shape.
- `_write_tsv()` writes `None` as blank cells.

<a id="sdrf-handlers"></a>
### `magetab/sdrf/constructor.py`

<a id="sdrf-dataclasses"></a>
#### SDRF dataclasses

- `SDRFAttr`: companion attributes, including nested companion attributes.
- `SDRFNode`: visible primary SDRF node columns.
- `SDRFEdge`: visible `Protocol REF` columns.
- `SDRFPath`: one logical row path.
- `ColumnGroup`: planned column plus companions.
- `SDRFAudit`: warnings, dropped values, validation errors.

<a id="sdrfconstructor"></a>
#### SDRFConstructor

`class SDRFConstructor`

`__init__(insdc_fetcher=None)`

- Accepts an INSDC fetcher dependency for tests or custom network behavior.
- Defaults to `INSDCWebfetcher()`.

- `_add_sdrf_to_idf()` appends an in-memory SDRF row to IDF rows; this remains for compatibility but `AEConstructor` now coordinates insertion.
- Generic, sequencing, bulk-sequencing, and array paths register and reference `sample.data_processing` as a `Data-Processing` protocol. AE parsing maps processing/normalization protocols back to that sample field and maps only extraction protocols to `channel.extract_protocol`; library-construction and sequencing protocols are not collapsed into extraction. Sequencing paths also reference a populated `sample.scan_protocol`.
- Single-cell read/barcode/isolation comment labels emitted by the SDRF handlers are recognized by `AEParser` and do not create self-generated unmapped-column warnings.
- `_miniml2sdrf(data, protocol_registry=None, technology_type=None)` uses the supplied AE technology key or detects one for compatibility, selects a handler, builds the table, stores `last_sdrf_audit`, and returns rows.
- `_detect_sdrf_technology(data)` delegates to `AEConstructor._detect_ae_technology(data)`.
- `_has_array_files(data)` delegates to `AEConstructor._has_array_files(data)`.
- `_lookup_sra(sra)` delegates SRA fetching/parsing to `self.insdc_fetcher.fetch_sra_runs()` and returns `[]` on request or XML parse errors.

<a id="sdrf-file-helpers"></a>
#### File helpers

- `normalized_extension(path)` parses URLs/paths, strips one compression suffix, and returns the lowercase extension.
- `classify_file(path)` returns `sequencing_raw`, `array_raw`, `matrix_or_derived`, or `supplementary`.

<a id="base-sdrf-handler"></a>
#### Base SDRF handler

`class _BaseSDRFHandler`

- Initializes samples, platform lookup, sample lookup, series accession, protocol registry, audit object, and factor tags.
- Uses enriched `sample.sra_run` when present; otherwise uses the parent constructor's injected `insdc_fetcher` for SRA accession extraction and run lookup.
- `build()` orchestrates path building, column planning, and rendering.
- `build_paths()` creates generic source/factor paths.
- `plan_columns()`, `merge_column_group()`, `render_paths()`, `column_labels()`, `column_values()`, `path_groups()`, `group_with_values()`, `attr_columns()`, `occurrence_key()`, and `render_value()` handle table planning and rendering; `render_value()` preserves quote characters while normalizing missing values to blank cells.
- `ordered_samples()` respects series sample refs before remaining samples.
- `channels()` normalizes missing channels to `[{}]` and warns for multi-channel samples.
- `source_node()`, `sample_comment_attrs()`, `characteristic_attrs()`, `organism_part_value()`, `provider()`, and `material_type()` build mapped source columns.
- `characteristic_attrs()` seeds required blank source characteristics for organism, organism part, developmental stage, disease, and genotype; the first matching JSON value fills the seeded column and repeated values remain as repeated columns.
- `factor_nodes()`, `_factor_tags()`, `factor_value()`, `characteristic_values()`, and `characteristic_value()` handle experimental factors.
- `extraction_edges()` and `protocol_edge()` register protocol refs and warn for required blank refs.
- All SDRF handlers add a required sample collection `Protocol REF`; sequencing handlers also add a required nucleic acid sequencing `Protocol REF`.
- `sample_accession()`, `biosample_accessions()`, `biosample_accession()`, `platform()`, `platform_accession()`, `instrument_model()`, `supplementary_files()`, `raw_files()`, `derived_files()`, `arrayexpress_ftp()`, `file_node()`, `sra_runs()`, and `clean()` provide common extraction utilities.

<a id="sequencing-handlers"></a>
#### Sequencing handlers

- `_SequencingSDRFHandler.build_paths()` builds source, sample collection protocol, extraction, extract, library protocol, assay, nucleic acid sequencing protocol, scan, and factor nodes for each sample/channel/run.
- `extract_node()` maps material type and library attributes.
- `library_attrs()` maps library layout, selection, source, and strategy; source is uppercased for `Comment[LIBRARY_SOURCE]`.
- `geo_first_value()` prefers GEO over conflicting SRA values and records warnings.
- `library_protocol_text()` combines extraction/library/SRA fields for protocol descriptions.
- `assay_node()` maps technology type, ENA/SRA identifiers, submitted file, MD5, and instrument model.
- `geo_first_instrument_model()` prefers GEO instrument model over SRA.
- `scan_node()` maps scan name and sequencing file attrs.
- `sequencing_file_attrs()` maps FASTQs and raw sequencing files; the former derived data comment block is intentionally left commented out.
- `_BulkSequencingSDRFHandler` is selected for ordinary non-single-cell sequencing and expands each sample/channel/run into one path per FASTQ URI.
- `_SingleCellSequencingSDRFHandler` resolves each sample/channel/run through `resolve_chemistry()`, adds supported annotations, and records operation-local diagnostics. Unlabelled read-length lists do not supply cDNA lengths.
- `_DropletSingleCellSequencingSDRFHandler` and the retained v2/v3 rendering selections share evidence-based annotation. No version, barcode, primer, strand or read-length presets are supplied. See [scoped chemistry](#scoped-library-chemistry).
- `_PlateSingleCellSequencingSDRFHandler` inherits the bulk per-FASTQ row behavior and adds source-level index and description comments.
- `_SpatialSequencingSDRFHandler` adds Visium library construction, read geometry, and read type/read index comments based on submitted filenames.

<a id="array-and-generic-handlers"></a>
#### Array and generic handlers

- `_ArraySDRFHandler.build_paths()` builds source, extraction, labeled extract, hybridization, assay, scan, file, and factor nodes.
- `array_extract_node()` builds array extract nodes.
- `labeled_extract_node()` maps channel labels.
- `array_assay_node()` maps array assay technology and Array Design REF.
- `array_file_nodes()` maps raw image/data files, derived matrix files, and supplementary files.
- `_GenericSDRFHandler` uses the base source/factor path behavior.

<a id="legacy-fallback-notes"></a>
#### Legacy fallback notes

- `_GEOFallbackComments` and `_SRAFallbackComments` are commented out as reference code.
- No greedy fallback comments are emitted at runtime.

<a id="harmonizers"></a>
### `metadata/ontology_mappings.py`

`class GEO2OLS`

- Ensures an `ontologies` dict exists.
- Registers EFO and OBI term source metadata.

`geoprotocols2efo(protocol_type: str) -> list`

- Maps known MAGE-TAB/GEO protocol labels to ontology term, source ref, and accession.
- Maps `Library-Construction-Protocol` to `nucleic acid library construction protocol`, EFO, and `EFO_0004184`; registered EFO term-source metadata uses release `3.90.0`.
- Raises `ValueError` for blank protocol type.
- Returns `[protocol_type, None, None]` for unknown non-blank protocol labels, allowing custom protocol labels to survive in IDF output.

### `metadata/ontology_mappings.py`

`class Pubmed2OLS`

- Ensures an `ontologies` dict exists.
- Registers EFO and MeSH term source metadata.

`pubstatus2efo(pub_status: str) -> list`

- Returns `[None, None, None]` for blank status.
- Splits composite statuses on `+` and maps the first token.
- Maps common PubMed statuses such as `ppublish`, `epublish`, `pubmed`, `medline`, and `retracted`.
- Returns `[original_status_label, None, None]` for unknown non-blank statuses.

### `metadata/ontology_mappings.py`

`class Harmonizer(Pubmed2OLS, GEO2OLS)`

- Combines PubMed and GEO ontology mappings through multiple inheritance.
- Initializes a shared `ontologies` dictionary before calling parent initializers.

<a id="json-helper"></a>
### `helpers/json_helper.py`

`class JSONHandler`

`_from_path(obj, path_str)`

- Resolves dotted paths through dict/list structures.
- Numeric path components are treated as list indexes.
- `*` expands over lists.
- Missing branches return `[None]`, preserving positional behavior for callers.

`_flatten_values(value)`

- Recursively flattens nested lists and dict values.
- Returns scalar values as a one-item list.

<a id="request-helper"></a>
### `helpers/request_helper.py`

`class NCBIApplicationIdentity`

- Validates a 1-64 character NCBI tool identifier, optional contact email, and
  optional nonblank API-key token. The secret is excluded from `repr`.
- `params() -> dict[str, str]` returns configured `tool`/`email`/`api_key`
  parameters used
  by both NCBI fetchers. Defaults identify the released MSC application and its
  public maintainer contact; callers may inject an approved replacement. A
  missing email emits one warning and does not relax conservative pacing.
- Request telemetry contains service/host/attempt/status only and never logs
  these parameters.

`class HostRequestCooldownDeferred`

- Signals that a persisted provider cooldown exceeds the caller's inline wait
  budget. `retry_at` is an absolute Unix timestamp suitable for a retryable
  checkpoint; it is not a request URL or provider payload.

`class HostRequestGate`

- Normalizes a provider/model key and stores only versioned request-start and
  cooldown timestamps in an owner-only state file selected by SHA-256.
- Uses `flock` across processes and keeps the lock through the bounded wait and
  timestamp update, so unrelated traces under the same Unix user share one
  conservative start schedule.
- `slot(...)` adds a distinct one-at-a-time cross-process lease around a caller's
  provider operation and applies the same start pacing after acquisition. It is
  intended for model keys that must not overlap; cancellation while waiting
  raises `InterruptedError`, and the lease always releases on context exit.
- Honors numeric and HTTP-date `Retry-After` values. A cooldown beyond
  `max_wait_seconds` raises `HostRequestCooldownDeferred` instead of sleeping
  past a worker's budget.
- Uses `SCIENTIFIC_PROVIDER_GATE_DIR` when explicitly configured, otherwise an
  XDG runtime directory when writable, with an owner-only per-user `/tmp`
  fallback. Explicit invalid, unowned, or symbolic-link paths fail closed.

`class RequestSettings`

- Stores request behavior: `timeout`, `request_delay`, `max_in_flight`,
  `max_retries`, retry HTTP statuses, exponential backoff base, and maximum
  backoff plus `max_inline_wait`. Invalid time, delay, concurrency, retry, or
  wait values fail at construction.
- Defaults retry HTTP statuses to `{403, 408, 425, 429, 500, 502, 503, 504}`;
  ordinary non-throttling 4xx responses are not retried.

`DEFAULT_REQUEST_SETTINGS`

- `ncbi_eutils`: timeout 30 seconds, request delay 0.5 seconds, two maximum
  in-flight requests, and 3 retries.
- `geo_ftp`, `biostudies`, and `ena_portal`: timeout 30 seconds, request delay
  1.0 seconds, two maximum in-flight requests, and 3 retries.

`class RateLimitedRequester`

- `get(url: str, **kwargs: Any) -> requests.Response` wraps `requests.get()`, applies a default timeout,
  enforces host-wide delay and in-flight limits, retries configured statuses, and returns a response
  or raises the exhausted HTTP/transport error.
- `reset_service_state()` is a class-level test/operations hook that clears
  shared limiter timestamps; it mutates process-global requester state.
- Maintains process-local in-flight limits and delegates every actual request
  attempt to `HostRequestGate`, so separate processes and libraries targeting
  the same host share start pacing and cooldowns. Different hosts remain
  independent.
- Retries transient HTTP statuses. Numeric or HTTP-date `Retry-After` controls
  the persisted cooldown; otherwise full jitter selects a value from zero to
  `min(0.5 * (2 ** attempt), 8.0)`.
- Retries `ConnectionError`, `Timeout`, and `ChunkedEncodingError` with the same
  bounded full-jitter policy and exact configured attempt count.
- Raises the exhausted retry response through `response.raise_for_status()`.
- Exposes cumulative `provider_attempts`, `retry_count`, and
  `rate_wait_seconds` counters for the lifetime of the requester. Gate waits are
  counted once at the actual provider-attempt boundary; cached work never
  changes them.

<a id="pubmed-fetcher"></a>
### `sources/pubmed.py`

`class PubmedWebFetcher`

`__init__(requester=None, request_settings=None, ncbi_identity=None,
resource_profile="standard", resource_overrides=None)`

- Defaults to `RateLimitedRequester(service="ncbi_eutils")`.
- Accepts a custom requester, NCBI E-utilities request settings, or validated
  application identity.

`fetch_pubmed_summary(pubmed_id: str) -> ET.Element`

- Calls NCBI PubMed ESummary for one PubMed ID through the `ncbi_eutils`
  requester with the configured tool/contact parameters.
- Raises for HTTP errors and returns the parsed XML root.

`pubmed_summary(pubmed_id: str) -> tuple`

- Parses DOI, author list, title, and PubMed publication status from ESummary XML.
- Maps publication status through `Harmonizer().pubstatus2efo()`.
- Returns the existing IDF tuple shape: DOI, author string, title, mapped status, term source ref, and term accession.

<a id="insdc-fetcher"></a>
### `sources/insdc.py`

`class INSDCWebfetcher`

`__init__(ncbi_requester=None, ena_requester=None,
ncbi_request_settings=None, ena_request_settings=None, ncbi_identity=None,
resource_profile="standard", resource_overrides=None)`

- Defaults to `RateLimitedRequester(service="ncbi_eutils")` for NCBI SRA EFetch.
- Defaults to `RateLimitedRequester(service="ena_portal")` for ENA Portal file reports.
- Accepts custom requesters, per-service request settings, or a validated NCBI
  application identity. EFetch receives its tool/contact parameters; ENA calls
  do not receive NCBI-specific fields.

`_extract_sra(sra: str)`

- Extracts SRA/ENA/DDBJ-style accessions matching `[SED]R[RXSP]` plus digits.
- Matching is case-insensitive.

`_ncbi_nrx(nrx: str)`

- Calls NCBI SRA EFetch with `retmode=xml` through the `ncbi_eutils` requester.
- Raises for HTTP errors and returns the parsed XML root.

`fetch_sra_runs(accession: str) -> list`

- Calls `_ncbi_nrx()` and ENA Portal `filereport`, then parses and merges run metadata into the run dictionaries consumed by SDRF handlers.
- Preserves the existing run record shape and adds `study`: library layout/source/strategy/selection, SRA/ENA study/sample/run IDs, GEO sample ID, BioSample ID, instrument model, submitted FASTQ filename, MD5, read lengths, and per-FASTQ `filename`/`uri`/`md5` records.
- ENA `fastq_ftp` links are preferred for FASTQ `uri`; if ENA links are absent or unavailable, original NCBI SRA XML `url`/`Alternatives` links are used as fallback.

ENA file report helpers:

- `fetch_ena_file_report()` calls `https://www.ebi.ac.uk/ena/portal/api/filereport` through the `ena_portal` requester with `result=read_run`, FASTQ fields, and JSON output.
- `fetch_ena_fastq_files()` groups parsed ENA FASTQ records by `run_accession`.
- `_parse_ena_fastq_report()`, `_split_ena_file_field()`, `_normalize_ena_ftp_uri()`, and `_filename_from_uri()` parse semicolon-delimited ENA file fields.

SRA XML helper methods:

- `_parse_sra_library()`, `_parse_sra_sample_ids()`, `_parse_sra_instrument_model()`, `_parse_sra_fastqs()`, `_element_accession()`, `_find_text()`, `_strip_ns()`, and `_clean_sdrf_text()` support SRA parsing.

- Because it returns `None`, normal validation currently raises `AssertionError`.

<a id="maintenance-notes"></a>
## Maintenance Notes

- Documentation should describe implemented behavior, not planned burndown items.
- The SDRF system intentionally maps known GEO/SRA values and leaves the old greedy fallback comment code disabled.
- Protocol rows in the IDF should be driven by the same `ProtocolRegistry` used while building SDRF paths.
- PubMed and SRA lookups normally happen during parsed MINiML enrichment through injected fetchers, so unit tests should mock or fake `MINiMLEnricher`, `PubmedWebFetcher`, and `INSDCWebfetcher` when exercising conversion behavior.
- Parser-specific behavior is documented in this handoff under GEO Parse Flow and GEO parser.
- Generated outputs under `.dev/` and `output/` are examples/debug artifacts, not package code.
- H5AD sources remain immutable; normalization writes new files and records source hashes and provenance.
- nf-core revisions are pinned in `NFCoreRunner.REVISIONS`; revision changes require command/output-discovery tests and documentation updates.
- Rootless Compose deliberately mounts only `.out/json2h5ad`. Inputs required by nested nf-core containers must be copied or generated below that path.
- Socket access still grants full control of the dedicated rootless daemon; keep `nfcore-runner` locked and deny it unrelated files and credentials.

<a id="harmonization-overrides"></a>
## Harmonization overrides

`JSONPackageSource` recognizes an Agentic Curator result envelope containing
`miniml_json` and an optional sibling `harmonization_overrides` profile. The
profile remains attached to each dataset group and is inactive unless a JSON
consumer receives `use_harmonization_overrides=True` or the corresponding CLI
flag. Native MINiML and Atlas v1 behavior is unchanged.

The shared resolver validates profile schema version `1.0`, fixed MSC
destinations or `characteristics.<normalized_tag>`, and ordered canonical
source fields. The first populated typed annotation whose `field` matches a
configured source wins; its value, term source, term accession, and hierarchy
depth drive the resolved view. Invalid profiles warn and fall back to the
complete raw package. Resolution uses a deep copy, so destination replacement
never mutates the canonical typed input.

H5AD/obs publish canonical schema-3 columns, harmonization provenance columns,
and `uns["msc_harmonization"]`; `uns["msc_miniml"]` retains the untouched source.
TSV/CSV publish the same canonical and provenance view. MAGE-TAB replaces its
semantic destination and emits standard ontology companions where available;
private `hz_*` columns are never rendered. ECTO exposure and PCL cell-state
annotations pass through json2ae/ae2json, json2tsv, json2h5ad, and json2obs as
typed annotations without format-specific ontology code.

Importable implementation symbols are
`meta_standards_converter.metadata.harmonization_overrides.HarmonizationSelection`,
`meta_standards_converter.metadata.harmonization_overrides.HarmonizationResolution`,
`meta_standards_converter.metadata.harmonization_overrides.resolve_harmonization_overrides`,
and
`meta_standards_converter.metadata.harmonization_overrides.validate_harmonization_overrides`.

<a id="test-plan"></a>
## Test Plan

Run the full suite:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/matplotlib-meta-standards \
  pytest -p no:cacheprovider -q
```

Install `.[test]` first. Pytest is canonical because unittest discovery skips
pytest functions, fixtures, parametrization, and subtests. The autouse guard in
`tests/conftest.py` blocks network and child processes fail-closed. Only
reviewed bounded-PATH shell tests carry the `fake_process` marker.

Important test coverage:

- `tests/miniml/test_geo_parser.py`: parser package scoping, cardinality, namespace handling, empty cleanup, related-series traversal, and fixture-backed parsing with `tests/GSE328265_family.xml`.
- `tests/converters/test_geo2ae.py`: converter orchestration, related-series forwarding, enrichment, stage logging, and `remove_empty` forwarding.
- `tests/converters/test_geo2json.py`: JSON converter orchestration, optional enrichment, JSON file writing, and stage logging.
- `tests/converters/test_json2ae.py`: object/list loading, validation, default and skipped enrichment, MAGE-TAB writing, safe logging, independent fixture expectations, and extension restoration.
- `tests/converters/test_ae2json.py`: IDF/SDRF mapping, typed protocol/declaration/assay-path capture, model edit authority, assay multiplicity, units/ontology, sidecar/fingerprint creation, unchanged lossless reuse, edited-core precedence, keyed IDF/SDRF overlay union, occurrence-aware duplicate headers, harmonized `hz_*` columns, ambiguity-safe row alignment, multiple SDRFs, conflicts, unmapped restoration, frozen strict E-MTAB-6486 normalization, and output writing.
- `tests/sources/test_ae_webfetcher.py`: bounded local and streamed HTTPS resolution, typed profile propagation, explicit host policy, explicit SDRF overrides, paginated BioStudies discovery/download calls, in-memory remote content, and invalid source metadata.
- `tests/expression/test_json2h5ad.py`: asset precedence/manifests/downloads, canonical H5AD schema 1 metadata, normalized multivalue rows, smart observation IDs, opaque source-column preservation, real dictionary reference scoping, artifact-relative provenance, MINiML enrichment and publication filtering, count/TPM matrices, catalogue-only output, fail-closed compatibility evidence, Entrez/symbol separation, canonical/generic study splitting, correlation-safe partial-failure summaries, partial results, and raw-output reintegration.
- `tests/test_retrieval.py`: host/address/redirect policy, cache integrity,
  byte/disk/aggregate ceilings, and bounded NCBI range fallback behavior.
- `tests/test_atlas_v1_reader.py`: producer-owned golden fixture consumption, harmonized-state adaptation, structural validation, v1 cutover failure, and no-ThematicAtlases dependency proof.
- `tests/test_json_source.py`: native MINiML and Atlas v1 grouping, harmonized-status filtering, source diagnostics, and duplicate conflict handling.
- `tests/converters/test_json2tsv.py`: neutral default columns, direct Atlas aggregation, injected neutral metadata services, replacement projectors, collisions, and validation behavior.
- `tests/test_miniml_model_authority.py`: sample-bound typed protocol/material projection, exact ontology preservation, fallback ordering, and shared H5AD semantics.
- `tests/test_miniml_stabilization.py`: deterministic MINiML migration, validation, and captured-index ordering without quadratic equality scans.
- `tests/magetab/test_magetab_miniml_v2.py`: legacy and canonical IDF companion-label parsing, typed ontology alignment, and canonical semantic MAGE-TAB regeneration.
- `tests/metadata/test_metadata_projector.py`: generic sample projector and ignored legacy combined-hook
  lifecycle, scalar broadcasting, axis-length validation, collision rejection,
  warning/error propagation, fail-closed output, invalid-output opt-in, and
  bundle rollback fault injection.
- `tests/test_external_guard.py`: fail-closed network/process guard self-tests and bounded fake-process opt-in.
- `tests/e2e/`: independent public converter and CLI workflows against reviewed fixtures; `tests/expression/test_lazy_scientific_imports.py` checks processed H5AD conversion without importing Scanpy.
- `tests/expression/test_h5ad_pipeline.py`: reference/annotation combinations, GFF3 conversion and reuse, FASTQ samplesheets, mixed modality grouping, pinned commands, warning extraction, output discovery, and workflow failure logs.
- `tests/expression/test_h5ad_pipeline.py`: rootless enforcement also covers accepted, rootful, and unreachable Docker daemons.
- `tests/test_docker_artifacts.py`: pinned runtime tooling, rootless-only Compose mounts, hardening, and provisioning/runner script syntax.
- `tests/cli/test_cli_geo2ae.py`: CLI defaults, multiple accession order, aliases, keep-empty behavior, out directory forwarding, logging controls, file logging, and failure continuation.
- `tests/cli/test_cli_geo2json.py`: JSON CLI defaults, enrichment toggle, aliases, logging controls, file logging, and failure continuation.
- `tests/cli/test_cli_json2ae.py`: JSON-to-MAGE-TAB CLI defaults, enrichment toggle, multiple input ordering, output forwarding, logging, and failure continuation.
- `tests/cli/test_cli_ae2json.py`: MAGE-TAB-to-JSON CLI defaults, repeated SDRF overrides, source validation, multiple input ordering, logging, and failure continuation.
- `tests/cli/test_cli_json2h5ad.py`: H5AD CLI defaults, workflow/reference/asset flags, partial status, multiple input order, logging, and failure continuation.
- `tests/cli/test_cli_json2tsv.py`: TSV/CSV input order and partial exit status.
- `tests/test_project_scripts.py`: console script registration.
- `tests/policy/test_documentation_policy.py`: stable documentation anchors, required README Guide structure including configuration, complete Mermaid platform-handler hierarchy coverage, interface-specific quickstart links, live-parser coverage for every documented CLI argument and alias, console-script mentions, docs links, and author-header policy.
- `tests/magetab/test_ae_constructor.py`: IDF rows, merged and source-aligned secondary accessions, protocol registry behavior, AE constructor sequencing, SDRF row insertion, file normalization, and protocol ref consistency.
- `tests/magetab/test_ae_sdrf_handlers.py`: SDRF graph rendering, source/comment/characteristic behavior, file classification, sequencing/array/single-cell/spatial handlers, SRA precedence warnings, and disabled greedy fallback comments.
- `tests/metadata/test_miniml_enricher.py`: additive PubMed/SRA enrichment fields, deduplication, and fetch error tolerance.
- `tests/test_request_helper.py`: timeout forwarding, cross-process host pacing
  and model-slot serialization, persisted cooldowns, retry statuses,
  `Retry-After`, exponential full jitter, exhausted retry errors, and request,
  retry, and rate-wait counters.
- `tests/sources/test_geo_webfetcher.py`: GEO URL handling, requester delegation, and MINiML tarball extraction.
- `tests/sources/test_insdc_webfetcher.py`: SRA accession extraction, NCBI/ENA requester delegation, parsed SRA run records, and ENA fallback behavior.
- `tests/sources/test_pubmed_webfetcher.py`: PubMed ESummary requester delegation, parsing, publication status mapping, and IDF constructor delegation.

<a id="live-api-provider-contracts"></a>
### Live API provider contracts

`tests/live_api/test_public_provider_contracts.py` owns PubMed ESummary, NCBI
SRA EFetch plus ENA file-report, and paginated BioStudies IDF/SDRF contracts.
It uses 10-second timeouts, no retries/delays, a 12-send ceiling, and declared
NCBI/EBI hosts. Run it with
`RUN_LIVE_API_TESTS=1 python -m pytest tests/live_api -m live_api -vv`.

After documentation edits, also check:

```bash
rg -n "^(#|##|###) " docs/codebase.md docs/index.md
```

Then retrieve representative anchors with the commands in `docs/index.md` to confirm each index entry resolves to its intended section.

<a id="durable-artifact-publication"></a>
## Durable artifact publication

`JSON2TSVConverter` and `JSON2OBSConverter` publish related tabular and AnnData metadata
artifacts with `DurableArtifactBundlePublisher`. Each commit copies staged files
into a new immutable generation, records size and SHA-256 values in a schema-1.0
generation manifest, fsyncs every file and directory, refreshes direct-file
compatibility views, and atomically replaces one `current.json` pointer.
`resolve_current_bundle` validates the pointer, manifest, contained paths, and
content digests. Pointer-aware consumers therefore see one complete generation;
legacy filenames remain a transition view without bundle-wide atomic visibility.

If publication fails before pointer commit, prior compatibility views are
restored. A restoration failure raises `ArtifactRecoveryError` with the original
publication error, all recovery errors, and surviving recovery paths, which are
deliberately retained. The public surface also exports
`PublishedArtifactBundle` and `resolve_current_bundle`; result envelopes expose
the additive `bundle_pointer_path` field without changing H5AD metadata schema
1.0 or Atlas v1 inputs. Qualified public symbols are
`meta_standards_converter.converters.DatasetBundleRecoveryError` and
`meta_standards_converter.expression.catalogue.DatasetBundleRecoveryError`,
`meta_standards_converter.artifact_bundle.ArtifactRecoveryError`,
`meta_standards_converter.artifact_bundle.DurableArtifactBundlePublisher`,
`meta_standards_converter.artifact_bundle.PublishedArtifactBundle`, and
`meta_standards_converter.artifact_bundle.resolve_current_bundle`.
[source](../src/meta_standards_converter/artifact_bundle.py)

The publisher also backs up the previous pointer before replacement. A failure
while making the new pointer durable restores both that pointer and the legacy
views, preventing consumers from observing different generations through the
two interfaces.

<a id="neutral-ae-construction-state"></a>
## Neutral AE construction state

`magetab/technology.py` owns `ProtocolRegistry`, normalized file-extension
classification, array-file detection, and platform technology selection. Both
`AEConstructor` and `SDRFConstructor` import these contracts in one direction.
This removes their prior mutual import and method-local constructor imports
without changing handler keys, protocol references, technology decisions, or
the public `ae_constructor.ProtocolRegistry` import path.
Importable symbols are
`meta_standards_converter.magetab.technology.ProtocolRegistry`,
`meta_standards_converter.magetab.technology.detect_ae_technology`,
`meta_standards_converter.magetab.technology.has_array_files`, and
`meta_standards_converter.magetab.technology.normalized_extension`, and
`meta_standards_converter.magetab.technology.series_identity`. Study
identity follows the Python model: a usable `series.iid` takes precedence over
the first usable accession, and generated MAGE-TAB protocol identifiers use the
same value. GEO-shaped `GSE` identities remain numeric-only.

<a id="operational-events-v1"></a>
## Operational events v1

`meta_standards_converter.operational_events.OperationalEventEmitter`
implements the cross-repository schema-1.0 envelope for structured logging and
optional JSONL. It recursively redacts credential-like keys and URL query
strings, aggregates event/status counts and durations, and atomically exports
aggregate JSON or Prometheus text. `RateLimitedRequester` accepts an optional
emitter and reports terminal provider success/failure with safe host, attempt,
status-code, and duration fields; request parameters and response bodies never
enter events. Qualified symbols are
`meta_standards_converter.operational_events.EVENT_SCHEMA_VERSION`,
`meta_standards_converter.operational_events.OperationalEventEmitter`, and
`meta_standards_converter.operational_events.redact`.

<a id="runtime-contracts-v2"></a>
## Runtime contracts v2

`meta_standards_converter.runtime_contracts.OperationStatusV2` is the clean
status-2.0 boundary. It keeps execution, completeness, evidence confidence,
validation, and publication disposition independent, preserves stable terminal
reason codes, and rejects legacy/unversioned mappings. A successful lookup with
no scientific match is therefore representable as operationally succeeded and
complete while evidence is insufficient and publication requires review.
Failed/invalid/incomplete work cannot claim a publishable disposition, and
`safe_to_publish` is true only for succeeded, complete, valid, explicitly
publishable results. `aggregate()` applies deterministic worst-axis rules while
retaining item-safe errors.

`meta_standards_converter.runtime_contracts.SafeErrorEnvelope` is the mandatory
persistence-safe error shape. `from_exception()` records only the exception
class, provider, URL-without-userinfo/query/fragment or path basename, HTTP
status, retry category, stage/item, and correlation ID. Raw exception text is
not serialized. `meta_standards_converter.cli.common.record_safe_cli_error`
is the single CLI boundary: every converter logs only this safe metadata, never
a traceback or raw exception message, and structured CLI summaries persist the
same envelope. The independent enums and error/status types are:

- `meta_standards_converter.runtime_contracts.ExecutionStatus`;
- `meta_standards_converter.runtime_contracts.CompletenessStatus`;
- `meta_standards_converter.runtime_contracts.EvidenceConfidence`;
- `meta_standards_converter.runtime_contracts.ValidationStatus`;
- `meta_standards_converter.runtime_contracts.PublicationDisposition`;
- `meta_standards_converter.runtime_contracts.RetryCategory`.

`meta_standards_converter.runtime_contracts.ResourceProfile` owns typed
standard/large ceilings. Standard uses 3 redirects, 10/60-second connect/read
timeouts, 128 MiB XML, 256 MiB compressed and 1 GiB expanded archives, 3 GiB
per ontology file, 100 GiB per matrix/H5AD, an 8 GiB in-memory matrix ceiling,
200 GiB aggregate downloads,
250 GiB cache, four network workers, and one ontology-build worker. Large uses
5 redirects, 30/300 seconds, 512 MiB XML, 1/4 GiB archive limits, 5 GiB per
ontology file, 500 GiB matrix/H5AD, a 32 GiB in-memory matrix ceiling, 1 TiB
aggregate downloads, 500 GiB cache,
eight network workers, and two ontology-build workers. Ontology and matrix byte
ceilings remain disk/download object sizes; the separate in-memory ceilings
are RAM admission limits. Normal conversion also applies 70% of currently
available host/cgroup RAM, whichever is lower, and resume-only force admission
is hard-bounded at 90%. Both profiles require 10% disk headroom.
`--resource-profile`, repeated `--resource-override`, and the
Python constructors expose explicit selection/overrides. `geo2ae`, `geo2json`,
`ae2json`, and `json2h5ad` expose the shared CLI options; default network
collaborators receive the same immutable configured profile. Passing an
already-configured `ResourceProfile` preserves all existing overrides unless
new explicit overrides replace named fields.

The shared CLI parsing/configuration symbols are
`meta_standards_converter.cli.common.add_resource_profile_arguments`,
`meta_standards_converter.cli.common.parse_resource_override`, and
`meta_standards_converter.cli.common.configured_resource_profile`.

Qualified resource/disk symbols are
`meta_standards_converter.runtime_contracts.DiskBudget`,
`meta_standards_converter.runtime_contracts.DiskBudgetError`,
`meta_standards_converter.runtime_contracts.get_resource_profile`, and
`meta_standards_converter.runtime_contracts.require_disk_headroom`.

<a id="secure-retrieval-and-xml"></a>
## Secure retrieval and XML boundaries

`meta_standards_converter.retrieval.RetrievalPolicy` accepts only configured
schemes and exact/provider-suffix hosts, rejects URL userinfo and non-public
IPv4/IPv6 answers, and revalidates each same-scheme redirect.
`meta_standards_converter.retrieval.RetrievalService` streams assets with
connect/read timeouts and object, aggregate-run, cache, and disk-headroom
checks. An exclusive cache lock covers verification, capacity reservation,
streaming, and atomic publication. Before consuming a body, one cache snapshot
and one disk preflight reserve the declared response size, or the full object
ceiling when length is unknown; chunk processing updates only byte counters and
digests and never rescans the directory. Successful reuse atomically refreshes
the sidecar's `last_used_at`. Retention reporting is locked and dry-run by
default, verifies candidate SHA-256/size/sidecar contracts, preserves explicit
active references and a configurable newest minimum, and can only quarantine
with recovery paths—never delete. For the exact configured
`ftp.ncbi.nlm.nih.gov` host, an HTTP 403 may
activate bounded 16 MiB HTTPS range requests; every request repeats URL/DNS
validation, redirects remain disabled, `Content-Range` must be contiguous and
truthful, and the same object/run/cache/disk limits apply. Cache publication
records sanitized origin, byte count, SHA-256,
optional MD5, fetch time, and profile; reuse fails closed if bytes or sidecar
do not match. `meta_standards_converter.retrieval.AssetDownloader` preserves
the former supported facade while delegating to this service. Failure types are
`meta_standards_converter.retrieval.RetrievalError`,
`meta_standards_converter.retrieval.RetrievalSecurityError`,
`meta_standards_converter.retrieval.RetrievalSizeError`, and
`meta_standards_converter.retrieval.CacheIntegrityError`.

`meta_standards_converter.xml_safety.read_limited_response` and
`meta_standards_converter.xml_safety.stream_limited_response` enforce declared
and actual decoded-body limits. A compressed transport's `Content-Length`
describes its encoded body and is therefore not compared with decoded bytes;
the declared and decoded ceilings are still independently enforced.
`meta_standards_converter.xml_safety.parse_xml` accepts and removes one ordinary
external SYSTEM/PUBLIC DTD without resolving it, then uses standard-library
parsing. Entity declarations, internal subsets, malformed/multiple DTDs, and
DTDs outside the prolog are rejected. Errors are typed as
`meta_standards_converter.xml_safety.XMLSafetyError`,
`meta_standards_converter.xml_safety.XMLSizeLimitError`, and
`meta_standards_converter.xml_safety.UnsafeXMLDocumentError`. GEO retrieval
streams the compressed response to a disk-preflighted temporary archive,
requires exactly one root `{GSE}_family.xml`, permits safe auxiliary regular
files/directories, validates member count and aggregate expanded/XML size, and
rejects unsafe paths, duplicate names, links, special files, and unexpected
XML without calling `extractall`. SRA/PubMed XML uses the same bounded parser,
and PubMed uses HTTPS. `AEWebFetcher` applies the same `RetrievalPolicy` host, credential,
public-address, and redirect checks to BioStudies API/file URLs and explicit
IDF/SDRF URLs; it bounds UTF-8 API JSON and MAGE-TAB text per file and across
the run. Provider suffixes are trusted by default. Additional exact explicit
source hosts require `ae2json --source-host HOST` or an injected policy.

<a id="msc6-source-services"></a>
## MSC 6 source services

`GEOSource` owns retrieval and related-series traversal; `GEOParser` parses supplied
XML without network calls. Converter-specific enrichment defaults and guarded
parent-publication inheritance remain unchanged. `JSONPackageSource` owns JSON
recognition and dataset grouping. INSDC accepts explicit clients with
`fetch_sra_xml(nrx)` and `fetch_ena_file_report(accession)` methods; ThematicAtlases
implements checkpoint interception by composition. `metrics()` returns cumulative
request snapshots without exposing nested requesters. Existing MINiML, source
evidence and checkpoint serialization remain unchanged.

### Added source service callable inventory
- `meta_standards_converter.sources.contracts.RequestMetrics`: [contracts.py](../src/meta_standards_converter/sources/contracts.py#L14).
- `meta_standards_converter.sources.contracts.MetricsProvider`: [contracts.py](../src/meta_standards_converter/sources/contracts.py#L19).
- `meta_standards_converter.sources.contracts.INSDCClient`: [contracts.py](../src/meta_standards_converter/sources/contracts.py#L22).
- `meta_standards_converter.sources.contracts.PubMedClient`: [contracts.py](../src/meta_standards_converter/sources/contracts.py#L26).
- `meta_standards_converter.sources.contracts.request_metrics`: [contracts.py](../src/meta_standards_converter/sources/contracts.py#L29).
- `meta_standards_converter.sources.geo.RelatedSeriesParseResult`: [geo.py](../src/meta_standards_converter/sources/geo.py#L210).
- `meta_standards_converter.sources.geo.GEOSource`: [geo.py](../src/meta_standards_converter/sources/geo.py#L235).

<a id="msc6-service-architecture"></a>
## MSC 6 service architecture and migration

The converter API is a coordinated breaking release. CLI commands and serialized scientific contracts remain unchanged. Consumers import types from their owning packages; retired converter modules and lowercase classes are removed.

Import converter classes from `meta_standards_converter.converters` and projection types such as `AnnDataMetadataProjection` and `TabularMetadataProjection` from `meta_standards_converter.metadata.projection`. The README's custom-projector examples use these public exports.

```text
cli -> converters (one workflow per module)
          |-> sources -> request policy / retrieval -> providers
          |-> miniml.geo_parser -> typed MINiML packages
          |-> metadata.enrichment -> injected PubMed / INSDC clients
          |-> magetab.constructor -> SDRF handlers / IDF -> writer
          |-> sources.json -> dataset groups -> metadata.projection.tabular
          `-> expression.planning -> selected assets
                -> injected reader / nfcore -> normalization -> projection
                -> version-specific checkpoints -> catalogue publication
json2obs -> JSON2H5ADConverter -> AnnDataComponentExporter -> artifact bundle
```

`SourcePlanner` retains precedence and selection policy; inject `AssetDiscovery` to supply candidates. `ProcessedAssetReader.read(asset, orientation=..., localize=...)` owns selected-asset reading, with scientific imports delayed until use. `JSON2H5ADConverter` injects readers and projectors, delegates metadata normalization to `AnnDataNormalizer`, projection execution to `AnnDataProjectorRunner`, processed state to `ProcessedCheckpointStore`, and sample catalogue publication to `CataloguePublisher`. Observation export aggregates rows without integrating sample matrices.

`GEOParser` is network-free. `GEOSource` separately retrieves original XML and collects related studies. `INSDCClient.fetch_sra_xml` and `fetch_ena_file_report` permit checkpoint-aware interception. Public cumulative `RequestMetrics` snapshots replace nested requester traversal. ThematicAtlases saves original source evidence before parsing and retains incremental enrichment checkpoints; its checkpoint path does not gain ordinary converter parent-publication inheritance.

Protocol registries, SDRF handler audits and caches remain scoped to each build. Shared request gates and retrieval caches retain their existing ownership. Technology inheritance remains intact. MINiML models, ontology policy, Atlas lifecycle, runtime contracts, XML safety, and artifact publication retain their existing responsibilities.

Processed checkpoint fingerprints retain the existing payload and SHA-256 algorithm. Destinations now include a hash of the MSC package version. A different version cannot reuse or overwrite historical processed checkpoints. Stored Atlas envelopes are unchanged; production resume still requires a supported, validated transition or a new run. No live cutover is performed by this refactor.

### Owning service symbol reference

- `meta_standards_converter.expression.components.AnnDataComponentExporter`: `AnnDataComponentExporter()`; [source](../src/meta_standards_converter/expression/components.py).
- `meta_standards_converter.expression.normalization.AnnDataNormalizer`: `AnnDataNormalizer(self, *, metadata_service, planner, localize, package_version, combination_policy=None)`; [source](../src/meta_standards_converter/expression/normalization.py).
- `meta_standards_converter.metadata.projection.anndata.AnnDataProjectorRunner`: `AnnDataProjectorRunner(self, projectors=())`; [source](../src/meta_standards_converter/metadata/projection/anndata.py).
- `meta_standards_converter.expression.planning.AssetDiscovery`: `AssetDiscovery()`; [source](../src/meta_standards_converter/expression/planning.py).
- `meta_standards_converter.expression.readers.AssetReader`: `AssetReader()`; [source](../src/meta_standards_converter/expression/readers.py).
- `meta_standards_converter.expression.catalogue.CataloguePublisher`: `CataloguePublisher()`; [source](../src/meta_standards_converter/expression/catalogue.py).
- `meta_standards_converter.magetab.sdrf.model.ColumnGroup`: `ColumnGroup()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.expression.planning.DefaultAssetDiscovery`: `DefaultAssetDiscovery()`; [source](../src/meta_standards_converter/expression/planning.py).
- `meta_standards_converter.sources.contracts.GEOXMLParser`: `GEOXMLParser()`; [source](../src/meta_standards_converter/sources/contracts.py).
- `meta_standards_converter.sources.contracts.MAGETabSourceResolver`: `MAGETabSourceResolver()`; [source](../src/meta_standards_converter/sources/contracts.py).
- `meta_standards_converter.magetab.writer.MAGETabWriter`: `MAGETabWriter()`; [source](../src/meta_standards_converter/magetab/writer.py).
- `meta_standards_converter.metadata.enrichment.MetadataEnrichment`: `MetadataEnrichment()`; [source](../src/meta_standards_converter/metadata/enrichment.py).
- `meta_standards_converter.sources.contracts.PackageLoader`: `PackageLoader()`; [source](../src/meta_standards_converter/sources/contracts.py).
- `meta_standards_converter.expression.readers.ProcessedAssetReader`: `ProcessedAssetReader()`; [source](../src/meta_standards_converter/expression/readers.py).
- `meta_standards_converter.expression.checkpoints.ProcessedCheckpointStore`: `ProcessedCheckpointStore(self, package_version)`; [source](../src/meta_standards_converter/expression/checkpoints.py).
- `meta_standards_converter.magetab.protocols.ProtocolRegistry`: `ProtocolRegistry(self, series_accession: str)`; [source](../src/meta_standards_converter/magetab/protocols.py).
- `meta_standards_converter.magetab.sdrf.model.SDRFAttr`: `SDRFAttr()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.magetab.sdrf.model.SDRFAudit`: `SDRFAudit()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.magetab.sdrf.model.SDRFEdge`: `SDRFEdge()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.magetab.sdrf.model.SDRFNode`: `SDRFNode()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.magetab.sdrf.model.SDRFPath`: `SDRFPath()`; [source](../src/meta_standards_converter/magetab/sdrf/model.py).
- `meta_standards_converter.magetab.sdrf.renderer.SDRFRenderer`: `SDRFRenderer()`; [source](../src/meta_standards_converter/magetab/sdrf/renderer.py).
- `meta_standards_converter.expression.assets.classify_asset`: `classify_asset(path: str | None)`; [source](../src/meta_standards_converter/expression/assets.py).
- `meta_standards_converter.magetab.sdrf.handlers.base.classify_file`: `classify_file(path: str)`; [source](../src/meta_standards_converter/magetab/sdrf/handlers/base.py).
- `meta_standards_converter.expression.readers.read_h5ad`: `read_h5ad(anndata, path: str)`; [source](../src/meta_standards_converter/expression/readers.py).
- `meta_standards_converter.expression.readers.scanpy_module`: `scanpy_module()`; [source](../src/meta_standards_converter/expression/readers.py).
- `meta_standards_converter.expression.components.scientific_modules`: `scientific_modules()`; [source](../src/meta_standards_converter/expression/components.py).
- `meta_standards_converter.expression.readers.scientific_modules`: `scientific_modules()`; [source](../src/meta_standards_converter/expression/readers.py).
- `meta_standards_converter.expression.readers.underlying_suffix`: `underlying_suffix(path: str)`; [source](../src/meta_standards_converter/expression/readers.py).

### Explicit MAGE-TAB evidence orchestration

`meta_standards_converter.metadata.enrichment.MAGETabEvidenceResolver` accepts injected PubMed and INSDC clients. `AEConstructor` creates an operation-local SDRF handler, resolves its ordered sequencing run evidence, renders protocols and paths, validates IDF prefix rows, resolves missing publication details, and constructs the remaining IDF rows. Enriched run/publication evidence suppresses retrieval; array and generic handlers do not trigger SRA retrieval. The SRA fallback still catches only request and XML parsing errors. IDF and SDRF builders no longer create network clients. [Evidence resolver](../src/meta_standards_converter/metadata/enrichment.py), [orchestration](../src/meta_standards_converter/magetab/constructor.py).

Converter-focused tests mirror `sources`, `miniml`, `magetab`, `metadata`, `expression`, `converters`, and `cli` boundaries under `tests/`. Cross-service observable contracts live in `tests/e2e/` and `tests/test_service_boundaries.py`; public import smoke checks live in `tests/test_public_api.py`.

<a id="converter-test-contracts"></a>
## Converter end-to-end test contracts

`tests/e2e/` independently exercises all seven converters and their CLI entrypoints against stored inputs and complete reviewed outputs. The corpus uses reduced public GEO evidence, an unmodified public MAGE-TAB study, and a four-cell six-gene public PBMC count slice, plus explicitly synthetic edge cases. [Fixture provenance and review](../tests/fixtures/README.md) explains exact comparisons, source checks, permitted execution normalization, offline boundaries and the resolved 5-prime/3-prime discrepancy (passing regression MSC-TEST-001).

Run `.venv/bin/python -m pytest tests/e2e -q`; run `.venv/bin/python -m pytest -q` for the complete offline suite. Tests compare full artifacts, exercise injected collaborators, verify old checkpoint bytes survive a version change, and deliberately corrupt outputs to prove scientific drift is rejected. The project-local installed MSC version must match the checkout. No production code or consumer APIs are changed by the test audit.

[Test audit and dispositions](../tests/AUDIT.md) records retained, strengthened, consolidated, relocated and removed checks. Documentation policy lives under `tests/policy/`, separately from workflow correctness. `tests/test_documented_imports.py` verifies executable MSC imports in Python examples.


<a id="scoped-library-chemistry"></a>
## Scoped library chemistry

`magetab.chemistry.resolve_chemistry(sample, channel=None, run=None, series=None)` is a pure service exported from `magetab` with `ChemistryResult`, `ChemistryEvidence` and `ChemistryDiagnostic`. It returns manufacturer, family, version candidates, library role, index configuration, ordered renderable attributes, field-level source evidence and diagnostic codes. It performs no retrieval and changes no input models.

Flow: GEO retrieval/parser/enrichment or JSON package loading → AEConstructor sample/library technology decisions → shared retrieval and protocol registration → per-sample/channel/run SDRF handler → scoped chemistry resolution → ordered SDRF attributes and operation audit → existing IDF/SDRF publication. Retained v2/v3 handler options select rendering, not scientific presets.

The resolver examines extraction protocols, explicit sample/library descriptions, and channel characteristics named `singlecell_type`, `chemistry`, or `library_chemistry` (case and space/underscore/hyphen variants). Exact documented SC3Pv1/v2/v3/v4, SC3Pv3HT, SC5P-PE/R2, SC5P-PE-v3/R2-v3 and SC5PHT identifiers supply family/version candidates, never read geometry. Unknown codes retain raw evidence and an `unknown_identifier` diagnostic; `auto` supplies no identity. A channel-specific identifier narrows compatible version alternatives within one preparation phrase, while independent contradictory preparations remain ambiguous. Narrowing never crosses channels; all original evidence is retained. It normalizes prime spellings and version punctuation only for matching, binds versions to preparation phrases, excludes unsupported software/fixation/compatibility/negated clauses, and keeps original evidence text and paths. Shared overall-design sentences require an explicit sample identifier or an all-samples/all-libraries preparation statement. General summaries and neighbouring samples do not supply chemistry.

Fields resolve independently: 5-prime v1.1/v2 retains 5-prime and ambiguous versions; competing 3-prime and 5-prime preparation phrases retain neither as the end bias. Flex and Multiome are separate families. Reported transcript read lengths, paired barcode/UMI descriptions, explicit offsets and index cycle counts can be emitted; a bare read-length list and vendor-recommended recipes cannot. Multiple labelled library recipes without a matching library identity remain unresolved. Diagnostics are deduplicated within the handler's `last_sdrf_audit.warnings` and include source paths. Existing imported assay-path attributes remain authoritative during semantic overlay.

This is a deliberately bounded grammar, not a general natural-language parser or a vendor layout registry. Complex/unsupported wording can remain unresolved. Existing spatial barcode presets remain outside this fix; spatial rendering also stops assigning cDNA length from an unlabelled list. No raw matrices, schemas, ontology policy or consumer code change.

`tests/magetab/test_chemistry.py` covers resolution, scope, conflicts, unknown kits, read recipes, forced handlers and mixed samples. `tests/e2e/test_chemistry.py` exercises both actual output converters using seven source-excerpt fixtures in explicitly constructed GEO envelopes, compares complete SDRFs, checks diagnostics and reparses repeated chemistry comments. `test_geo_chemistry_regression.py` exercises the original frozen real GEO record. Fixture provenance and expected-value review are in [chemistry fixtures](../tests/fixtures/chemistry/README.md).

Public source symbols:

- `meta_standards_converter.magetab.chemistry.resolve_chemistry`: pure resolver with the signature and scoping contract above.
- `meta_standards_converter.magetab.chemistry.ChemistryResult`: immutable result with ordered attributes, source evidence and diagnostics.
- `meta_standards_converter.magetab.chemistry.ChemistryEvidence`: immutable field, value, source path and original text.
- `meta_standards_converter.magetab.chemistry.ChemistryDiagnostic`: immutable code, affected field and source paths.


<a id="sample-library-routing"></a>
## Sample and library technology routing

`meta_standards_converter.magetab.technology.resolve_technology(sample, channel=None, run=None, *, data=None)` is a pure decision service, exported from `magetab`. It returns a handler key, original evidence with paths, and diagnostics. `detect_ae_technology(data)` supplies an IDF summary key, not a study-wide SDRF assignment.

Evidence precedence is library/channel identity (including documented chemistry identifiers), sample title/description, preparation evidence, general library source, then unambiguous study fallback. Explicit scRNA/snRNA identity therefore survives a shared Visium preparation description. A generic "single cell" source label does not contradict explicit Visium identity. Equally applicable scRNA and spatial identities yield `sequencing` with `ambiguous_technology`, rather than an arbitrary winner. Negated/compatibility clauses, processing text and filenames do not supply method assignments. Study statements naming another GEO sample or other/subset libraries do not apply to the current sample. Unsupported phrasing remains unresolved; this is bounded deterministic matching.

Automatic construction retrieves run evidence once using the operation cache (including sequencing samples inside mixed array/sequencing studies), then refines routing for each ordered sample/channel/run. Homogeneous studies retain their original handler build and audit behavior. Mixed handlers operate on scoped views of the original source objects, share one `ProtocolRegistry` and audit, preregister protocols in source order and collect paths before a single global column plan/render. IDF uses generic sequencing for mixed sequencing technologies, retaining every registered protocol. Explicit `platform_handler` options retain their dispatch contract. Exceptions are not caught or reclassified by dispatch.

Spatial construction names Visium only from applicable sample/library evidence; general shared study descriptions cannot name a particular library. Existing explicitly spatial barcode presets remain unchanged. Imported explicit metadata still overlays generated annotations through the existing semantic retention path. No MINiML schema, consumer API, matrix-processing or live-run change is involved.

Public types and service:

- `meta_standards_converter.magetab.technology.resolve_technology`: pure scoped routing with the signature above.
- `meta_standards_converter.magetab.technology.TechnologyDecision`: immutable handler, ordered evidence and diagnostics.
- `meta_standards_converter.magetab.technology.TechnologyEvidence`: source path, original text and technology candidates.
- `meta_standards_converter.magetab.technology.TechnologyDiagnostic`: diagnostic code and contributing paths.

`tests/magetab/test_scoped_routing.py` verifies the identifier map, conflicts, channel isolation, mixed routing, forced dispatch, shared SRA lookup and protocol references, plus unchanged homogeneous audits. `tests/e2e/test_scoped_routing.py` runs both actual converters on reduced GSM5388031/GSM9254695 records and a constructed mixed study, comparing complete stored IDF/SDRF and audits and reparsing repeated comments and original protocols. See [fixture provenance and scientific expectation review](../tests/fixtures/chemistry/routing/README.md). Scientific deltas were reviewed against isolated source at `6170907`; previous fixture expectations remain unchanged.
