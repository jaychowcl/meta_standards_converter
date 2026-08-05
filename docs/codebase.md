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
- **Decision:** `ae2json` stores a typed `mage_tab.model` plus fingerprints and source tables; `json2ae` restores unchanged tables or overlays eligible mapped edits.
- **Rationale:** Not documented.
- **Consequences:** Lossless unchanged round trips coexist with editable mapped fields, but fingerprint and keyed-union rules are compatibility-sensitive.
- **Affected components:** `AEParser`, `ae_roundtrip`, `ae_model`, and `AEConstructor`.
- **Evidence:** [`ae_parser.py`](../src/meta_standards_converter/ae_handlers/ae_parser.py), [`ae_roundtrip.py`](../src/meta_standards_converter/ae_handlers/ae_roundtrip.py), and [`ae_model.py`](../src/meta_standards_converter/ae_handlers/ae_model.py).

<a id="decision-expression-source-planning"></a>
### AD-003: Select expression assets before normalizing AnnData

- **Status:** Observed
- **Decision:** `json2h5ad` resolves explicit and discovered sources into a plan, then dispatches processed assets or raw nf-core execution before a common normalization/projection stage.
- **Rationale:** Not documented.
- **Consequences:** Source precedence is centralized; incompatible samples remain as per-sample H5ADs and make the aggregate result partial.
- **Affected components:** `AssetManifest`, `SourcePlanner`, `NFCoreRunner`, and `JSON2H5ADConverter`.
- **Evidence:** [`json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py) and [`test_h5ad_pipeline.py`](../tests/test_h5ad_pipeline.py).

<a id="decision-fail-closed-projectors"></a>
### AD-004: Extend emitted metadata through fail-closed projector protocols

- **Status:** Observed
- **Decision:** AnnData and delimited outputs accept injected projector protocols; explicit tabular projectors replace the default MSC projection.
- **Rationale:** Not documented.
- **Consequences:** Downstream adapters can add organization-specific fields without coupling them into MSC, while collisions, invalid vector lengths, and projector errors stop unsafe output unless tabular invalid-output mode is explicit.
- **Affected components:** `AnnDataMetadataProjector`, `TabularMetadataProjector`, H5AD normalization, and delimited converters.
- **Evidence:** [`json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py), [`json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py), and their projector tests.

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
- Unchanged MAGE-TAB round trips remain lossless. Regeneration preserves
  typed-model structure and only overlays mapped fields when identity is
  unambiguous.
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
- Converter-owned nf-core input/output/reference parameters override additional
  params JSON; Nextflow config is for infrastructure and resources.
- Rootless Compose refuses a daemon without rootless security mode and confines
  writable state to the configured output tree.

<a id="proposed-enriched-miniml-core"></a>
## Enriched MINiML-compatible core

`ae2json` validates and publishes MAGE-TAB-only semantics under the namespaced
`mage_tab.model` and preserves exact source tables and fingerprints under
`mage_tab.roundtrip`. An unchanged single-SDRF package can therefore be
restored exactly; edited packages regenerate from the typed model and overlay
only unambiguous mapped core values. Removing `mage_tab` is not lossless.

The fixed MINiML-compatible core has no faithful location for arbitrary
protocol graphs and their performers, hardware, software, parameters, or
accessions; independent assay paths and ordered node/`Protocol REF` chains;
typed attributes with unit and ontology companions; QC, replicate, and
normalization declarations; arbitrary IDF/SDRF properties; or multiple source
document boundaries. Flattening these concepts into existing Series, Sample,
or Platform keys would lose ordering, multiplicity, identity, or
column-occurrence information needed to reconstruct MAGE-TAB.

Schema version 1 contains protocols, assay paths, typed attributes,
declarations, generic properties, and document boundaries. Typed attributes
may carry additive `hz_value*`, `hz_field`, and `hz_unit*` annotations.
`json2ae` publishes these as adjacent reserved `Comment[hz_*]` columns without
replacing raw SDRF cells, and `ae2json` reattaches them on parsing.

The contract is additive: existing core names, types, and meanings remain
unchanged, packages without the extension remain valid, and raw round-trip
evidence is never harmonized. Tabular consumers add sample-bound
`msc.mage_tab.parameter.*` summaries; H5AD also publishes the lossless
occurrence table in `uns["msc_mage_tab"]`.

**Current-state evidence:** [`ae_parser.py`](../src/meta_standards_converter/ae_handlers/ae_parser.py),
[`ae_model.py`](../src/meta_standards_converter/ae_handlers/ae_model.py),
[`ae_roundtrip.py`](../src/meta_standards_converter/ae_handlers/ae_roundtrip.py),
and [`ae_constructor.py`](../src/meta_standards_converter/ae_handlers/ae_constructor.py).

<a id="component-relationships-and-data-flow"></a>
## Component relationships and data flow

| Source | Destination | Interface and direction | Data/lifecycle | Failure behavior |
| --- | --- | --- | --- | --- |
| CLI module | converter | `convert(...)` control call | One converter instance per command; inputs processed in order | Logs per-input exception and records non-zero status |
| GEO converters | `GEOWebFetcher` → `GEOParser` | internal calls around HTTP/XML | GSE becomes Series-scoped package dictionaries | Fetch/parse errors propagate |
| GEO/JSON converters | `MINiMLEnricher` | optional internal mutation | Adds PubMed and SRA/ENA evidence to a package | Enricher records service-specific misses where implemented |
| GEO/JSON converters | `AEConstructor` | `miniml2magetab` then optional `magetab2file` | Package becomes IDF/SDRF row collections and files | Validation/handler/write errors propagate |
| `ae2json` | `AEWebFetcher` → `AEParser` | resolve and parse | IDF/SDRF text becomes package + typed sidecar | Invalid source cardinality or MAGE-TAB fails |
| tabular converters | `JSONPackageSource` → projectors | load/group then `project_sample` | Dataset groups become ordered row maps | Collisions/invalid output fail closed by default |
| `JSON2H5ADConverter` | planner/downloader/runner | plan, localize, or process | Per-conversion state becomes sample AnnData | Per-sample failures retained; aggregate can be partial |
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
- direct Python converter classes, the sixteen formal exports from
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
- `JSONDataOutputOrchestrator` is the public JSON-origin facade. Its manifest,
  H5AD, and AnnData-metadata methods back the three thin CLI wrappers.
  `JSON2H5ADConverter` owns the expression conversion lifecycle beneath it.
  `SourcePlanner`, `AssetManifest`, and `AssetDownloader` resolve inputs;
  `ReferenceResolver`, `AnnotationConverter`, and `NFCoreRunner` own raw-data
  execution; result dataclasses expose complete and partial outcomes.
<a id="orchestrator-json2delimited-converter"></a>
- `JSON2DelimitedConverter` owns JSON grouping, sample iteration, projection,
  column ordering, validation policy, and file output. `JSON2TSVConverter`
  selects TSV or CSV through `output_format`.
<a id="core-rate-limited-requester"></a>
- `RateLimitedRequester` is the shared external-call boundary.
  `GEOWebFetcher`, `AEWebFetcher`, `PubmedWebFetcher`, and `INSDCWebfetcher`
  apply repository-specific URL and response semantics.
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

The formal support boundary is the twenty-one names in
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

- **Signature:** `JSON2H5ADConverter(planner=None, pipeline_runner=None, downloader=None, metadata_projectors=None, package_source=None)`; `convert(...)` and `convert_source(...)` own the documented expression workflow.
- Dataset, study, sample, checkpoint, and output identities are validated as safe single path components before publication. `series.iid` is the canonical native-package study identity, so an ArrayExpress IID is not displaced by an earlier GEO secondary accession.
- `DatasetBundleRecoveryError` preserves the original publication error, rollback errors, and surviving recovery paths when an overwrite cannot be fully restored.
- **Inputs:** native MINiML or Atlas v1 JSON, source/reference/runtime options,
  and optional public collaborators.
- **Outputs:** single or batch conversion results plus a transactional H5AD
  artifact bundle.
- **Failures:** structural, validation, filesystem, external process, and
  publication failures follow the strict/permissive and batch contracts in
  [JSON to H5AD](#workflow-json2h5ad).
- **Side effects:** may read/download assets, run nf-core, and publish staged
  H5AD/provenance artifacts.
- **Support:** formal export and composition boundary.
- **Source:** [`converters/json2h5ad.py`](../src/meta_standards_converter/converters/json2h5ad.py).

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

- **Signature:** structural `Protocol` with `project_sample(*, adata, context) -> AnnDataMetadataProjection` and `project_combined(*, adata, contexts) -> AnnDataMetadataProjection`.
- **Inputs:** current AnnData plus one context, or combined AnnData plus a context sequence.
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
### `JSONDataOutputOrchestrator`

- **Signature:** `JSONDataOutputOrchestrator(tabular_projectors=None, h5ad_converter=None)` with `export_manifest`, `export_h5ad`, and `export_anndata_metadata`.
- **Inputs:** parsed MINiML/Atlas JSON, output directory, format/component controls, and the H5AD asset/reference/pipeline options.
- **Outputs:** `TabularConversionResult`, H5AD results, or `AnnDataMetadataExportResult`/`AnnDataMetadataBatchResult`.
- **Failures:** source, validation, asset, pipeline, serialization, collision, and atomic-publication failures propagate or enter batch failures.
- **Side effects:** publishes operation-owned artifact bundles and JSON result manifests.
- **Support:** formal export and preferred Python interface for JSON-origin outputs.
- **Source:** [`converters/json_outputs.py`](../src/meta_standards_converter/converters/json_outputs.py).

<a id="api-anndata-metadata-export-result"></a>
### `AnnDataMetadataExportResult`

- **Signature:** `AnnDataMetadataExportResult(dataset_id, obs, var, uns, obs_path, var_path, uns_path, manifest_path, warnings=(), errors=(), partial=False)`.
- **Inputs:** one dataset's combined typed AnnData metadata, optional published paths, and diagnostics.
- **Outputs:** in-memory `obs`, optional `var`/`uns`, artifact locations, and a compact `to_dict()` summary.
- **Failures:** construction performs no custom validation; the orchestrator validates and serializes its data.
- **Side effects:** none.
- **Support:** formal export and successful per-dataset result for `json2obs`.
- **Source:** [`converters/json_outputs.py`](../src/meta_standards_converter/converters/json_outputs.py).

<a id="api-anndata-metadata-batch-result"></a>
### `AnnDataMetadataBatchResult`

- **Signature:** `AnnDataMetadataBatchResult(source, conversions=(), warnings=(), failures=())`.
- **Inputs:** source path, completed dataset exports, cross-dataset warnings, and keyed failures.
- **Outputs:** immutable batch status, `partial` state, and compact `to_dict()` summary.
- **Failures:** construction performs no custom validation; conversion failures are retained in `failures`.
- **Side effects:** none.
- **Support:** formal export and multi-dataset result for `json2obs`.
- **Source:** [`converters/json_outputs.py`](../src/meta_standards_converter/converters/json_outputs.py).

<a id="api-json2tsv-converter"></a>
### `JSON2TSVConverter`

- **Signature:** `JSON2TSVConverter(metadata_projectors=None, package_source=None, *, output_format="tsv")`; inherited `convert_source(source, destination, *, allow_invalid=False, overwrite=False) -> TabularConversionResult`.
- **Inputs:** parsed MINiML JSON or canonical Atlas v1, TSV/CSV format, and destination.
- **Outputs:** selected delimited file and result metadata.
- **Failures:** invalid formats, source, projector, collision, fail-closed diagnostic, and protected-output errors propagate.
- **Side effects:** creates the destination parent and writes TSV or CSV.
- **Support:** formal export.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

<a id="api-msc-metadata-projector"></a>
### `MSCMetadataProjector`

- **Signature:** `MSCMetadataProjector()` and `project_sample(*, context: TabularMetadataContext) -> TabularMetadataProjection`.
- **Inputs:** one sample context.
- **Outputs:** canonical dotted `msc.*` values and ordered base columns.
- **Failures:** converter validation applies to the returned projection.
- **Side effects:** none.
- **Support:** formal export and default projector when no explicit projectors are supplied.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

<a id="api-tabular-conversion-result"></a>
### `TabularConversionResult`

- **Signature:** `TabularConversionResult(row_count, columns, dataset_ids, warnings=(), errors=(), output_path=None)`.
- **Inputs:** immutable output summary values.
- **Outputs:** frozen result; `partial` is `True` exactly when `errors` is non-empty.
- **Failures:** no custom validation.
- **Side effects:** none.
- **Support:** formal export; `partial` is its public property.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

<a id="api-tabular-metadata-context"></a>
### `TabularMetadataContext`

- **Signature:** `TabularMetadataContext(package, sample, dataset_id, study_accession, sample_accession, base_metadata)`.
- **Inputs:** group/package/sample mappings, identifiers, and base metadata.
- **Outputs:** frozen context passed to every tabular projector.
- **Failures:** no custom validation.
- **Side effects:** none.
- **Support:** formal export; these are all current fields.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

<a id="api-tabular-metadata-projection"></a>
### `TabularMetadataProjection`

- **Signature:** `TabularMetadataProjection(values, columns=(), warnings=(), errors=())`.
- **Inputs:** projected values, preferred column order, and diagnostics.
- **Outputs:** frozen per-sample projection.
- **Failures:** converter rejects wrong types/collisions; errors raise `TabularProjectionError` unless `allow_invalid=True`.
- **Side effects:** none.
- **Support:** formal export.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

<a id="api-tabular-metadata-projector"></a>
### `TabularMetadataProjector`

- **Signature:** structural `Protocol` with `project_sample(*, context: TabularMetadataContext) -> TabularMetadataProjection`.
- **Inputs:** one immutable sample context.
- **Outputs:** projected values, order, warnings, and errors.
- **Failures:** projector exceptions propagate; converter validates type and collisions.
- **Side effects:** none required.
- **Support:** formal export/injection extension point; not runtime-checkable.
- **Source:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py).

Supported production symbols are inventoried below by qualified name. Their
signatures, inputs, outputs, exceptions, side effects, and decisive internal or
external calls are documented in the linked legacy callable sections that
follow this canonical overview.

- MAGE-TAB:
  `meta_standards_converter.ae_handlers.ae_constructor.validate_platform_handler`,
  `meta_standards_converter.ae_handlers.ae_constructor.ProtocolRegistry`,
  `meta_standards_converter.ae_handlers.ae_constructor.AEConstructor`,
  `meta_standards_converter.ae_handlers.ae_idf_handlers.IDFConstructor`,
  `meta_standards_converter.ae_handlers.ae_model.MAGETabModelError`,
  `meta_standards_converter.ae_handlers.ae_model.build_model`,
  `meta_standards_converter.ae_handlers.ae_model.validate_model`,
  `meta_standards_converter.ae_handlers.ae_model.render_model`,
  `meta_standards_converter.ae_handlers.ae_model.overlay_core`,
  `meta_standards_converter.ae_handlers.ae_parser.normalized_label`,
  `meta_standards_converter.ae_handlers.ae_parser.AEParser`,
  `meta_standards_converter.ae_handlers.ae_roundtrip.semantic_sha256`,
  `meta_standards_converter.ae_handlers.ae_roundtrip.model_sha256`,
  `meta_standards_converter.ae_handlers.ae_roundtrip.build_roundtrip`,
  `meta_standards_converter.ae_handlers.ae_roundtrip.unchanged_magetab`,
  `meta_standards_converter.ae_handlers.ae_roundtrip.restore_extensions`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFAttr`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFNode`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFEdge`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFPath`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.ColumnGroup`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFAudit`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.SDRFConstructor`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.normalized_extension`,
  `meta_standards_converter.ae_handlers.ae_sdrf_handlers.classify_file`,
  `meta_standards_converter.ae_handlers.ae_webfetcher.TextResource`,
  `meta_standards_converter.ae_handlers.ae_webfetcher.MAGETabInput`, and
  `meta_standards_converter.ae_handlers.ae_webfetcher.AEWebFetcher`.
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
  `meta_standards_converter.converters.ae2json.ae2json`,
  `meta_standards_converter.converters.geo2ae.geo2ae`,
  `meta_standards_converter.converters.geo2json.geo2json`,
  `meta_standards_converter.converters.json2ae.json2ae`,
  `meta_standards_converter.converters.json2h5ad.Asset`,
  `meta_standards_converter.converters.json2h5ad.MetadataProjectionContext`,
  `meta_standards_converter.converters.json2h5ad.AnnDataMetadataProjection`,
  `meta_standards_converter.converters.json2h5ad.AnnDataProjectionError`,
  `meta_standards_converter.converters.json2h5ad.AnnDataMetadataProjector`,
  `meta_standards_converter.converters.json2h5ad.AssetManifest`,
  `meta_standards_converter.converters.json2h5ad.AssetDownloader`,
  `meta_standards_converter.converters.json2h5ad.PipelineRun`,
  `meta_standards_converter.converters.json2h5ad.ConversionResult`,
  `meta_standards_converter.converters.json2h5ad.BatchConversionResult`,
  `meta_standards_converter.converters.json2h5ad.RawProcessingResult`,
  `meta_standards_converter.converters.json2h5ad.ReferenceResolver`,
  `meta_standards_converter.converters.json2h5ad.AnnotationConverter`,
  `meta_standards_converter.converters.json2h5ad.NFCoreRunner`,
  `meta_standards_converter.converters.json2h5ad.SourcePlanner`,
  `meta_standards_converter.converters.json2h5ad.JSON2H5ADConverter`, and
  `meta_standards_converter.converters.json2h5ad.json2h5ad`.
- Tabular and JSON source:
  `meta_standards_converter.converters.json2tabular.TabularMetadataContext`,
  `meta_standards_converter.converters.json2tabular.TabularMetadataProjection`,
  `meta_standards_converter.converters.json2tabular.TabularMetadataProjector`,
  `meta_standards_converter.converters.json2tabular.TabularConversionResult`,
  `meta_standards_converter.converters.json2tabular.TabularProjectionError`,
  `meta_standards_converter.converters.json2tabular.MSCMetadataProjector`,
  `meta_standards_converter.converters.json2tabular.JSON2DelimitedConverter`,
  `meta_standards_converter.converters.json2tabular.JSON2TSVConverter`,
  `meta_standards_converter.converters.json2tabular.json2tsv`,
  `meta_standards_converter.converters.json_outputs.JSONDataOutputOrchestrator`,
  `meta_standards_converter.converters.json_outputs.AnnDataMetadataExportResult`,
  `meta_standards_converter.converters.json_outputs.AnnDataMetadataBatchResult`,
  `meta_standards_converter.converters.json_source.DatasetPackageGroup`,
  `meta_standards_converter.converters.json_source.SourceLoadResult`, and
  `meta_standards_converter.converters.json_source.JSONPackageSource`.
- Fetch, parse, enrich, harmonize, and helpers:
  `meta_standards_converter.enrichers.miniml_enricher.MINiMLEnricher`,
  `meta_standards_converter.geo_handlers.geo_parser.GEOParser`,
  `meta_standards_converter.geo_handlers.geo_webfetcher.GEOWebFetcher`,
  `meta_standards_converter.harmonizers.geo2ols.GEO2OLS`,
  `meta_standards_converter.harmonizers.harmonizers.Harmonizer`,
  `meta_standards_converter.harmonizers.pubmed2ols.Pubmed2OLS`,
  `meta_standards_converter.helpers.json_helper.JSONHandler`,
  `meta_standards_converter.helpers.request_helper.RequestSettings`,
  `meta_standards_converter.helpers.request_helper.RateLimitedRequester`,
  `meta_standards_converter.insdc_handlers.insdc_webfetcher.INSDCWebfetcher`,
  `meta_standards_converter.meta_store.meta_store.MetaStore`, and
  `meta_standards_converter.pubmed_handlers.pubmed_webfetcher.PubmedWebFetcher`.

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

**Evidence:** [`converters/geo2ae.py`](../src/meta_standards_converter/converters/geo2ae.py), [`cli/geo2ae.py`](../src/meta_standards_converter/cli/geo2ae.py), and [`geo_webfetcher.py`](../src/meta_standards_converter/geo_handlers/geo_webfetcher.py).

<a id="workflow-geo2json"></a>
### `geo2json`: GEO to parsed JSON

```text
GSE -> fetch --failure--> exception
    -> parse --failure--> exception
    -> [enrich?] yes -> enrich each --failure--> exception
                  no --------------------------> retain parsed packages
    -> [out?] JSON file / else return only
```

1. The CLI calls `geo2json.convert` once per accession and continues after failures.
2. GEO fetch and parsing are shared with `geo2ae`.
3. `enrich=False` bypasses PubMed/SRA enrichment.
4. `json2file` creates the output directory and writes `{GSE}.json` when requested.

Pseudocode: `packages = parse(fetch(gse)); [enrich packages]; [write]; return`.

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
6. Round-trip evidence may restore source tables; mapped edits use overlay
   rules. Enriched `hz_*` attributes become adjacent reserved Comment columns
   while raw value/unit cells remain unchanged. Additive characteristic
   columns are ontology-agnostic, including ECTO `hz_exposure_name*` and PCL
   `hz_cell_state_name*` annotations.
7. `out` controls writing; construction errors propagate.

Pseudocode: `validate(flatten(source.load(path))); warn(skipped); for package:
[enrich] -> construct -> [write]; return`.

**Evidence:** [`converters/json2ae.py`](../src/meta_standards_converter/converters/json2ae.py), [`ae_constructor.py`](../src/meta_standards_converter/ae_handlers/ae_constructor.py), and [`ae_roundtrip.py`](../src/meta_standards_converter/ae_handlers/ae_roundtrip.py).

<a id="workflow-ae2json"></a>
### `ae2json`: MAGE-TAB to parsed JSON

```text
local/HTTP/accession -> resolve IDF + SDRF(s) --failure--> exception
                     -> parse/model/round-trip --failure--> exception
                     -> [out?] accession JSON / else return package list
```

1. `AEWebFetcher.resolve` accepts a local/HTTP IDF or BioStudies accession and optional SDRF overrides.
2. Resolution requires exactly one IDF and at least one SDRF.
3. `AEParser.parse` maps core fields and retains typed/source round-trip evidence.
4. `out` writes a sanitized accession filename; otherwise no file is created.

Pseudocode: `resolved = fetcher.resolve(source); package = parser.parse(resolved); [write]; return [package]`.

**Evidence:** [`converters/ae2json.py`](../src/meta_standards_converter/converters/ae2json.py), [`ae_webfetcher.py`](../src/meta_standards_converter/ae_handlers/ae_webfetcher.py), and [`ae_parser.py`](../src/meta_standards_converter/ae_handlers/ae_parser.py).

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
          -> per-sample failure retained; successes continue
          -> sample H5AD -> [compatible?] combined H5AD
          -> provenance manifest / partial result
```

1. `convert(json_path, ...)` returns `ConversionResult` for one group or
   `BatchConversionResult` for multiple groups.
2. `convert_source(json_path, ...)` always returns `BatchConversionResult`,
   including for one group.
3. With multiple groups, `_convert_groups` places each group in an output-root
   child directory named for `dataset_id`; one group uses the root directly.
4. Planning applies manifest, explicit, then discovered asset precedence.
5. Processed assets normalize directly; ordinary H5AD and delimited paths do
   not import Scanpy, while 10x HDF5/MTX branches import it lazily. Raw assets
   call reference resolution and Nextflow/nf-core through `NFCoreRunner`.
6. Per-sample failures, incompatibility, and allowed projector errors can make `ConversionResult.partial`;
   per-group exceptions are caught in `BatchConversionResult.failures`.
7. Invalid path/source/no-group/no-sample and unsafe dataset-ID conditions raise before aggregation; successful
   groups and diagnostics survive later failures.
8. Bound MAGE-TAB Parameter Values are projected to dotted `obs` columns and
   every typed attribute occurrence is retained in `uns["msc_mage_tab"]`.
9. `--processed-checkpoint-dir` writes each normalized, projected sample H5AD
   atomically with a fingerprint over the source JSON, sample, asset,
   orientation, and converter version. With `--resume`, matching checkpoints
   are loaded instead of downloading and converting that sample again;
   incomplete, corrupt, or stale checkpoint pairs are ignored.

Pseudocode: `load -> if one and convert: convert_packages; else for group: try convert_packages into child; except record; return batch`.

**Evidence:** [`JSON2H5ADConverter`](../src/meta_standards_converter/converters/json2h5ad.py), [`JSONPackageSource`](../src/meta_standards_converter/converters/json_source.py), and [`cli/json2h5ad.py`](../src/meta_standards_converter/cli/json2h5ad.py).

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

**Evidence:** [`converters/json2tabular.py`](../src/meta_standards_converter/converters/json2tabular.py), [`converters/json_source.py`](../src/meta_standards_converter/converters/json_source.py), and [`cli/json2tsv.py`](../src/meta_standards_converter/cli/json2tsv.py).

<a id="workflow-json2obs"></a>
### `json2obs`: MINiML/Atlas JSON to AnnData metadata

```text
JSON + expression assets -> shared H5AD assembly
       -> combined AnnData unavailable -------------------> failure
       -> obs with named cell_id --------------------------> .obs.csv
       -> include_var -------------------------------------> .var.csv
       -> include_uns -------------------------------------> typed .uns.json
       -> atomic component bundle + JSON result manifest
```

1. Asset resolution, raw processing, normalization, projectors, smart IDs, and combination exactly match `json2h5ad`.
2. The combined `obs` is exported with `cell_id`; source and canonical names remain unchanged.
3. Optional `var` uses `feature_id`; optional `uns` uses tagged JSON for nested mappings, arrays, and DataFrames.
4. Every batch conversion publishes beneath its dataset-ID directory, even when other groups fail and only one dataset succeeds.
5. The CLI prints only a compact JSON summary to stdout, sends logs to stderr, and returns status `1` for partial/failure outcomes.

Pseudocode: `assemble -> read combined AnnData -> serialize selected components -> atomic publish -> result`.

**Evidence:** [`converters/json_outputs.py`](../src/meta_standards_converter/converters/json_outputs.py) and [`cli/json2obs.py`](../src/meta_standards_converter/cli/json2obs.py).

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

<a id="project-purpose-and-layout"></a>
## Project Purpose And Layout

`meta_standards_converter` converts biological study metadata between repository standards. The current package focuses on GEO MINiML to ArrayExpress/MAGE-TAB-style output.

```text
src/meta_standards_converter/
├── cli/
│   ├── common.py                 # shared CLI logging helpers
│   ├── geo2ae.py                 # geo2ae command-line entrypoint
│   ├── geo2json.py               # geo2json command-line entrypoint
│   ├── json2ae.py                 # parsed JSON-to-MAGE-TAB command-line entrypoint
│   ├── ae2json.py                 # MAGE-TAB-to-JSON command-line entrypoint
│   ├── json2h5ad.py              # multi-source JSON-to-H5AD command-line entrypoint
│   ├── json2tsv.py               # JSON/Atlas-to-TSV-or-CSV manifest entrypoint
│   └── json2obs.py               # JSON/Atlas plus data assets to obs/var/uns entrypoint
├── converters/
│   ├── geo2ae.py                 # top-level GEO to AE orchestration
│   ├── geo2json.py               # top-level GEO to JSON orchestration
│   ├── json2ae.py                 # parsed JSON validation and AE orchestration
│   ├── ae2json.py                 # MAGE-TAB resolution and JSON orchestration
│   ├── json2h5ad.py              # asset planning, AnnData conversion, and nf-core orchestration
│   ├── json2tabular.py           # injectable TSV/CSV projection orchestration
│   ├── json_outputs.py           # shared manifest/H5AD/obs output orchestration
│   └── json_source.py            # MINiML and Atlas v1 package grouping
├── atlas_v1/
│   └── reader.py                 # standalone Atlas v1 validation and adaptation
├── geo_handlers/
│   ├── geo_webfetcher.py         # GEO MINiML URL building and download
│   └── geo_parser.py             # MINiML XML to JSON-ready per-Series packages
├── pubmed_handlers/
│   └── pubmed_webfetcher.py      # PubMed ESummary lookup and parsed publication metadata
├── enrichers/
│   └── miniml_enricher.py        # Adds PubMed/SRA records to parsed MINiML JSON
├── ae_handlers/
│   ├── ae_idf_handlers.py        # IDF row construction
│   ├── ae_constructor.py         # MAGE-TAB coordination, protocol registry, file writing
│   ├── ae_parser.py              # IDF/SDRF to MINiML-compatible package mapping
│   ├── ae_roundtrip.py           # source-table sidecar, fingerprint, and restoration
│   ├── ae_webfetcher.py          # local, HTTP, and BioStudies MAGE-TAB resolution
│   └── ae_sdrf_handlers.py       # SDRF graph model and technology handlers
├── harmonizers/
│   ├── geo2ols.py                # GEO protocol type ontology mapping
│   ├── pubmed2ols.py             # PubMed status ontology mapping
│   └── harmonizers.py            # combined harmonizer
├── helpers/
│   ├── json_helper.py            # dotted-path JSON helpers
│   └── request_helper.py         # service-specific rate limiting and retries
├── insdc_handlers/
│   └── insdc_webfetcher.py       # SRA accession extraction, NCBI lookup, and run parsing
└── meta_store/
    └── meta_store.py             # placeholder metadata validation store
```

Tests cover parser packaging, converter orchestration, CLI flags, AE constructor composition, IDF behavior, and SDRF rendering:

```text
tests/test_geo_parser.py
tests/test_geo2ae.py
tests/test_geo2json.py
tests/test_json2ae.py
tests/test_ae2json.py
tests/test_ae_webfetcher.py
tests/test_json2h5ad.py
tests/test_cli_geo2ae.py
tests/test_cli_geo2json.py
tests/test_cli_json2ae.py
tests/test_cli_ae2json.py
tests/test_cli_json2h5ad.py
tests/test_project_scripts.py
tests/test_ae_constructor.py
tests/test_ae_sdrf_handlers.py
tests/test_miniml_enricher.py
tests/test_request_helper.py
tests/test_geo_webfetcher.py
tests/test_insdc_webfetcher.py
tests/test_pubmed_webfetcher.py
tests/GSE328265_family.xml
```

<a id="runtime-behavior"></a>
## Runtime Behavior

- Distribution version `1.0.0` is the initial unified JSON-output release. It uses
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
- The `geo2ae`, `geo2json`, `json2ae`, `ae2json`, `json2h5ad`, `json2tsv`, and `json2obs` console scripts point to their matching modules under `meta_standards_converter.cli`.
- Network calls are owned by platform fetchers and routed through `RateLimitedRequester`: `GEOWebFetcher` handles GEO FTP MINiML tarballs and related-series traversal, `AEWebFetcher` handles BioStudies discovery and HTTP(S) MAGE-TAB text, `INSDCWebfetcher` handles NCBI SRA EFetch plus ENA Portal file reports, and `PubmedWebFetcher` handles NCBI PubMed ESummary publication metadata.
- Default request settings are selected per service but enforced across the
  process by normalized hostname: `ncbi_eutils` uses timeout 30s, delay 0.5s,
  at most two in flight, and 3 retries; `geo_ftp`, `biostudies`, and
  `ena_portal` use timeout 30s, delay 1.0s, at most two in flight, and 3 retries.
- Library logging propagates safe structured telemetry to caller handlers.
  DEBUG records service/host, attempt, status, timeout, and duration without URL
  queries or request parameters. INFO records retries, GEO fetch sizes/duration,
  MINiML structural counts, related-series progress, and enrichment hit/failure
  totals. XML, parsed metadata, publication content, tokens, and credentials are
  never logged.
- `geo2ae.convert()` keeps parsed and enriched GEO metadata in memory for MAGE-TAB construction.
- `geo2json.convert()` returns parsed GEO package JSON, enriched by default, and can write `{accession}.json`.
- `json2ae.convert()` loads one parsed package object or a non-empty package list, enriches it by default, and returns or writes MAGE-TAB outputs.
- `ae2json.convert()` resolves one IDF and one or more SDRFs, returns one MINiML-compatible package in a list, and can write `{accession}.json`.
- `json2h5ad.convert()` selects per-sample H5AD, matrix, or raw FASTQ sources; normalizes them into AnnData; and writes per-sample plus compatible combined H5AD outputs.
- `json2tsv --format {tsv,csv}` emits the neutral MSC sample projection; `json2obs` exports the combined AnnData metadata view. Both are orchestrator methods with CLI wrappers.
- When `out` is supplied, `geo2ae.convert()` writes `{accession}.idf.txt` and `{accession}.sdrf.txt`.
- `geo2ae` `out` controls MAGE-TAB output only; use `geo2json` for parsed JSON snapshots.
- Processed `json2h5ad` conversion requires the `h5ad` extra. Raw processing directly on the host additionally requires Nextflow, Java, and a supported execution profile/runtime. The project image includes Java 21, pinned Nextflow, the Docker CLI, and `.[h5ad]`.
- `MetaStore._validate_investigation_metadata_structure()` is a `pass` placeholder, so `validate_investigation_metadata()` currently asserts for normal input.

<a id="end-to-end-geo2ae-flow"></a>
## End-To-End geo2ae Flow

```text
main(argv)
  -> parse CLI args
  -> instantiate geo2ae()
  -> for each GSE accession:
       geo2ae.convert(gse, related_series, remove_empty, out, platform_handler)
       continue to the next accession if a conversion fails
  -> return 1 if any accession failed, else 0

geo2ae.convert(gse, related_series, remove_empty, out, platform_handler=None)
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

`geo2json.convert(gse, related_series, remove_empty, enrich, out)` follows the same GEO fetch and parse stages, optionally enriches each parsed package through `MINiMLEnricher`, writes `{gse}.json` when `out` is truthy, and returns the list of JSON packages without invoking AE/MAGE-TAB construction.

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
  -> instantiate json2ae()
  -> for each JSON path:
       json2ae.convert(json_path, enrich, out, platform_handler)
       continue to the next path if conversion fails
  -> return 1 if any path failed, else 0

json2ae.convert(json_path, out, enrich=True, platform_handler=None)
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
  -> parse one or more sources, optional repeated --sdrf, --out, and logging flags
  -> require exactly one source when --sdrf is present
  -> for each source:
       ae2json.convert(source, out, sdrf_sources)
       continue to the next source if conversion fails
  -> return 1 if any source failed, else 0

ae2json.convert(source, out=None, sdrf_sources=None)
  -> AEWebFetcher.resolve(source, sdrf_sources)
       existing path: read the IDF and relative/local/HTTP SDRF references
       HTTP(S) URL: fetch IDF text and resolve relative/HTTP SDRF references
       accession: paginate BioStudies files, query study info, then fetch IDF/SDRF text
  -> AEParser.parse(resolved_input)
       parse the IDF and rectangular SDRF tables
       map known investigation, publication, contributor, protocol, sample,
       platform, factor, characteristic, SRA/FASTQ, and array-file metadata
       merge repeated sample/platform records in first-seen order
       keep the first conflicting scalar and record a warning
       identify the root as magetabv1.1 with the BioStudies specification URL
       set series.iid from the primary ArrayExpress/investigation accession
       preserve unmapped metadata, typed MAGE-TAB entities, and source tables under mage_tab
       fingerprint the public package fields for unchanged-table detection
  -> if out, write [{package}] to {study_accession}.json
  -> return [package]
```

IDF labels are matched case- and whitespace-insensitively. The parser accepts general MAGE-TAB inputs rather than only files emitted by this project. Its output uses the existing MINiML-compatible top-level shape (`database`, `organization`, `contributor`, `platform`, `sample`, and `series`), plus a namespaced `mage_tab` extension containing source/version metadata, unmapped data, warnings, an editable typed model, and an optional lossless round-trip sidecar. AE-origin packages use root `version = "magetabv1.1"`, root `schema_location = "https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf"`, and an ArrayExpress/investigation accession as `series.iid`; `mage_tab.version` separately retains the original IDF value such as `1.1`. An unchanged one-SDRF package reproduces its parsed source rows exactly. Typed-model edits regenerate MAGE-TAB without merging unsupported values into MINiML fields; exact core edits overlay the corresponding modeled output.

Accession resolution calls `GET /api/v1/files/{accession}` to discover exactly one IDF and at least one SDRF, calls `GET /api/v1/studies/{accession}/info` for the HTTP base, and downloads only those metadata files beneath `Files/`. Remote content is decoded as UTF-8 with optional BOM and remains in memory. FTP sources and referenced assay data downloads are not supported.

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
| Protocols | Native GEO protocol elements | Recognized descriptions in the core; complete independent records in `mage_tab.model.protocols` |
| Platform/contributors | Referenced GEO records with their available MINiML fields | IDF/SDRF projection, usually a smaller record |
| PubMed/SRA enrichment | Enabled by default and optional with `enrich=False`/`--no-enrich` | No remote enrichment stage; publication, ENA, run, and file metadata come from IDF/SDRF fields |
| Empty fields | Removed by default or retained with `remove_empty=False`/`--keep-empty` | Omitted unless a mapped source value exists; required extension structure remains present |
| Assay-row multiplicity | MINiML samples plus optionally enriched `sra_run` lists | Core samples/runs may consolidate rows; every original SDRF row remains an independent `mage_tab.model.assay_paths` record |
| Unsupported metadata | Remains available when it exists as an XML element/attribute | Stored independently in the typed model, unmapped lists, and raw round-trip tables |
| Lossless MAGE-TAB round trip | Not applicable; no `mage_tab` extension | `mage_tab.model` is editable and `mage_tab.roundtrip` retains exact source IDF/SDRF tables and fingerprints |

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
  "version": "magetabv1.1",
  "schema_location": "https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf",
  "series": {"iid": "E-MTAB-1", "accession": [{"value": "E-MTAB-1", "database": "ArrayExpress"}]},
  "sample": [],
  "platform": [],
  "contributor": [],
  "mage_tab": {
    "version": "1.1",
    "model": {"schema_version": 1, "protocols": [], "declarations": {}, "assay_paths": []},
    "roundtrip": {"schema_version": 1, "semantic_sha256": "...", "model_sha256": "..."}
  }
}
```

The shared core makes downstream processing reusable; it does not imply field-for-field parity between repositories. Consumers that need complete MAGE-TAB semantics must retain `mage_tab`, while consumers using only common study/sample metadata can read the core fields from either converter.

<a id="json2h5ad-flow"></a>
## End-To-End json2h5ad Flow

```text
json2h5ad.convert(json_path, out, asset_manifest, asset_specs, force_reprocess, ...)
  -> load ordinary parsed MINiML JSON or validate a canonical Atlas v1 document
  -> group packages by dataset; fail if no convertible groups
  -> one group: convert directly
  -> multiple groups: convert each below out/{dataset_id}, recording group exceptions
  -> AssetManifest loads explicit CSV/TSV and CLI mappings
  -> SourcePlanner discovers sample/study assets and selects per sample:
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
  -> write one normalized H5AD per sample
  -> combine compatible samples using an outer sparse feature join
  -> invoke ordered combined-study metadata projector callbacks
  -> write optional combined study H5AD and JSON provenance manifest
  -> return ConversionResult for one group or BatchConversionResult for multiple
```

`AssetDownloader` streams HTTP(S)/FTP processed assets into an output-local cache and verifies an MD5 when supplied. Source files are never modified. Gzip-compressed H5AD assets are expanded into a temporary `.h5ad` only while AnnData reads them; the cached download remains compressed. Study-level H5AD splitting recognizes canonical `msc.sample.accession` and the external generic columns `geo_accession`, `sample_id`, `sample`, and `gsm_accession`.

<a id="h5ad-metadata-schema-v3"></a>
<a id="h5ad-metadata-schema-v1"></a>
### H5AD metadata schema 1.0

Converter-owned observation metadata uses only dotted names grouped under `msc.sample`, `msc.series`, `msc.platform`, `msc.archive`, `msc.library`, `msc.instrument`, `msc.protocol`, `msc.database`, `msc.asset`, `msc.expression`, `msc.characteristics`, `msc.observation`, and `msc.combination`. The converter does not generate underscore aliases. Existing underscore-style columns from an input H5AD remain opaque source columns: normalization preserves but neither interprets nor validates them. Custom projectors retain ownership of their injected names.

Stable observation fields cover sample/study accessions, title/description, organism and taxid, organism part, developmental stage, disease, genotype, biological source, material/provider/molecule, platform, SRA/ENA/BioSample/run accessions, library fields, instrument, modality, asset provenance, and database identity. Organism resolution evaluates each channel independently, preferring non-empty `hz_organism` values before raw `organism`; original taxids remain separate. Database identity uses `public_id`, then `iid`, then `name`. Every characteristic becomes `msc.characteristics.<normalized_tag>`, including harmonized `hz_*` tags. Missing values are empty in `obs`; repeated values are case-insensitively de-duplicated in source order and displayed with `; ` separators.

`uns["msc_metadata"]` declares schema version `1.0` and contains the authoritative normalized `sample_values` DataFrame with `sample_accession`, `field`, `ordinal`, `value`, and `value_type`. It stores one row per non-empty canonical value, so embedded semicolons and list cardinality remain recoverable without parsing the display string. Sample H5ADs contain their sample rows; combined H5ADs contain every sample plus `msc.combination.batch`. `uns["msc_miniml"]` remains the complete typed source ledger at schema 1.0. H5AD provenance and manifests separately declare the H5AD metadata schema version.

Normalization copies each incoming index into `msc.observation.original_id`. An identifier is already sample-qualified when its accession occurs case-insensitively as a token bounded by the start/end or `-`, `_`, `.`, or `:`. Qualified identifiers are preserved; other identifiers receive `-{sample_accession}`. Repeated candidates receive source-order numeric suffixes. A final cross-sample pass qualifies any remaining collision before concatenation and fails if uniqueness cannot be established. Sample H5ADs and the combined H5AD therefore expose identical, globally unique observation identifiers.

`uns["msc_miniml"]` contains schema/policy metadata, the source JSON path and SHA-256, and a typed long-form `fields` DataFrame (`package_index`, `entity_type`, `entity_id`, `path`, `value`, `value_type`). GSM files contain the sample plus its series and transitively referenced platform, contributor, and database records without following `sample_ref`; both scalar references and real MINiML `{"ref": "..."}` objects are resolved. GSE files contain all package entities. Protocol descriptions remain in this table, while `msc.protocol.types`, source refs, and accessions reuse `Harmonizer.geoprotocols2efo()` for the established treatment, growth, extraction, labeling, hybridization, scan, and data-processing paths. Publication records are whitelisted to PubMed ID, DOI, title, authors, status, and status ontology fields; abstracts, full text, article bodies, sections, and other publication content are not embedded. GEO series summary and overall design remain experiment metadata.

Persisted local provenance paths in the manifest and H5AD metadata are relative to the containing artifact and declare `path_base = artifact_parent`; parallel scope fields distinguish internal, external (`../...`), and remote locations. Remote URLs remain unchanged. Recorded Nextflow command path arguments are also relative, while generated runtime configs retain the absolute paths required by Nextflow resume. `ConversionResult` continues to return absolute paths in memory. ANSI-stripped, de-duplicated Nextflow warnings are stored on each pipeline run and promoted to the manifest's top-level warnings without changing a successful return code.

Before writing nf-core samplesheets, `NFCoreRunner` upgrades raw `ftp://` URLs from the known ENA and NCBI archive hosts to their equivalent `https://` endpoints. This avoids truncated Java FTP transfers through rootless container networking while leaving unknown FTP servers unchanged.

Raw processing pins `nf-core/scrnaseq` 4.2.0 and `nf-core/rnaseq` 3.26.0 by default. The scrnaseq 4.2.0 floor includes the upstream strict-syntax fixes required by the pinned Nextflow 26 runtime; 4.1.0 contains a reference to a missing `conf/test_multiome.config` and fails during config parsing. The runner requires Nextflow, Java, and the selected Docker/Podman/Apptainer/Singularity runtime. `scrnaseq` discovers H5ADs across the results tree, associates outputs when the sample accession appears in the filename or containing directories, and prefers CellBender-filtered, then filtered, then raw output. This includes QCATCH names such as `GSM1_filtered_quants.h5ad`, not only `*_matrix.h5ad`. For each sample, `rnaseq` selects its named numeric column from the merged count and TPM tables, tolerates the standard text `gene_name` column, writes counts to sparse `X`, writes aligned TPM values to `layers["tpm"]`, and preserves `gene_name` in `var`.

When `META_STANDARDS_REQUIRE_ROOTLESS_DOCKER` is truthy and the Docker profile is selected, `NFCoreRunner._preflight()` queries `docker info` before creating workflow files. An unreachable daemon or security options without `rootless` abort the conversion before Nextflow starts. Other deployments retain the existing runtime-presence checks.

Combination preserves successful per-sample outputs when expression modalities, organisms, declared reference builds, or feature namespaces are incompatible. The result is marked partial, no combined H5AD is written, and the CLI returns status `1`.

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
JSONDataOutputOrchestrator.export_manifest(source, outdir, output_format)
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

The deterministic suite was refreshed on 2026-08-02 and reported
`464 passed, 3 skipped`. The public wire contract is Atlas document schema 1.0
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

`GEOParser.parse()` returns `list[dict]`, with one self-contained package per top-level MINiML `Series`.

```python
[
    {
        "version": str | None,
        "database": list[dict],
        "organization": list[dict],
        "contributor": list[dict],
        "platform": list[dict],
        "sample": list[dict],
        "series": dict,
    }
]
```

`AEParser.parse()` uses the same core package vocabulary but identifies its source dialect with `version = "magetabv1.1"` and the BioStudies MAGE-TAB specification URL in `schema_location`. Its `series.iid` is the explicit ArrayExpress accession, then an ArrayExpress-form investigation accession, then another ArrayExpress-classified accession, with the investigation accession as fallback. GEO secondary accessions remain in `series.accession` and do not displace an available ArrayExpress IID.

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

`GEOParser.parse_related_series(miniml, remove_empty=False, strict=True)` uses the same traversal logic but returns only related packages, excluding the input packages. When `strict=True`, fetch or parse failures raise. When `strict=False`, failed related accessions are skipped.

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

- Normal `geo2ae.convert()` enrichment happens after `GEOParser.parse()` and before `AEConstructor.miniml2magetab()`.
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
- `json2obs` accepts the same metadata and raw/processed data inputs and processed-checkpoint controls as `json2h5ad`, writes a required combined `obs.csv` with an explicit `cell_id` column, and can add typed `var.csv` and `uns.json` sidecars.

`main(argv=None) -> int`

- Parses arguments, configures `meta_standards_converter` logging to stdout and optional file output, creates the command's converter, and converts each accession or JSON file in order.
- Default logging emits `WARNING+`; `-v` emits `INFO+`, `-vv` emits `DEBUG+`, and `--quiet` emits `ERROR+`.
- On conversion failure, logs the exception traceback, marks the run failed, and continues.
- Success/progress messages are logged rather than printed; normal success output appears with `-v`.
- Returns `1` if any accession failed, otherwise `0`.

<a id="converter"></a>
### `converters/geo2ae.py`, `converters/geo2json.py`, `converters/json2ae.py`, `converters/ae2json.py`, `converters/json2h5ad.py`, and `converters/json2tabular.py`

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

- Accepts injectable `SourcePlanner`, `NFCoreRunner`, `AssetDownloader`, and
  ordered `AnnDataMetadataProjector` collaborators.
- `convert(..., allow_invalid=False) -> ConversionResult | BatchConversionResult` accepts ordinary
  parsed MINiML JSON or a canonical Atlas v1 document. It returns the
  single-group result directly and aggregates multiple groups.
- `convert_source(json_path, out=None, allow_invalid=False, **options) -> BatchConversionResult`
  always aggregates groups. Per-group exceptions populate `failures` and do
  not discard successful conversions.
- `ConversionResult` exposes `combined_h5ad`, `sample_h5ads`, retained pipeline files, pipeline commands, warnings/errors/failures, `primary_h5ad`, and `partial`.
- `AssetManifest` loads CSV/TSV mappings or `ACCESSION=PATH` CLI specifications. Manifest entries outrank CLI entries, which outrank discovered JSON assets.
- `AssetManifest.load(path: str) -> list[Asset]` reads CSV/TSV, requires
  `scope_id`/`path`, groups raw members, and raises `ValueError` for blank or
  unsupported records. `parse_spec(spec: str) -> Asset` parses the compact
  CLI form and raises on malformed specifications; neither writes files.
- `AssetDownloader.localize(value: str, md5: str | None = None) -> str`
  returns local paths unchanged or streams HTTP(S)/FTP into its cache, verifies
  an optional digest, and raises on transport/checksum failure; downloading is
  its filesystem/network side effect.
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
before MINiML attachment/writing; combined projectors run after
`anndata.concat()` and before combined MINiML attachment/writing. Scalar
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
### `enrichers/miniml_enricher.py`

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
### `geo_handlers/geo_webfetcher.py`

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
- Extracts `{gse}_family.xml` from the tarball and returns it as UTF-8 text.

<a id="ae-web-fetcher"></a>
### `ae_handlers/ae_webfetcher.py`

`TextResource(name, text, origin)` and `MAGETabInput(idf, sdrfs, source, source_kind)` are immutable transport records used between resolution and parsing.

`class AEWebFetcher`

- `__init__(requester=None, request_settings=None)` defaults to `RateLimitedRequester(service="biostudies")`.
- `resolve(source, sdrf_sources=None) -> MAGETabInput` dispatches existing paths, HTTP(S) URLs, and accession tokens. Missing path-like inputs raise `FileNotFoundError` instead of becoming accession lookups.
- Local and HTTP IDFs use explicit SDRF overrides when supplied; otherwise every `SDRF File` value is resolved relative to the IDF.
- Accession lookup requires exactly one discovered IDF and at least one SDRF. Explicit overrides are rejected for accession sources.
- Remote metadata is fetched as text and never persisted by the fetcher.

<a id="ae-parser"></a>
### `ae_handlers/ae_parser.py`

`class AEParser`

- `parse(source: MAGETabInput) -> dict` parses one IDF plus all SDRFs into the existing MINiML-compatible package shape.
- Root metadata uses the normalized MAGE-TAB format/version and specification URL. `series.iid` prefers `Comment[ArrayExpressAccession]`, then an ArrayExpress-form investigation/classified accession, then the investigation accession; `mage_tab.version` retains the exact IDF version.
- IDF rows are normalized by case and whitespace. Repeated row values remain ordered and feed investigation, accessions, design/factor, status, publication, contributor, database, and protocol records.
- SDRF headers map source/sample identities, characteristics, factors, protocol refs, platforms, technology, SRA/ENA runs, FASTQ metadata, and array raw/derived files. Repeated sample rows merge without duplicating list values.
- Conflicting scalar values keep the first value and append a warning. Unknown IDF rows and SDRF columns are preserved verbatim under `package["mage_tab"]` and also generate warnings.
- `build_model(...)` records complete protocol columns, QC/replicate/normalization declarations, every SDRF assay path, ordered nodes and protocol references, comments/files, and per-value unit/ontology companions under `package["mage_tab"]["model"]`.
- Malformed non-rectangular SDRF rows fail with a filename and column-count error.

<a id="ae-roundtrip"></a>
### `ae_handlers/ae_roundtrip.py`

- `semantic_sha256(package)` hashes every top-level public field except `mage_tab` using stable JSON encoding.
- `build_roundtrip(...)` stores schema version 1, core and typed-model fingerprints, complete parsed IDF rows, and named SDRF row tables under `mage_tab.roundtrip`.
- `unchanged_magetab(package)` returns the preserved source payload when the fingerprint still matches and exactly one SDRF is present.
- `restore_extensions(package, magetab)` runs after normal rendering for edited packages. Generated JSON fields win; unsupported IDF rows and SDRF columns are restored by row count or source/sample/assay identity when safe, otherwise a warning is logged.
- When a typed model is combined with a changed mapped core, `overlay_core()` performs a keyed union rather than an intersection-only replacement: missing allowlisted IDF rows and missing non-structural SDRF columns are inserted at core-order anchors. Model assay rows, material/assay node columns, and `Protocol REF` columns are never synthesized from the core projection.
- Packages without the optional sidecar follow the ordinary GEO/JSON rendering path unchanged. Multiple source SDRFs use semantic consolidation rather than the one-SDRF exact fast path.
- The fixed GEO/MINiML-compatible core still has no generic entities for arbitrary protocol graphs, assay/hybridization/scan identities, performers, protocol hardware/software/parameters, QC/replicate declarations, per-value units/ontology annotations, or custom MAGE-TAB fields. These are editable through `mage_tab.model` and backed by the raw sidecar; removing `mage_tab` is intentionally not lossless.

<a id="typed-mage-tab-model"></a>
### `ae_handlers/ae_model.py`

- `MAGETabModelError` is the public validation failure and
  `validate_model(model)` enforces schema version 1 collections, unique SDRF
  and assay identities, references, step shapes, and scalar harmonization
  annotations.
- `build_model(idf_rows, sdrfs)` creates the version-1 `mage_tab.model` extension without modifying the fixed MINiML projection.
- `protocols` contains one position-stable record per IDF protocol, including name, arbitrary type, ontology, description, hardware, software, parameters, contact, and performer.
- `declarations` independently stores aligned quality-control, replicate, and normalization terms with source/accession annotations.
- `assay_paths` contains one record per original SDRF data row. Ordered steps distinguish material/assay nodes, protocol references, annotated characteristics/factors/parameters, comments, files, and generic fields. This preserves array assay multiplicity and many-to-one sample relationships.
- Attribute steps keep `Unit`, `Term Source REF`, and `Term Accession Number` as independent fields; barcode/read geometry remains independent comment steps rather than being folded into protocol prose.
- `render_model(model)` regenerates one SDRF directly or consolidates multiple SDRFs by header plus occurrence. `overlay_core(model_rows, core_rows)` unions eligible core fields into that rendering while retaining model-only protocols, identities, annotations, rows, and structural graph columns.
- IDF matching uses normalized row labels and inserts only rows in the mapped allowlist. SDRF matching uses `(normalized header, occurrence)` keys, so repeated characteristics remain position-stable. Missing core columns are inserted relative to the nearest core-order neighbor; independent curator fields such as `Characteristics[hz_cell_type]`, `Characteristics[hz_cell_type_id]`, and `Characteristics[hz_cell_type_onto]` remain separate rather than being reinterpreted as native ontology companions.
- SDRF values align through the available `Sample Name`, `Source Name`, and `Comment[ENA_RUN]` identities. A core value replaces or populates a model cell only when all matching core rows agree on exactly one value. An unmatched model row keeps its existing value; an ambiguous newly inserted cell remains blank. Core-only rows are not added or broadcast as new assay paths.
- Editing `mage_tab.model` invalidates raw-table reuse through `model_sha256`. Old packages without a typed model or model hash continue through the existing sidecar or ordinary renderer.

<a id="geo-parser"></a>
### `geo_handlers/geo_parser.py`

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

`parse_related_series(miniml, remove_empty=False, strict=True) -> list[dict]`

- Parses the input MINiML, seeds a queue from related-series relations, and returns only fetched related packages.
- Deduplicates GSE accessions.
- Raises on fetch/parse failures in strict mode; skips failures in non-strict mode.
- Applies empty cleanup to related packages when requested.

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
### `ae_handlers/ae_idf_handlers.py`

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
### `ae_handlers/ae_constructor.py`

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

- `_detect_ae_technology()` chooses `bulk_sequencing`, `plate_single_cell_sequencing`, `droplet_single_cell_sequencing`, `tenx_v2_droplet_single_cell_sequencing`, `tenx_v3_droplet_single_cell_sequencing`, `spatial_sequencing`, `array`, or `generic`.
- `_has_array_files()` detects array-like files from platform/sample/series supplementary data and raw data.
- `_normalize_magetab_rows()` accepts row lists, comma-delimited legacy strings, and legacy `"SDRF file", sdrf` pairs.
- Legacy `_strip_quotes()` remains private but is no longer used by construction or file writing.
- `_sdrf_row_index()` finds the SDRF row case-insensitively.
- `_magetab_accession()` searches ArrayExpress, investigation, then secondary accession rows.
- `_safe_filename_token()` removes path separators from accession-derived filenames.
- `_is_table()` validates row-table shape.
- `_write_tsv()` writes `None` as blank cells.

<a id="sdrf-handlers"></a>
### `ae_handlers/ae_sdrf_handlers.py`

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
- `_SingleCellSequencingSDRFHandler` adds library construction, cDNA read size, technical replicate group, and study text helpers.
- `_DropletSingleCellSequencingSDRFHandler` adds 10x/droplet read geometry and isolation comments.
- `_TenXV2DropletSingleCellSequencingSDRFHandler` and `_TenXV3DropletSingleCellSequencingSDRFHandler` inherit the droplet path and emit fixed 10x chemistry library attributes such as cDNA read, cDNA read offset/size, barcode read/offset/size, end bias, input molecule, library construction, primer, strand, single-cell isolation, spike-in, and UMI geometry. v2 emits `Comment[library construction] = 10xV2`, cDNA read size `98`, and UMI barcode size `10`; v3 emits `10xV3`, cDNA read size `91`, and UMI barcode size `12`.
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
### `harmonizers/geo2ols.py`

`class GEO2OLS`

- Ensures an `ontologies` dict exists.
- Registers EFO and OBI term source metadata.

`geoprotocols2efo(protocol_type: str) -> list`

- Maps known MAGE-TAB/GEO protocol labels to ontology term, source ref, and accession.
- Maps `Library-Construction-Protocol` to `nucleic acid library construction protocol`, EFO, and `EFO_0004184`; registered EFO term-source metadata uses release `3.90.0`.
- Raises `ValueError` for blank protocol type.
- Returns `[protocol_type, None, None]` for unknown non-blank protocol labels, allowing custom protocol labels to survive in IDF output.

### `harmonizers/pubmed2ols.py`

`class Pubmed2OLS`

- Ensures an `ontologies` dict exists.
- Registers EFO and MeSH term source metadata.

`pubstatus2efo(pub_status: str) -> list`

- Returns `[None, None, None]` for blank status.
- Splits composite statuses on `+` and maps the first token.
- Maps common PubMed statuses such as `ppublish`, `epublish`, `pubmed`, `medline`, and `retracted`.
- Returns `[original_status_label, None, None]` for unknown non-blank statuses.

### `harmonizers/harmonizers.py`

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

`class RequestSettings`

- Stores request behavior: `timeout`, `request_delay`, `max_in_flight`,
  `max_retries`, retry HTTP statuses, exponential backoff base, and maximum
  backoff. Invalid time, delay, concurrency, or retry values fail at construction.
- Defaults retry HTTP statuses to `{429, 500, 502, 503, 504}`.

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
- Maintains shared per-host limiter state, so separate fetcher instances and
  different service labels targeting the same host respect the most conservative
  registered delay and concurrency ceiling. Different hosts do not block one another.
- Retries transient HTTP statuses. Numeric `Retry-After` headers control retry sleep; otherwise fallback delay is `min(0.5 * (2 ** attempt), 8.0)`.
- Retries `ConnectionError`, `Timeout`, and `ChunkedEncodingError` with the same deterministic exponential schedule and exact configured attempt count.
- Raises the exhausted retry response through `response.raise_for_status()`.

<a id="pubmed-fetcher"></a>
### `pubmed_handlers/pubmed_webfetcher.py`

`class PubmedWebFetcher`

`__init__(requester=None, request_settings=None)`

- Defaults to `RateLimitedRequester(service="ncbi_eutils")`.
- Accepts a custom requester or NCBI E-utilities request settings.

`fetch_pubmed_summary(pubmed_id: str) -> ET.Element`

- Calls NCBI PubMed ESummary for one PubMed ID through the `ncbi_eutils` requester.
- Raises for HTTP errors and returns the parsed XML root.

`pubmed_summary(pubmed_id: str) -> tuple`

- Parses DOI, author list, title, and PubMed publication status from ESummary XML.
- Maps publication status through `Harmonizer().pubstatus2efo()`.
- Returns the existing IDF tuple shape: DOI, author string, title, mapped status, term source ref, and term accession.

<a id="insdc-fetcher"></a>
### `insdc_handlers/insdc_webfetcher.py`

`class INSDCWebfetcher`

`__init__(ncbi_requester=None, ena_requester=None, ncbi_request_settings=None, ena_request_settings=None)`

- Defaults to `RateLimitedRequester(service="ncbi_eutils")` for NCBI SRA EFetch.
- Defaults to `RateLimitedRequester(service="ena_portal")` for ENA Portal file reports.
- Accepts custom requesters or per-service request settings.

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

<a id="metastore"></a>
### `meta_store/meta_store.py`

`class MetaStore`

`validate_investigation_metadata(investigation_metadata: dict) -> bool`

- Calls `_validate_investigation_metadata_structure()` and asserts it is truthy.
- Returns `True` only if validation passes.

`_validate_investigation_metadata_structure(investigation_metadata: dict) -> bool`

- Placeholder with `pass`.
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

The shared resolver validates schema version `1.0`, fixed MSC destinations or
`characteristics.<normalized_tag>`, and ordered canonical source fields. The
first populated `hz_<source>` wins. Scalar, numbered, container, and
characteristic representations carry value, ID, ontology, and hierarchy depth.
Invalid profiles warn and fall back to the complete raw view. Resolution uses a
deep copy: destination fields change only in the conversion view and every
`hz_*` field remains.

H5AD/obs publish canonical schema-3 columns, harmonization provenance columns,
and `uns["msc_harmonization"]`; `uns["msc_miniml"]` retains the untouched source.
TSV/CSV publish the same canonical and provenance view. MAGE-TAB replaces its
semantic destination, emits standard ontology companions where available, and
retains separate `Characteristics[hz_*]` columns. ECTO exposure and PCL
provisional-state annotations therefore pass through json2ae/ae2json,
json2tsv, json2h5ad, and json2obs without format-specific ontology code.

Importable implementation symbols are
`meta_standards_converter.converters.harmonization_overrides.HarmonizationSelection`,
`meta_standards_converter.converters.harmonization_overrides.HarmonizationResolution`,
`meta_standards_converter.converters.harmonization_overrides.resolve_harmonization_overrides`,
and
`meta_standards_converter.converters.harmonization_overrides.validate_harmonization_overrides`.

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

- `tests/test_geo_parser.py`: parser package scoping, cardinality, namespace handling, empty cleanup, related-series traversal, and fixture-backed parsing with `tests/GSE328265_family.xml`.
- `tests/test_geo2ae.py`: converter orchestration, related-series forwarding, enrichment, stage logging, and `remove_empty` forwarding.
- `tests/test_geo2json.py`: JSON converter orchestration, optional enrichment, JSON file writing, and stage logging.
- `tests/test_json2ae.py`: object/list loading, validation, default and skipped enrichment, MAGE-TAB writing, safe logging, independent fixture expectations, and extension restoration.
- `tests/test_ae2json.py`: IDF/SDRF mapping, typed protocol/declaration/assay-path capture, model edit authority, assay multiplicity, units/ontology, sidecar/fingerprint creation, unchanged lossless reuse, edited-core precedence, keyed IDF/SDRF overlay union, occurrence-aware duplicate headers, harmonized `hz_*` columns, ambiguity-safe row alignment, multiple SDRFs, conflicts, unmapped restoration, and output writing.
- `tests/test_ae_webfetcher.py`: local and HTTP relative resolution, explicit SDRF overrides, paginated BioStudies discovery/download calls, in-memory remote content, and invalid source metadata.
- `tests/test_json2h5ad.py`: asset precedence/manifests/downloads, canonical H5AD schema 3 metadata, normalized multivalue rows, smart observation IDs, opaque source-column preservation, real dictionary reference scoping, artifact-relative provenance, MINiML enrichment and publication filtering, count/TPM matrices, sparse combination, canonical/generic study splitting, partial results, and raw-output reintegration.
- `tests/test_atlas_v1_reader.py`: producer-owned golden fixture consumption, harmonized-state adaptation, structural validation, v1 cutover failure, and no-ThematicAtlases dependency proof.
- `tests/test_json_source.py`: native MINiML and Atlas v1 grouping, harmonized-status filtering, source diagnostics, and duplicate conflict handling.
- `tests/test_json2tabular.py`: neutral default columns, direct Atlas aggregation, replacement projectors, collisions, and validation behavior.
- `tests/test_metadata_projector.py`: generic sample/combined projector
  lifecycle, scalar broadcasting, axis-length validation, collision rejection,
  warning/error propagation, fail-closed output, invalid-output opt-in, and
  bundle rollback fault injection.
- `tests/test_external_guard.py`: fail-closed network/process guard self-tests and bounded fake-process opt-in.
- `tests/test_public_e2e.py`: offline public GEO/JSON/MAGE-TAB round trips and processed-H5AD bundle conversion without nf-core.
- `tests/test_h5ad_pipeline.py`: reference/annotation combinations, GFF3 conversion and reuse, FASTQ samplesheets, mixed modality grouping, pinned commands, warning extraction, output discovery, and workflow failure logs.
- `tests/test_h5ad_pipeline.py`: rootless enforcement also covers accepted, rootful, and unreachable Docker daemons.
- `tests/test_docker_artifacts.py`: pinned runtime tooling, rootless-only Compose mounts, hardening, and provisioning/runner script syntax.
- `tests/test_cli_geo2ae.py`: CLI defaults, multiple accession order, aliases, keep-empty behavior, out directory forwarding, logging controls, file logging, and failure continuation.
- `tests/test_cli_geo2json.py`: JSON CLI defaults, enrichment toggle, aliases, logging controls, file logging, and failure continuation.
- `tests/test_cli_json2ae.py`: JSON-to-MAGE-TAB CLI defaults, enrichment toggle, multiple input ordering, output forwarding, logging, and failure continuation.
- `tests/test_cli_ae2json.py`: MAGE-TAB-to-JSON CLI defaults, repeated SDRF overrides, source validation, multiple input ordering, logging, and failure continuation.
- `tests/test_cli_json2h5ad.py`: H5AD CLI defaults, workflow/reference/asset flags, partial status, multiple input order, logging, and failure continuation.
- `tests/test_cli_json2tabular.py`: TSV/CSV input order and partial exit status.
- `tests/test_project_scripts.py`: console script registration.
- `tests/test_docs_index.py`: stable documentation anchors, required README Guide structure including configuration, complete Mermaid platform-handler hierarchy coverage, interface-specific quickstart links, live-parser coverage for every documented CLI argument and alias, console-script mentions, docs links, and author-header policy.
- `tests/test_ae_constructor.py`: IDF rows, merged and source-aligned secondary accessions, protocol registry behavior, AE constructor sequencing, SDRF row insertion, file normalization, and protocol ref consistency.
- `tests/test_ae_sdrf_handlers.py`: SDRF graph rendering, source/comment/characteristic behavior, file classification, sequencing/array/single-cell/spatial handlers, SRA precedence warnings, and disabled greedy fallback comments.
- `tests/test_miniml_enricher.py`: additive PubMed/SRA enrichment fields, deduplication, and fetch error tolerance.
- `tests/test_request_helper.py`: timeout forwarding, shared service delays, retry statuses, `Retry-After`, exponential backoff, and exhausted retry errors.
- `tests/test_geo_webfetcher.py`: GEO URL handling, requester delegation, and MINiML tarball extraction.
- `tests/test_insdc_webfetcher.py`: SRA accession extraction, NCBI/ENA requester delegation, parsed SRA run records, and ENA fallback behavior.
- `tests/test_pubmed_webfetcher.py`: PubMed ESummary requester delegation, parsing, publication status mapping, and IDF constructor delegation.

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

`JSONDataOutputOrchestrator` publishes related tabular and AnnData metadata
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
`meta_standards_converter.converters.json2h5ad.DatasetBundleRecoveryError`,
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

`ae_handlers/ae_common.py` owns `ProtocolRegistry`, normalized file-extension
classification, array-file detection, and platform technology selection. Both
`AEConstructor` and `SDRFConstructor` import these contracts in one direction.
This removes their prior mutual import and method-local constructor imports
without changing handler keys, protocol references, technology decisions, or
the public `ae_constructor.ProtocolRegistry` import path.
Importable symbols are
`meta_standards_converter.ae_handlers.ae_common.ProtocolRegistry`,
`meta_standards_converter.ae_handlers.ae_common.detect_ae_technology`,
`meta_standards_converter.ae_handlers.ae_common.has_array_files`, and
`meta_standards_converter.ae_handlers.ae_common.normalized_extension`.

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
