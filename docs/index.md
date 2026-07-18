<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# codebase.md Index

Read this file first, then retrieve only the relevant anchored sections from `docs/codebase.md`.
Use the `anchor` values below as stable header references. The codebase guide includes method/class responsibilities, pseudocode, internal calls, and external service calls.

## Agent Retrieval Template

```bash
anchor='<anchor-id>'
rg -n "<a id=\"${anchor}\"></a>|^#{2,6} " docs/codebase.md

rg -n "^#{1,6} " docs/codebase.md docs/index.md
```

## Main Sections

- id: project-purpose-and-layout
  title: Project Purpose And Layout
  anchor: project-purpose-and-layout
  keywords: purpose, layout, package tree, tests

- id: runtime-behavior
  title: Runtime Behavior
  anchor: runtime-behavior
  keywords: dependencies, console scripts, network, outputs, logging, DEBUG, INFO, request duration, parse stats, enrichment stats, safe telemetry

- id: end-to-end-geo2ae-flow
  title: End-To-End geo2ae Flow
  anchor: end-to-end-geo2ae-flow
  keywords: geo2ae, geo2json, json2h5ad, CLI, fetch, parse, MAGE-TAB, pseudocode, external APIs

- id: parsed-miniml-data-shape
  title: Parsed MINiML Data Shape
  anchor: parsed-miniml-data-shape
  keywords: MINiML, parser, JSON shape, enrichment fields

- id: workflow-details
  title: Workflow Details
  anchor: workflow-details
  keywords: workflows, related series, IDF, SDRF, enrichment

- id: json2h5ad-flow
  title: End-To-End json2h5ad Flow
  anchor: json2h5ad-flow
  keywords: H5AD, h5ad.gz, gzip, AnnData, matrix, FASTQ, nf-core, scrnaseq, rnaseq, assets, Nextflow

- id: rootless-json2h5ad-runtime
  title: Rootless json2h5ad Runtime
  anchor: rootless-json2h5ad-runtime
  keywords: Docker Compose, rootless, nfcore-runner, socket, security, provisioning, containers, ACL, permissions, 0660, NXF_OPTS, nextflow-tmp, noexec

- id: public-api-and-callable-reference
  title: Public API And Callable Reference
  anchor: public-api-and-callable-reference
  keywords: API, classes, methods, constructors, helpers, internal calls

- id: maintenance-notes
  title: Maintenance Notes
  anchor: maintenance-notes
  keywords: caveats, maintainers, generated output

- id: test-plan
  title: Test Plan
  anchor: test-plan
  keywords: tests, unittest, acceptance, coverage

## Workflow Sections

- id: geo-parse-flow
  title: GEO Parse Flow
  anchor: geo-parse-flow
  keywords: GEOParser, XML, package scoping, references

- id: related-series-flow
  title: Related-Series Flow
  anchor: related-series-flow
  keywords: related, superseries, subseries, queue, dedupe

- id: idf-and-mage-tab-construction-flow
  title: IDF And MAGE-TAB Construction Flow
  anchor: idf-and-mage-tab-construction-flow
  keywords: IDF, AEConstructor, ProtocolRegistry, composition

- id: sdrf-graph-and-rendering-flow
  title: SDRF Graph And Rendering Flow
  anchor: sdrf-graph-and-rendering-flow
  keywords: SDRF, graph, paths, columns, render

- id: technology-handler-selection
  title: Technology Handler Selection
  anchor: technology-handler-selection
  keywords: technology, handlers, sequencing, array, generic, 10x, v2, v3

- id: sequencing-sdrf-flow
  title: Sequencing SDRF Flow
  anchor: sequencing-sdrf-flow
  keywords: sequencing, SRA, ENA, FASTQ, derived files disabled, LIBRARY_SOURCE

- id: array-sdrf-flow
  title: Array SDRF Flow
  anchor: array-sdrf-flow
  keywords: array, hybridization, raw files, derived files disabled

- id: base-sdrf-behavior
  title: Base SDRF Behavior
  anchor: base-sdrf-behavior
  keywords: base handler, factors, protocols, sample ordering, characteristics, organism, MINiML value, required blank columns

- id: sra-pubmed-and-ontology-enrichment
  title: SRA, PubMed, And Ontology Enrichment
  anchor: sra-pubmed-and-ontology-enrichment
  keywords: SRA, PubMed, ontology, harmonizer

## API Sections

- id: cli
  title: CLI
  anchor: cli
  keywords: geo2ae, geo2json, json2h5ad, flags, logging, main

- id: converter
  title: Converter
  anchor: converter
  keywords: geo2ae, geo2json, json2h5ad, convert, JSON writing, H5AD, AnnData, nf-core

- id: miniml-enricher
  title: MINiML enricher
  anchor: miniml-enricher
  keywords: enrichment, PubMed, SRA, sample metadata

- id: geo-web-fetcher
  title: GEO web fetcher
  anchor: geo-web-fetcher
  keywords: GEOWebFetcher, MINiML, GEO FTP, tarball

- id: geo-parser
  title: GEO parser
  anchor: geo-parser
  keywords: GEOParser, parse helpers, related series

- id: ae-idf-handlers
  title: AE IDF handlers
  anchor: ae-idf-handlers
  keywords: IDFConstructor, platform IDF, secondary accession, GEO, ENA, SRA, DRA, publications, dates, protocols

- id: ae-constructor
  title: AE constructor
  anchor: ae-constructor
  keywords: AEConstructor, ProtocolRegistry, technology detection, file writing

- id: sdrf-handlers
  title: SDRF handlers
  anchor: sdrf-handlers
  keywords: SDRF, handlers, graph, file classification

- id: harmonizers
  title: Harmonizers
  anchor: harmonizers
  keywords: Harmonizer, GEO2OLS, Pubmed2OLS, ontology

- id: json-helper
  title: JSON helper
  anchor: json-helper
  keywords: JSONHandler, dotted paths, flattening

- id: request-helper
  title: Request helper
  anchor: request-helper
  keywords: RateLimitedRequester, timeout, retries, backoff, requests.get, external APIs

- id: pubmed-fetcher
  title: PubMed fetcher
  anchor: pubmed-fetcher
  keywords: PubMed, ESummary, DOI, publication status

- id: insdc-fetcher
  title: INSDC fetcher
  anchor: insdc-fetcher
  keywords: SRA, INSDC, ENA, FASTQ, runs

- id: metastore
  title: MetaStore
  anchor: metastore
  keywords: MetaStore, validation, placeholder

## Parser Callables

- id: geoparser-class-and-parse-methods
  title: GEOParser class and parse methods
  anchor: geoparser-class-and-parse-methods
  keywords: repeated_children, parse, cleanup

- id: parser-reference-resolution
  title: Reference resolution
  anchor: parser-reference-resolution
  keywords: package, samples, platforms, contributors, databases

- id: parser-generic-xml-mapping
  title: Generic XML mapping
  anchor: parser-generic-xml-mapping
  keywords: _parse_element, snake_case, XML text

- id: parser-related-series-helpers
  title: Related-series helpers
  anchor: parser-related-series-helpers
  keywords: GSE extraction, relation matching

- id: parser-cleanup-and-helpers
  title: Cleanup and helpers
  anchor: parser-cleanup-and-helpers
  keywords: remove_empty, namespace, list normalization

## SDRF Callables

- id: sdrf-dataclasses
  title: SDRF dataclasses
  anchor: sdrf-dataclasses
  keywords: SDRFAttr, SDRFNode, SDRFEdge, SDRFPath, SDRFAudit

- id: sdrfconstructor
  title: SDRFConstructor
  anchor: sdrfconstructor
  keywords: SDRFConstructor, technology, SRA lookup

- id: sdrf-file-helpers
  title: File helpers
  anchor: sdrf-file-helpers
  keywords: classify_file, normalized_extension

- id: base-sdrf-handler
  title: Base SDRF handler
  anchor: base-sdrf-handler
  keywords: paths, columns, source, factors

- id: sequencing-handlers
  title: Sequencing handlers
  anchor: sequencing-handlers
  keywords: sequencing, bulk, single-cell, droplet, 10x, v2, v3, spatial

- id: array-and-generic-handlers
  title: Array and generic handlers
  anchor: array-and-generic-handlers
  keywords: array, generic, files, hybridization

- id: legacy-fallback-notes
  title: Legacy fallback notes
  anchor: legacy-fallback-notes
  keywords: fallback, disabled, comments
