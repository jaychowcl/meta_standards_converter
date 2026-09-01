<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# ENA converter reference snapshot

This directory is a source-faithful, offline reference bundle for designing a future `ena2json` converter. It documents provider capabilities and examples; it does not implement retrieval or conversion. Provider bytes are unchanged. [manifest.json](manifest.json) records URLs, the 2026-09-01 retrieval timestamp, content types, available HTTP validators, SHA-256 digests, local paths, and fixture accessions.

## Contents

| Path | Contents |
|---|---|
| `schemas/ENA.*.xsd` | Project, checklist, assembly, taxonomy, sample-group, EMBL, root, and Webin schemas from the ENA INSDC SRA 1.5 directory |
| `openapi/` | Portal, Browser, Xref, and Taxonomy OpenAPI snapshots |
| `catalogues/results.tsv` | Portal result types and snapshot record counts/update dates |
| `catalogues/return-fields/` | Return-field contracts for nine converter-relevant result types |
| `catalogues/search-fields/` | Search-field contracts for the same result types |
| `catalogues/controlled-vocab/` | All 18 controlled vocabularies declared by the Portal API at snapshot time |
| `catalogues/checklists/` | Every checklist XML successfully returned by the documented Browser endpoint |
| `fixtures/` | Browser XML plus Portal search/file-report, Xref, and Taxonomy responses for two linked public datasets |

## Snapshotted services

| Service | Principal route | Information exposed |
|---|---|---|
| Portal API | `/search`, `/results`, `/returnFields`, `/searchFields`, `/controlledVocab` | Indexed metadata search, result/field discovery, and controlled values |
| Portal file report | `/filereport` | Fast access to run/analysis file URLs, hashes, and byte sizes |
| Browser API | `/xml/{accession}` and XML search routes | Full archive XML for studies, samples, experiments, runs, analyses, assemblies, and other records |
| Xref API | `/json/search?accession=...` | Cross-reference records where the Xref index has matches; an empty JSON list is valid |
| Taxonomy API | `/tax-id/{taxId}` | Scientific name, rank, lineage, division, genetic codes, alternative names, and flags |

Portal responses contain indexed columns only; arbitrary submitter metadata can require Browser XML. Treat `returnFields`/`searchFields` as live endpoint contracts, not a universal guarantee that a value is populated for every record.

## Schemas and dependency resolution

The ENA XSDs are retained exactly as published. `ENA.assembly.xsd`, `ENA.checklist.xsd`, `ENA.project.xsd`, and `ENA.sample_group.xsd` import `SRA.common.xsd`; `ENA.root.xsd` and `ENA.webin.xsd` include multiple SRA object schemas. Resolve those names against `../sra/schemas/` in an XML catalog or a temporary validation directory—do not edit provider `schemaLocation` values.

`ENA.root.xsd` also retains historical FTP include locations. `ENA.webin.xsd` includes EGA schemas that are deliberately outside this bundle: the approved scope excludes unrelated `EGA.*.xsd` files. Consequently, XML syntax is tested offline, while complete Webin-schema compilation requires those external optional dependencies.

## Checklist snapshot caveat

The Portal controlled vocabulary declared 47 checklist identifiers on 2026-09-01. The documented Browser route returned XML for 31 (`ERC000011` through `ERC000041`) and HTTP 404 for 16 identifiers. All 47 declarations remain preserved in `controlled-vocab/checklist.json`; no XML was fabricated for unavailable records. [checklist-availability.md](checklist-availability.md) records the exact gap. In particular, `ERC100001` is the BioSamples minimal checklist described by BioSamples as JSON-schema validation rather than an ENA Browser XML record.

## Fixtures

| Fixture | Chain and representative behavior |
|---|---|
| `SRX017289` | `SRP002056` → `SRS011830` → `SRX017289` → `SRR037073`; `PRJNA123835`, `SAMN00009557`, GEO/PubMed-linked RNA-Seq, single FASTQ |
| `SRX7812918` | `SRP250911` → `SRS6225446` → `SRX7812918` → `SRR11192680`; `PRJNA609050`, `SAMN14218700`, human-gut amplicon, paired layout and two semicolon-aligned FASTQs, no publication |

Each fixture has Browser XML at study, sample, experiment, and run level, minimal Portal searches for all four levels, a file report, an Xref response, and a taxonomy response. Both Xref fixtures are valid empty arrays; embedded XML links and Portal accessions still establish the hierarchy.

## Integrity policy

Do not reformat JSON, TSV, XML, XSD, or OpenAPI snapshots. Refresh provider files from their manifest URLs and update hashes and retrieval metadata together. Project-authored Markdown keeps the canonical author header. Offline verification lives in `tests/test_provider_reference_material.py`.
