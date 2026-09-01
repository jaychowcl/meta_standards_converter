<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# ENA checklist availability on 2026-09-01

The Portal endpoint `controlledVocab?field=checklist` declared 47 values. Requests to the documented Browser route `https://www.ebi.ac.uk/ena/browser/api/xml/<accession>` returned parseable XML for 31 values: every accession from `ERC000011` through `ERC000041`. The following declared values returned HTTP 404 and therefore have no vendored XML:

| Accession | Snapshot outcome | Interpretation |
|---|---|---|
| `ERC000004` | Browser XML HTTP 404 | Historical checklist identifier; ENA announced deprecation of an older sample-checklist XML in 2014 |
| `ERC000042` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000043` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000044` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000045` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000047` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000048` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000049` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000050` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000051` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000052` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000053` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000055` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000056` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC000058` | Browser XML HTTP 404 | Still declared by Portal; no provider XML available from the documented route |
| `ERC100001` | Browser XML HTTP 404 | BioSamples documents this as its minimal checklist and validates it as JSON Schema, not ENA Browser XML |

This is provider evidence, not a normalization decision. A future refresh may change the split. Converter code should use the declared checklist ID as provenance, tolerate unavailable definitions, and avoid treating “listed by Portal” as proof that Browser XML exists.

References: [ENA sample-programming documentation](https://ena-docs.readthedocs.io/en/latest/submit/samples/programmatic.html), [ENA 2014 checklist update](https://www.ebi.ac.uk/about/news/updates-from-data-resources/update-ena-sample-checklist/), and [BioSamples release notes](https://www.ebi.ac.uk/biosamples/docs/releasenotes).
