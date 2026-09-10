<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# Scoped chemistry evidence and reviewed expectations

This corpus contains seven **source excerpts**, not seven complete deposited studies.
`manifest.json` records the URL, review/acquisition date, reduction and SHA-256
of each excerpt JSON and constructed XML envelope. Source text was transcribed
from the provider pages reviewed on 2026-09-10; no FASTQs or matrices were acquired.
An excerpt checksum seals the retained excerpt, not an unretained full webpage.

The XML files use synthetic GSE90000x/GSM90000x identities, platform, organism and
source information solely to exercise actual GEO parsing, enrichment and both
MAGE-TAB converters. These envelopes must never be interpreted as real biological
records. Only their protocol/processing excerpts come from the cited sources.
The independently retained real GSE328265/GSM9651991 fixture remains under
`studies/GSE328265` and has its own acquisition provenance.

## Expected scientific fields

- GSM6975354: 3-prime v2. Cell Ranger v3 in processing is not chemistry evidence.
- GSM8288494: both 3-prime v3.1 and 5-prime v2 preparation phrases. End bias and
  construction version are omitted. Family/version ambiguity is audited.
- PBMC 5-prime v2: explicit transcript R2=90; R1 barcode=16 and UMI=10; both index
  reads have 10 cycles. No offsets are inferred.
- ALL PBMC 5-prime v3: the 3-prime/5-prime fixation reference does not identify
  the library chemistry. Both indices have 10 reported cycles. R2=90 is reported
  without a read role in this excerpt, so it is **not** labelled cDNA read size.
  The shared GEX/VDJ description produces an ambiguous library-role diagnostic.
- PBMC 3-prime v3.1: explicit transcript R2=90; R1 barcode=16 and UMI=12; both index
  reads have 10 cycles. No offsets are inferred.
- GSM1520438: bacterial 5-prime end mapping is not a 10x assay. No 10x attributes.
- Flex: explicitly separate family/version. No conventional transcript geometry
  or 3-prime/5-prime end bias is inferred from its version.

Complete SDRFs were inspected by named column against these expectations,
including absent fields and repeated index1/index2 occurrences. Non-chemistry
columns describe only the constructed envelope. The runtime tests never accept
or regenerate expectations; the root checksum inventory seals the full corpus.
The source-backed checks in `tests/magetab/test_chemistry.py` independently assert
family, version alternatives and explicitly reported read/barcode fields.

## Round-trip boundaries

Both GEO2AE and JSON2AE must publish exactly the stored SDRFs. Real AE2JSON parsing
must retain the original published protocol descriptions and every extract-node
comment, including repeated index attributes, in order. Subsequent JSON2AE output
must retain the same ordered comments. Existing protocol registration flattens
whitespace for IDF, and the existing semantic renderer groups characteristics
before source comments; byte-identical regeneration is not claimed. Neither
behavior was changed by the chemistry correction.

Normal tests block network access and replace provider responses only. Unknown
kits, excluded wording, missing family/version, conflicting recipes, mixed
samples, forced handler selections, explicit offsets and library-name scope
also have constructed unit cases. This bounded grammar is not a general text
understanding system, and a passing corpus is not an estimated population-wide
accuracy rate.
