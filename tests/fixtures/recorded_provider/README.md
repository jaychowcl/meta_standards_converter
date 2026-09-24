<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->


# Recorded provider audit fixtures

These fixtures are reviewed reductions of source-faithful scientific evidence,
not a replacement for the full 18-study provider audit.

* `success.manifest.json` exercises the deposited E-MTAB-16847 IDF/SDRF through
  JSON, MAGE-TAB and TSV. Expected factor cells come directly from deposited
  SDRF; sample/organism expectations come from retained, previously reviewed
  output at `1cec34cb3ee31fada67fbf4d2474ec1e2a7dfefe`.
* `conflict.manifest.json` preserves a recorded E-MTAB-16253 organism conflict.
  Every destination must remain execution-succeeded, validation-valid and
  completeness-partial, with the exact expected cause. This is parser/merge
  replay, not SRA HTTP discovery or a fresh successful E-MTAB-16253 retrieval.
* `row-evidence.manifest.json` hashes five reduced packages and complete
  projected row ledgers. Each file records its original immutable evidence
  path/hash. Reduced packages retain one source sample and run, that run's
  authored/native paths, sample-only supplementary paths and scoped conflicts.
  Other samples/runs, redundant raw-data collections and raw provider extension
  envelopes are omitted explicitly. Scientific values are not rewritten.
* `http403/manifest.json` retains the original E-MTAB-16253 failed SDRF response
  and its three source catalogue/IDF responses. Four explicit offline reuses of
  that one 403 body exercise the default initial attempt plus three retries;
  this is not evidence of four separate captured live requests. The failure
  cannot be replaced by a later 200 response.
* `http403/budget.manifest.json` is explicitly a synthetic budget fault at a
  recorded request identity. Its purpose is to verify expected exhaustion and
  harness-error separation; it is not a captured provider failure.

Row ledgers preserve original column index/header and exact cell strings,
including repeated headers and authored rows. Counts are 16, 1558, 104, 51 and
24 for E-MTAB-16253, E-MTAB-14566, E-MTAB-14416, E-GEOD-52564 and E-GEOD-60450.
The 17/12 GEO supplementary branches use `Comment[SAMPLE_FILE_*]` cells.
These source-hashed ledgers support accounting; they do not themselves prove
that a fresh complete provider replay has passed.

Run the self-contained gate with:

```sh
.venv/bin/python scripts/audit_recorded_provider.py tests/fixtures/recorded_provider/success.manifest.json
.venv/bin/python scripts/audit_recorded_provider.py tests/fixtures/recorded_provider/conflict.manifest.json
```

The runner never regenerates expected files. Manifest version 1.0 selects exact
prepared request identities, ordered response/exception attempts and repeats.
Response paths and input/golden paths have SHA-256 hashes. No folder precedence
or live fallback exists. A swallowed miss remains a harness error and exits 2.
External manifests require an explicit `--evidence-root` and report missing
evidence as an unmet prerequisite. Full raw HTTP and expanded GeoMx packages
remain in immutable development evidence, outside this tracked CI corpus.
Every scenario declares route intent, nonempty formats, one logical outcome,
and exact package/artifact counts. Production-standard routes retain default
enrichment services even when the default preset argument is omitted. Only
linked-only routes may explicitly select local-AE and no-publication services.
Optional audit artifact capture uses fresh files and retains failed-call payloads;
the first exact-golden mismatch includes a bounded JSON-pointer location.
