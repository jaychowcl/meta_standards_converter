<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
<a id="provider-concurrency-benchmark"></a>
# Public-provider concurrency benchmark policy

This benchmark is an explicit, bounded live operation. It never invokes a
language model or embedding provider. It uses fixed public identifiers, one
warm-up, three measured sequential runs, and three measured two-worker runs.

The default ceilings intentionally leave margin below documented or observed
provider capacity:

- NCBI E-utilities: two request starts per second and two requests in flight.
  NCBI documents three requests per second without an API key and ten with one.
- ENA Portal and BioStudies: one request start per second and two requests in
  flight. Their public documentation does not publish a numeric allowance, so
  the client uses the conservative common EBI ceiling.

Two workers are recommended only when median wall time improves by at least
20%, response content and order are identical, every response succeeds without
retryable status, request p95 is no more than 1.5 times the sequential p95, and
the measured start spacing stays within the configured ceiling. A failed
criterion leaves sequential execution as the default.

Run from the repository root:

```bash
PYTHONPATH=src python benchmarks/provider_http_concurrency.py \
  --provider all --runs 3 --out benchmarks/results/provider-http-YYYY-MM-DD.json
```

Provider references:

- [NCBI E-utilities usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
- [ENA Browser API](https://www.ebi.ac.uk/ena/browser/api/)
- [BioStudies API](https://www.ebi.ac.uk/biostudies/help#API)

## 2026-08-02 result

All requests succeeded, returned deterministic content in input order, and
respected start-rate ceilings. Median sequential versus two-worker wall times
were 1.297s versus 1.248s for NCBI (+3.8%), 2.118s versus 2.120s for ENA
(-0.1%), and 2.142s versus 2.195s for BioStudies (-2.5%). NCBI and ENA also
failed the tail-latency criterion. No provider met the 20% threshold, so no
production workflow enables provider concurrency from this benchmark. The raw
report is [`benchmarks/results/provider-http-2026-08-02.json`](../benchmarks/results/provider-http-2026-08-02.json).
