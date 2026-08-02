<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->

# Rootless nf-core Acceptance — 2026-07-31

This report preserves the portable outcome of the dedicated rootless raw-data
acceptance run. It is evidence for the process boundary documented in
[`codebase.md`](codebase.md#rootless-json2h5ad-runtime), not a command that the
deterministic test suite repeats.

| Workflow | Pinned version | Process result | Converter result |
| --- | --- | --- | --- |
| `nf-core/rnaseq` | 3.26.0 | return code 0 | non-partial H5AD |
| `nf-core/scrnaseq` | 4.2.0 | return code 0 | non-partial H5AD |

In compact release notation, the accepted workflows were rnaseq 3.26.0 and
scrnaseq 4.2.0.

The run used the repository's dedicated rootless runner boundary. No rootful
Docker socket was mounted. The ordinary test suite continues to replace
Nextflow and container execution with bounded fake processes; opt-in provider
tests remain separate from this operational acceptance evidence.
