<!--
Authors

Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
https://github.com/jaychowcl
https://saezlab.org
https://www.ebi.ac.uk/about/teams/functional-genomics/
-->
# Supplied-source fidelity expectations

Reviewed on 2026-11-25 against the original GEO XML and saved MINiML inputs.
The older `../reviewed` and `../../../chemistry/routing` expectations remain
unchanged as historical evidence. The current exact golden tests use these
additional expectations for the repaired supplied-identifier projection.

The only changes are additive evidence: NCBITaxon declarations, organism-local
9606 (human) or 10090 (mouse) qualifiers, and the supplied SRX32831930 experiment
on its explicitly owned assay. The latter is separate from its existing run
annotation. No files, rows, whitespace in source text, protocol definitions,
checksums or pre-existing scientific cells are changed. JSON type/order and full
IDF/SDRF equality remain tested; no comparison normalization hides the repair.

The seven original MSC-LIVE expected failures now pass as ordinary regressions.
URI-form NCBITaxon identifiers equivalent to a numeric taxid retain their authored
spelling and do not become a conflict or additional organism occurrence.
