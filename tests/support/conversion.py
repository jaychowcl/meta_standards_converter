# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Assertions shared by scientific integration tests using the public MSC API."""

def converted(converter, source, **request):
    batch = converter.convert(source, **request)
    assert len(batch.items) == 1, batch.to_dict()
    item = batch.items[0]
    assert item.execution == 'succeeded' and item.validation == 'valid', batch.to_dict()
    return item.payload
