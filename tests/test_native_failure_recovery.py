# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.sources.metadata import NativeMetadataSource
"""Native failures must remain retryable through the public converter facade."""
from meta_standards_converter.converters import PreparationPolicy
from types import SimpleNamespace
from unittest.mock import Mock
import pytest
import requests

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.ena2json import ENA2JSONConverter
from meta_standards_converter.helpers.request_helper import HostRequestCooldownDeferred
from meta_standards_converter.sources.ena import ENASource
from tests.test_native_archive_sources import HTTP


@pytest.mark.parametrize('error', [requests.HTTPError('server unavailable', response=Mock(status_code=500)),
    requests.exceptions.ChunkedEncodingError('truncated'),
    HostRequestCooldownDeferred('www.ebi.ac.uk/ena', 2000000000)])
def test_no_native_package_retains_transient_error_and_cooldown(error):
    def fail(*args): raise error
    service=NativeMetadataSource("ena", source=ENASource(http=HTTP(fail)))
    result=Converter(services={'ena_source': service}).convert('ERP1', in_type='ena_accession', out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.items[0].execution=='failed'
    assert result.items[0].error.retry_category.value=='retryable'
    if hasattr(error,'retry_at'):
        assert result.items[0].retry_at==error.retry_at


def test_permanent_not_found_is_not_whole_job_retryable():
    error=requests.HTTPError('not found',response=Mock(status_code=404))
    def fail(*args):raise error
    result=Converter(services={'ena_source': NativeMetadataSource('ena', source=ENASource(http=HTTP(fail)))}).convert('ERP1', in_type='ena_accession', out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.items[0].error.retry_category.value=='terminal'


def test_partial_native_package_retains_transient_failure_without_discarding_artifacts(tmp_path):
    from tests.test_ena_missing_record_recovery import missing_sample_records
    from meta_standards_converter.sources.archive_support import Resolution, StudySeed
    records=missing_sample_records()
    error=HostRequestCooldownDeferred('www.ebi.ac.uk/ena', 2000000000)
    records.failures.append(error)
    records.issues.append('sample XML unavailable')
    source=SimpleNamespace(resolve=lambda accession:Resolution([StudySeed('ERP1','PRJEB1')]),
                           fetch=lambda seed, **options:records)
    result=Converter(services={'ena_source': NativeMetadataSource('ena', source=source)}).convert('ERP1', in_type='ena_accession', out_type='json', outdir=tmp_path, preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    item=result.items[0]
    assert item.execution=='succeeded' and item.completeness=='partial'
    assert item.artifacts
    assert item.error.retry_category.value=='retryable'
    assert item.retry_at==error.retry_at
