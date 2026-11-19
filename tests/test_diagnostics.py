# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Actionable, correlated diagnostics without exception payload disclosure."""
import errno
import json
import logging
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.unified.contracts import Diagnostic
from meta_standards_converter.conversion_logging import operation, labels
from meta_standards_converter.provider_issues import collect_provider_issues, record
from meta_standards_converter.sources.errors import MetadataSourceError
from tests.support.conversion import FailingSource


def test_diagnostic_additions_preserve_legacy_shape():
    assert Diagnostic('test', 'message').to_dict() == dict(
        code='test', message='message', stage='input', severity='error')
    diagnostic = Diagnostic('test', 'message', action='Retry later', provider='ena',
                            operation='source_loading', location='study.xml', correlation_id='cause-1')
    assert diagnostic.to_dict()['action'] == 'Retry later'
    assert diagnostic.to_dict()['correlation_id'] == 'cause-1'


@pytest.mark.parametrize('error,code', [
    (FileNotFoundError('private exception payload'), 'missing_file'),
    (PermissionError('private exception payload'), 'permission_denied'),
    (OSError(errno.ENOSPC, 'private exception payload'), 'disk_full'),
    (json.JSONDecodeError('private exception payload', '{private payload', 1), 'invalid_json'),
    (ET.ParseError('private exception payload'), 'invalid_xml'),
])
def test_expected_failures_have_actions_and_share_error_identity(error, code):
    item = Converter(handlers=[FailingSource(error)]).convert('GSE1', out_type='json').items[0]
    diagnostic = item.diagnostics[-1]
    assert diagnostic.code == code
    assert diagnostic.action
    assert diagnostic.correlation_id == item.error.correlation_id
    assert 'private' not in json.dumps(item.to_dict())


def test_provider_event_uses_failed_provider_not_active_source(caplog):
    from meta_standards_converter.converters.unified.contracts import ConversionItemResult
    item = ConversionItemResult('input_1', 'E-MTAB-16245')
    with caplog.at_level(logging.DEBUG), operation('test', heartbeat=False), labels(provider='ena'), collect_provider_issues(item):
        record(MetadataSourceError('source_inaccessible', 'E-MTAB-16245', http_status=403), 'biostudies')
    event = next(r.msc_event for r in caplog.records if hasattr(r, 'msc_event') and r.msc_event['event_type']=='provider_issue')
    assert event['attributes']['provider'] == 'biostudies'
    assert item.diagnostics[0].correlation_id == item.diagnostics[0].provider_issue.error.correlation_id


def test_unknown_failure_has_safe_debug_stack(caplog):
    with caplog.at_level(logging.DEBUG):
        item = Converter(handlers=[FailingSource(RuntimeError('private_payload_canary'))]).convert('GSE1', out_type='json').items[0]
    events = [r.msc_event for r in caplog.records if hasattr(r, 'msc_event') and r.msc_event['event_type']=='internal_failure']
    assert len(events) == 1
    assert events[0]['attributes']['correlation_id'] == item.error.correlation_id
    assert events[0]['attributes']['causes'][0]['frames']
    assert 'private_payload_canary' not in caplog.text
    assert 'private_payload_canary' not in json.dumps(item.to_dict())


def test_wrapping_preserves_an_existing_nonprovider_cause():
    from meta_standards_converter.converters.unified.outcomes import failure
    with operation('test', heartbeat=False):
        try:
            raise PermissionError('private_payload_canary')
        except PermissionError as original:
            first, _ = failure(original, stage='publication', provider=None, item_id='one')
            try:
                raise RuntimeError('private wrapper') from original
            except RuntimeError as wrapped:
                second, _ = failure(wrapped, stage='conversion', provider=None, item_id='one')
        assert first.correlation_id == second.correlation_id


def test_same_exception_in_distinct_operations_gets_distinct_debug_evidence(caplog):
    error = RuntimeError('private_payload_canary')
    with caplog.at_level(logging.DEBUG):
        results = [Converter(handlers=[FailingSource(error)]).convert('GSE1', out_type='json').items[0] for _ in range(2)]
    assert results[0].error.correlation_id != results[1].error.correlation_id
    events = [r.msc_event for r in caplog.records if hasattr(r, 'msc_event') and r.msc_event['event_type']=='internal_failure']
    assert len(events) == 2


def test_broken_external_logger_cannot_replace_original_exception():
    from meta_standards_converter.operational_events import OperationalEventEmitter
    class Broken:
        def log(self, *args, **kwargs):
            raise OSError('private sink')
    emitter = OperationalEventEmitter(repository='test', logger=Broken())
    with pytest.raises(ValueError, match='original'):
        with operation('test', heartbeat=False, emitter=emitter):
            raise ValueError('original')
