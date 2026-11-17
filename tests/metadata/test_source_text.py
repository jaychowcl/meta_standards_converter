# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Only corroborated same-field reversible text transformations are selected."""
from meta_standards_converter.converters import PreparationPolicy
from meta_standards_converter.converters import Converter
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.converters.test_audit_source_precedence import package
from tests.metadata.test_portable_conflicts import import_tables

CLEAN = 'Manufacturer’s instructions for Müller cells'
BROKEN = CLEAN.encode('utf-8').decode('cp1252')


def test_uncorroborated_source_text_is_retained_and_flagged_even_off():
    data = package(); data['sample'][0]['description'] = BROKEN
    result = Converter().convert(data, out_type='json', preparation=PreparationPolicy(enrichment='off'))
    assert result.items[0].payload[0]['sample'][0]['description'] == BROKEN
    assert result.status == 'partial' and result.items[0].validation == 'valid'
    assert any('encoding' in d.message and '/sample/SRS1/description' in d.message for d in result.items[0].diagnostics)


def test_corroborated_readable_statement_and_original_survive_archive_roundtrips():
    # Clean lower-ranked evidence can corroborate a reversible transformation.
    a, b = package('MAGE-TAB'), package()
    a['sample'][0]['description'] = BROKEN; b['sample'][0]['description'] = CLEAN
    merged, issues = merge_archive_metadata(MINiMLCodec().decode(a).package, MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    assert merged.to_mapping()['sample'][0]['description'] == CLEAN
    data = merged
    for _ in range(2):
        records = data.to_mapping()['extensions']['insdc']['records']
        normalizations = [r for r in records if r['kind'] == 'text_normalization']
        assert len(normalizations) == 1
        fact = normalizations[0]['metadata']
        assert {s['value'] for s in fact['statements']} == {CLEAN, BROKEN}
        assert fact['path'] == '/sample/SRS1/description'
        tables = Converter().convert(data, out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
        data, _ = import_tables(tables)


def test_curators_retains_primary_text_and_unrelated_fields_do_not_corroborate():
    for mode in ['standard', 'curators']:
        a, b = package('MAGE-TAB'), package()
        a['sample'][0]['description'] = BROKEN; b['sample'][0]['title'] = CLEAN
        if mode == 'curators': b['sample'][0]['description'] = CLEAN
        merged, _ = merge_archive_metadata(MINiMLCodec().decode(a).package, MINiMLCodec().decode(b).package, prefer=True, mode=mode)
        assert merged.to_mapping()['sample'][0]['description'] == BROKEN


def test_valid_unicode_is_not_changed_or_flagged():
    data = package(); data['sample'][0]['description'] = CLEAN + ' – déjà vu, São Paulo'
    result = Converter().convert(data, out_type='json', preparation=PreparationPolicy(enrichment='off'))
    assert result.status == 'complete'
    assert result.items[0].payload[0]['sample'][0]['description'] == data['sample'][0]['description']


def test_provider_bytes_override_incorrect_http_charset_and_keep_inherited_corruption():
    from meta_standards_converter.sources.archive_support import ArchiveHTTP
    from requests import Response
    class Requester:
        def __init__(self, value): self.value = value
        def get(self, *args, **kwargs):
            response = Response(); response.status_code = 200
            response._content = ('<?xml version="1.0" encoding="UTF-8"?><SAMPLE><DESCRIPTION>' + self.value + '</DESCRIPTION></SAMPLE>').encode('utf-8')
            response._content_consumed = True
            response.encoding = 'latin1'  # Provider header must not reinterpret XML bytes.
            return response
    for literal in (CLEAN, BROKEN):
        root = ArchiveHTTP('ncbi_eutils', requester=Requester(literal)).get('https://example.test/metadata')
        assert root.findtext('DESCRIPTION') == literal


def test_unknown_normalization_comment_is_retained_without_interpretation():
    tables = Converter().convert(package(), out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
    literal = '{"schema_version":"9.0","record":"do not interpret"}'
    tables.append(['Comment[MSCTextNormalization]', literal])
    data, warnings = import_tables(tables)
    assert any('MSCTextNormalization' in w for w in warnings)
    rows = Converter().convert(data, out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
    assert literal in [v for r in rows if r[0] == 'Comment[MSCTextNormalization]' for v in r[1:]]
