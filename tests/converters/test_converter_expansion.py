# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Only verified hierarchy edges expand; enrichment remains independent."""
from copy import deepcopy
from xml.etree import ElementTree as ET

import pytest

from meta_standards_converter.converters import Converter
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.insdc_support import tree
from tests.converters.test_json2tsv import package


def geo(accession, children=()):
    value = package(accession, 'GSM' + accession[3:])
    value['series']['iid'] = accession
    value['source']['format'] = 'GEO MINiML'
    value['series']['relation'] = [{'type': 'SubSeries', 'target': a} for a in children]
    return value


class Geo:
    def __init__(self, documents):
        self.documents, self.calls = documents, []
    def convert(self, accession, *, enrich, related_series, **kwargs):
        assert enrich is related_series is False
        self.calls.append(accession)
        value = self.documents[accession]
        if isinstance(value, Exception):
            raise value
        return [MINiMLCodec().decode(value).package]


@pytest.mark.parametrize('preset', ['off', 'curators', 'standard'])
def test_transitive_family_is_separate_deterministic_and_cycle_safe(preset):
    documents = {'GSE1': geo('GSE1', ['GSE3', 'GSE2']),
                 'GSE2': geo('GSE2', ['GSE1', 'GSE4']),
                 'GSE3': geo('GSE3', ['GSE4']), 'GSE4': geo('GSE4')}
    source = Geo(documents)
    class Enricher:
        def enrich_selected(self, value, **kwargs):
            return value
    result = Converter(services={'geo2json': source, 'enricher': Enricher()}).convert(
        'GSE1', out_type='json', enrichment=preset)
    assert result.status == 'complete', result.to_dict()
    assert source.calls == ['GSE1', 'GSE2', 'GSE3', 'GSE4']
    assert result.items[0].dataset_ids == ('GSE1', 'GSE2', 'GSE3', 'GSE4')
    assert len(result.items[0].payload) == 4
    assert all(len(p.samples) == 1 for p in result.items[0].payload)
    assert len([r for r in result.items[0].preparation if r['operation'] == 'study_expansion' and r['status'] == 'completed']) == 3


def test_saved_package_expansion_retains_input_and_partial_success():
    value = geo('GSE1', ['GSE2', 'GSE3'])
    before = deepcopy(value)
    source = Geo({'GSE2': geo('GSE2'), 'GSE3': OSError('unavailable')})
    result = Converter(services={'geo2json': source}).convert(value, out_type='json', enrichment='off')
    assert result.status == 'partial', result.to_dict()
    assert result.items[0].execution == 'succeeded'
    assert result.items[0].dataset_ids == ('GSE1', 'GSE2')
    assert value == before
    assert any(r['status'] == 'failed' for r in result.items[0].preparation)


def test_expansion_can_be_disabled_and_never_follows_free_text():
    value = geo('GSE1', ['GSE2'])
    source = Geo({})
    converter = Converter(services={'geo2json': source})
    assert converter.convert(value, out_type='json', enrichment='off', options={'expand_studies': False}).status == 'complete'
    value['series']['relation'] = [{'type': 'citation', 'target': 'GSE2', 'comment': 'SubSeries GSE2'}]
    value['series']['summary'] = ['superseries GSE3']
    assert converter.convert(value, out_type='json', enrichment='off').status == 'complete'
    assert not source.calls


def test_expansion_rejects_wrong_retrieved_identity():
    source = Geo({'GSE2': geo('GSE999')})
    result = Converter(services={'geo2json': source}).convert(geo('GSE1', ['GSE2']), out_type='json', enrichment='off')
    assert result.status == 'partial', result.to_dict()
    assert result.items[0].dataset_ids == ('GSE1',)


def test_native_project_edges_require_verified_owner_and_explicit_relation():
    from meta_standards_converter.converters.unified.families import recorded_neighbors
    value = package('ERP1', 'ERS1')
    value['source']['format'] = 'ENA'
    value['series']['accession'].append({'value': 'PRJEB1', 'database': 'GEO'})
    value['extensions'] = {'insdc': {'records': [
        {'provider': 'ena', 'kind': 'PROJECT', 'accession': 'PRJEB1', 'metadata': tree(ET.fromstring('''
        <PROJECT accession="PRJEB1"><RELATED_PROJECTS><RELATED_PROJECT>
          <CHILD_PROJECT accession="PRJEB2"/><PARENT_PROJECT accession="PRJEB3"/>
          <RELATED_PROJECT accession="PRJEB4"/>
        </RELATED_PROJECT></RELATED_PROJECTS></PROJECT>'''))},
        {'provider': 'ena', 'kind': 'PROJECT', 'accession': 'PRJEB99', 'metadata': tree(ET.fromstring('<PROJECT accession="PRJEB99"><CHILD_PROJECT accession="PRJEB5"/></PROJECT>'))}
    ]}}
    assert recorded_neighbors(MINiMLCodec().decode(value).package, 'ena') == ['PRJEB2', 'PRJEB3']
    assert recorded_neighbors(MINiMLCodec().decode(value).package, 'sra') == []


def test_ae_accession_keeps_x2x_signature():
    calls = []
    class AE:
        def convert(self, source, out=None, sdrf_sources=None):
            calls.append(source)
            return [MINiMLCodec().decode(package('E-MTAB-1', 'sample1')).package]
    result = Converter(services={'ae2json': AE()}).convert('E-MTAB-1', out_type='json', enrichment='off')
    assert result.status == 'complete', result.to_dict()
    assert calls == ['E-MTAB-1']


def test_native_project_expansion_keeps_multiple_studies_from_one_project():
    from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
    def native(study, project, relatives=()):
        value = package(study, 'ERS' + study[3:])
        value['series']['iid'] = study
        value['source']['format'] = 'ENA'
        value['series']['accession'].append({'value': project, 'database': 'GEO'})
        xml = '<PROJECT accession="' + project + '">' + ''.join('<CHILD_PROJECT accession="' + r + '"/>' for r in relatives) + '</PROJECT>'
        value['extensions'] = {'insdc': {'records': [{'provider': 'ena', 'kind': 'PROJECT', 'accession': project, 'metadata': tree(ET.fromstring(xml))}]}}
        return MINiMLCodec().decode(value).package
    calls = []
    class ENA:
        def convert(self, accession, **kwargs):
            calls.append(accession)
            # An umbrella accession resolves to its descendants; it need not be
            # asserted as an equivalent identifier of each descendant study.
            return ArchiveImportResult(accession, [StudyImportOutcome(study, 'PRJEB3', 'complete', native(study, 'PRJEB3', ['PRJEB1'])) for study in ['ERP2', 'ERP3']])
    result = Converter(services={'ena2json': ENA()}).convert(native('ERP1', 'PRJEB1', ['PRJEB2']), out_type='json', enrichment='off')
    assert result.status == 'complete', result.to_dict()
    assert result.items[0].dataset_ids == ('ERP1', 'ERP2', 'ERP3')
    assert calls == ['PRJEB2']


def test_ncbi_verified_hierarchy_query_uses_native_injected_source():
    from meta_standards_converter.converters.unified.families import StudyFamilies
    from types import SimpleNamespace
    value = package('SRP1', 'SRS1')
    value['series']['accession'].append({'value': 'PRJNA1', 'database': 'GEO'})
    value['extensions'] = {'insdc': {'records': [{'provider': 'sra', 'kind': 'Project', 'metadata': tree(ET.fromstring('<Project><ProjectID><ArchiveID id="1" accession="PRJNA1"/></ProjectID></Project>'))}]}}
    calls = []
    class Source:
        def publication_links(self, dbfrom, db, ids, name):
            calls.append(name)
            return {'1': ['2'] if name.endswith('u2d') else ['3']}
        def project_xml(self, uid, result):
            return ET.fromstring('<RecordSet><Project><ProjectID><ArchiveID id="' + uid + '" accession="PRJNA' + uid + '"/></ProjectID></Project></RecordSet>')
    context = SimpleNamespace(service=lambda key, factory: SimpleNamespace(source=Source()))
    families = StudyFamilies(context)
    p = MINiMLCodec().decode(value).package
    assert families._neighbors(p, 'sra') == (['PRJNA2', 'PRJNA3'], [])
    assert families._neighbors(p, 'sra') == (['PRJNA2', 'PRJNA3'], [])
    assert calls == ['bioproject_bioproject_u2d', 'bioproject_bioproject_d2u']
