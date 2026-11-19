# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Public-API baseline backchecks for owned unknown links and duplicate views."""
from copy import deepcopy
import json
from pathlib import Path
from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.magetab.constructor import AEConstructor

SOURCE = json.loads((Path(__file__).parents[1] / 'fixtures/edge_cases/metadata-integrity/source-cases.json').read_text())


def table(data):
    document = AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package)
    return next(r[1] for r in document if r[0] == 'SDRF File')


def test_unknown_link_keeps_source_ownership():
    data = GEOParser().parse(SOURCE['geo']['GSE31356'])[0].to_mapping()
    uri = 'https://example.org/unknown.txt'
    data['sample'][0]['supplementary_data'] = [{'value': uri, 'type': 'TXT'}]
    rows = table(data)
    assert 'Comment[Unclassified File]' in rows[0]
    assert rows[1][rows[0].index('Comment[Unclassified File]')] == uri
    assert rows[1][rows[0].index('Source Name')] == data['sample'][0]['iid']


def test_duplicate_source_occurrences_share_one_download_presentation():
    data = GEOParser().parse(SOURCE['geo']['GSE52778'])[0].to_mapping()
    file = {'uri': 'ftp://example.org/SRR1039508_1.fastq.gz', 'filename': 'SRR1039508_1.fastq.gz',
            'md5': 'a' * 32, 'bytes': '100'}
    data['sample'][0]['sra_run'] = [{'run': 'SRR1039508', 'experiment': 'SRX384345', 'study': 'SRP033351',
                                  'library_layout': 'PAIRED', 'fastq_files': [file, deepcopy(file)]}]
    before = deepcopy(data)
    rows = table(data)
    assert len(rows) - 1 == 1
    assert data == before and len(data['sample'][0]['sra_run'][0]['fastq_files']) == 2
