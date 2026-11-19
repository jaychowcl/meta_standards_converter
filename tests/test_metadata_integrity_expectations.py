# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reviewed expectations retain sealed origins and independent publication facts."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def test_reviewed_expectations_leave_historical_bytes_unchanged():
    root = Path(__file__).parents[1]
    reviewed = root / 'tests/fixtures/edge_cases/metadata-integrity/reviewed'
    manifest = json.loads((reviewed / 'review.json').read_text())
    for item in manifest['files']:
        assert hashlib.sha256((root / item['original']).read_bytes()).hexdigest() == item['original_sha256']
        assert hashlib.sha256((reviewed / item['derivative']).read_bytes()).hexdigest() == item['sha256']


def test_new_json_original_filename_and_checksum_are_supported_by_sra():
    fixtures = Path(__file__).parent / 'fixtures'
    xml = ET.fromstring((fixtures / 'studies/GSE328265/provider_responses/sra.xml').read_bytes())
    originals = {(n.get('filename'), n.get('md5')) for n in xml.iter('SRAFile') if n.get('supertype') == 'Original'}
    data = json.loads((fixtures / 'edge_cases/metadata-integrity/reviewed/GSE328265/geo2json/GSE328265.json').read_text())[0]
    record = data['sample'][0]['sra_run'][0]
    assert (record['submitted_file_name'], record['md5']) in originals
    assert originals <= {(f.get('filename'), f.get('md5')) for f in record['files']}
