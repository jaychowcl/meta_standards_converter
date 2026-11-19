# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""The legacy review never accepts incompatible file facts or authored renames."""
from copy import deepcopy
import pytest
from scripts.audit_integrity_legacy import completed_native_files, placeholder_names, added_cels


def test_completed_native_checksum_requires_exact_run_uri_and_consistent_digest():
    data = {'sample': [{'iid': 'S', 'sra_run': [{'run': 'ERR1', 'files': [
        {'uri': 'ftp://host/a.fastq', 'format': 'FASTQ', 'md5': 'a' * 32}]}]}]}
    old = [['Source Name', 'Comment[ENA_RUN]', 'Comment[FASTQ_URI]'], ['S', 'ERR1', 'ftp://host/a.fastq']]
    new = [old[0] + ['Comment[FASTQ_MD5]'], old[1] + ['a' * 32]]
    assert completed_native_files(data, old, new)[0] == old
    wrong = deepcopy(new); wrong[1][-1] = 'b' * 32
    with pytest.raises(AssertionError): completed_native_files(data, old, wrong)
    wrong = deepcopy(new); wrong[1][1] = 'ERR2'
    with pytest.raises(AssertionError): completed_native_files(data, old, wrong)
    data['sample'][0]['sra_run'][0]['files'].append({'uri': 'ftp://host/a.fastq', 'format': 'fastq', 'md5': 'b' * 32})
    with pytest.raises(AssertionError): completed_native_files(data, old, new)


def test_only_unauthored_blank_required_protocol_names_can_change():
    old = [['Protocol Name', 'P-old'], ['Protocol Type', 'sample collection protocol'],
           ['Protocol Term Source REF', 'EFO'], ['Protocol Term Accession Number', 'EFO_0005518'], ['Protocol Description', '']]
    new = deepcopy(old); new[0][1] = 'P-new'
    assert placeholder_names({'series': {}}, old, new)[1] == {'P-old': 'P-new'}
    with pytest.raises(AssertionError): placeholder_names({'series': {'protocols': [{'name': 'P-old'}]}}, old, new)
    old[-1][1] = new[-1][1] = 'Supplied collection protocol'
    with pytest.raises(AssertionError): placeholder_names({'series': {}}, old, new)


def test_restored_cel_must_belong_to_the_exact_sample():
    data = {'sample': [{'iid': 'S1', 'supplementary_data': [{'value': '\nftp://host/a.CEL.gz\n', 'type': 'CEL'}]},
                       {'iid': 'S2', 'supplementary_data': [{'value': 'ftp://host/b.CEL.gz', 'type': 'CEL'}]}]}
    old = [['Source Name', 'Scan Name'], ['S1', 'A']]
    new = [old[0] + ['Array Data File', 'Comment[ArrayExpress FTP file]'], old[1] + ['ftp://host/a.CEL.gz'] * 2]
    assert added_cels(data, old, new)[0] == old
    new[1][0] = 'S2'
    with pytest.raises(AssertionError): added_cels(data, old, new)
