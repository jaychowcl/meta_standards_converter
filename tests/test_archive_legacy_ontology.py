# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import pytest
from meta_standards_converter.miniml.archive_entities import declare_ontologies
from tests.test_native_archive_enrichment import native
from tests.test_protocol_export import render


@pytest.mark.parametrize('uri,expected',[
 ('http://purl.org/obo/owl/NCBITaxon#NCBITaxon_4932','NCBITaxon'),
 ('https://purl.org/obo/owl/UBERON#UBERON_0000955','UBERON'),
 ('http://purlXorg/obo/owl/NCBITaxon#NCBITaxon_4932','EFO'),
 ('http://purl.org/obo/owl/NCBITaxon#UBERON_4932','EFO'),
 ('http://purlXobolibraryYorg/obo/NCBITaxon_4932','EFO')])
def test_only_explicit_matching_ontology_namespaces_are_corrected(uri,expected):
    data=native().to_mapping();group={'name':'organism','value':'yeast','term_source_ref':'EFO','term_accession_number':uri}
    data['sample'][0]['channel'][0]['characteristics']=[deepcopy(group)]
    declare_ontologies(data)
    assert data['sample'][0]['channel'][0]['characteristics'][0]['term_source_ref']==expected
    assert data['sample'][0]['channel'][0]['characteristics'][0]['term_accession_number']==uri
    if expected!='EFO':assert any(r['kind']=='annotation' and r['metadata']['annotation']==group for r in data['extensions']['insdc']['records'])
    before=deepcopy(data);declare_ontologies(data);assert data==before


def test_saved_native_export_repairs_legacy_term_source_without_input_mutation():
    data=native().to_mapping();uri='http://purl.org/obo/owl/NCBITaxon#NCBITaxon_4932'
    group={'name':'organism','value':'yeast','term_source_ref':'EFO','term_accession_number':uri}
    for path in data['series']['assay_paths']:path['steps'][0]['characteristics']=[deepcopy(group)]
    before=deepcopy(data);output=render(data);assert data==before
    table=next(r[1] for r in output if r[0]=='SDRF File')
    for row in table[1:]:
        i=row.index(uri);assert table[0][i]=='Term Accession Number';assert row[i-1]=='NCBITaxon'
