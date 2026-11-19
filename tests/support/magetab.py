# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared magetab builders and deterministic collaborators."""

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
import csv, io
from copy import deepcopy
from tests.support.archive import legacy


def render(data):
    return AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package)


def definition(name, description="20 ug RNA’s preparation", **kwargs):
    return dict(name=name, description=description, type={"value": "library construction protocol"}, hardware=["machine"], **kwargs)


def parse(header, rows, **kwargs):
    return AEParser().parse(MAGETabInput(
        TextResource('x.idf.txt', 'Investigation Accession\tE-GEOD-1\nInvestigation Title\tExample\nSDRF File\tx.sdrf.txt\n', 'memory:idf'),
        (TextResource('x.sdrf.txt', '\n'.join('\t'.join(r) for r in [header, *rows]), 'memory:sdrf'),),
        'E-GEOD-1', 'accession', **kwargs,
    ))


def sample(title='', protocol='', code=None, iid='GSM1'):
    channel={'source': 'source', 'extract_protocol': protocol}
    if code is not None:
        channel['characteristics']=[{'name':'singlecell_type','value':code}]
    return {'iid':iid,'title':title,'type':'SRA','library_strategy':'RNA-Seq',
            'accession':[{'database':'GEO','value':iid}], 'channel':[channel]}


def import_tables(tables):
    idf=[];sdrf=next(r[1] for r in tables if r[0]=='SDRF File')
    def text(rows):
        out=io.StringIO();csv.writer(out,delimiter='\t').writerows(rows);return out.getvalue()
    for row in tables:idf.append(['SDRF File','study.sdrf.txt'] if row[0]=='SDRF File' else row)
    source=MAGETabInput(TextResource('study.idf.txt',text(idf),'recorded'),
        (TextResource('study.sdrf.txt',text(sdrf),'recorded'),),'recorded','local')
    parser=AEParser();value=parser.parse(source)
    return value,parser.warnings


def records(data):
    if hasattr(data,'to_mapping'):data=data.to_mapping()
    return [r for r in data.get('extensions',{}).get('insdc',{}).get('records',[]) if r.get('kind')=='source_conflict']


def file(uri, format='fastq', **values):
    return {'kind':'array_data_file','name':uri.rsplit('/',1)[-1], 'link':{'value':uri,'type':format},
            'comments':[{'name':k,'value':str(v)} for k,v in {'FORMAT':format,**values}.items()]}


def package():
    data=legacy();base=deepcopy(data['series']['assay_paths'][0]['steps']);base=[s for s in base if s['kind']!='array_data_file']
    paths=[]
    for f in [file('ftp://reads/one_1.fastq.gz',MD5='abc',BYTES=12),file('ftp://reads/one_2.fastq.gz'),
              file('ftp://reads/one_I1.fastq.gz'),file('ftp://reads/original.srf','SRF',ROLE='SUBMISSION_FILE'),
              file('ftp://reads/one_1.fastq.gz')]:
        paths.append({'steps':deepcopy(base)+[f]})
    data['series']['assay_paths']=paths
    return data


def projected(data):
    from meta_standards_converter.magetab.native_files import project_native_files
    project_native_files(data)
    return data['series']['assay_paths']


def comments(path):
    return [(c['name'],c['value']) for s in path['steps'] for c in s.get('comments',[])]
