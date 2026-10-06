# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operation-local immutable construction checkpoints for native imports."""
from dataclasses import dataclass, replace
import hashlib
import json
import xml.etree.ElementTree as ET

from ..conversion_logging import measured
from . import insdc_support
from .codec import MINiMLCodec, MINiMLCompatibilityError


@dataclass(frozen=True)
class _NativeCheckpoint:
    provider: str
    fingerprint: str
    construction: str
    package: object


@measured('indexing')
def _fingerprint(records):
    digest=hashlib.sha256()
    digest.update(json.dumps([records.seed.study,records.seed.primary]).encode())
    for root in records.xml:
        for node in root:
            if node.tag in {'STUDY','PROJECT','SAMPLE','EXPERIMENT','RUN','EXPERIMENT_PACKAGE'}:
                digest.update(ET.tostring(node))
    digest.update(json.dumps({k:records.indexed[k] for k in ('study','sample','read_experiment','read_run')
                             if k in records.indexed},sort_keys=True).encode())
    return digest.hexdigest()


@measured('validation')
def _valid_native(package):
    # finish() already constructed this package through the canonical mapping
    # decoder. Re-run strict validation without serializing/decoding it again.
    diagnostics=package.validate()
    structural=tuple(item for item in diagnostics
                     if not MINiMLCodec._is_source_compatible_sample_title_warning(item,package))
    if structural:
        raise MINiMLCompatibilityError(structural)
    return package


def _checkpoint(parser,records):
    from .ena_parser import ENAParser
    from .sra_parser import SRAParser
    if type(parser) not in (ENAParser,SRAParser):
        return None
    provider='ena' if type(parser) is ENAParser else 'sra'
    construction=parser._construct(records)
    frozen=json.dumps(construction,separators=(',',':'))
    package=_valid_native(insdc_support.finish(provider,records,*construction))
    return _NativeCheckpoint(provider,_fingerprint(records),frozen,package)


@measured('supplementary_projection')
def _complete(checkpoint,parser,records):
    if _fingerprint(records)!=checkpoint.fingerprint:
        return parser.parse(records)
    construction=json.loads(checkpoint.construction)
    return insdc_support.finish(checkpoint.provider,records,*construction)


def _projection_failed(checkpoint,error):
    from copy import deepcopy
    from .model import _plain
    extensions=_plain(checkpoint.package.extensions)
    extensions.setdefault('insdc',{'version':'2.0','records':[]})['records'].append({
        'provider':checkpoint.provider,'kind':'retrieval_coverage',
        'accession':checkpoint.package.series.iid,'metadata':{'complete':False,
            'unfinished':[{'stage':'native_supplementary','reason':
                          f'{type(error).__name__}: invalid projection; valid core retained'}]}})
    package=replace(checkpoint.package,extensions=extensions)
    evidence=getattr(checkpoint.package,'_archive_source_records',None)
    if evidence is not None:
        object.__setattr__(package,'_archive_source_records',deepcopy(evidence))
    return _valid_native(package)
