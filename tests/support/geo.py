# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared geo builders and deterministic collaborators."""

from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.magetab.constructor import AEConstructor


def miniml_body(body: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML" version="0.5.4">
{body}
</MINiML>
"""


def geo_package():
    samples = ''.join(f'''<Sample iid="GSM{i}"><Title>Sample  {i}</Title>
      <Library-Strategy>RNA-Seq</Library-Strategy><Library-Source>transcriptomic</Library-Source>
      <Channel><Source>  whole  embryo </Source><Organism>Drosophila melanogaster</Organism>
      <Molecule>total RNA</Molecule><Characteristics tag="tissue">whole  embryo</Characteristics>
      <Characteristics tag="developmental stage">stage {i}</Characteristics>
      <Extract-Protocol>10x Chromium Single Cell 3' v3.1.  Author's method.</Extract-Protocol>
      </Channel><Supplementary-Data type="TXT">https://example.test/processed{i}.txt</Supplementary-Data>
      </Sample>''' for i in (1, 2))
    return GEOParser().parse(f'''<MINiML version="0.5.4">{samples}
      <Series iid="GSE1"><Title>Study</Title><Summary>Summary.</Summary><Overall-Design>Design.</Overall-Design>
      <Status><Submission-Date>May 01 2023</Submission-Date><Release-Date>May 02 2023</Release-Date></Status>
      <Sample-Ref ref="GSM1"/><Sample-Ref ref="GSM2"/></Series></MINiML>''')[0].to_mapping()


def rows(data, handler='single_cell_sequencing'):
    return AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package, platform_handler=handler)


def row(table, label):
    return next(r for r in table if r[0] == label)
