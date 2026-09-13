# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""SRA native import command."""
from .archive import parser_for, run_cli
from ..converters.sra2json import SRA2JSONConverter


def _parser():
    return parser_for('SRA')


def main(argv=None):
    return run_cli(_parser(), SRA2JSONConverter, argv)


if __name__ == '__main__':
    raise SystemExit(main())
