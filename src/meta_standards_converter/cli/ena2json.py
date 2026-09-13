# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""ENA native import command."""
from .archive import parser_for, run_cli
from ..converters.ena2json import ENA2JSONConverter


def _parser():
    return parser_for('ENA')


def main(argv=None):
    return run_cli(_parser(), ENA2JSONConverter, argv)


if __name__ == '__main__':
    raise SystemExit(main())
