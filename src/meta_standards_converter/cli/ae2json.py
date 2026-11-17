# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""ae2json: named entrypoint for the shared metadata conversion pipeline."""
from .convert import main as convert_main


def main(argv=None):
    return convert_main(argv, command='ae2json', output_type='json', default_in_type=None)


if __name__ == "__main__":
    raise SystemExit(main())
