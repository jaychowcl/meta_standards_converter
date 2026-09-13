"""ENA native import command."""
from .archive import parser_for, run_cli
from ..converters.ena2json import ENA2JSONConverter


def _parser():
    return parser_for('ENA')


def main(argv=None):
    return run_cli(_parser(), ENA2JSONConverter, argv)


if __name__ == '__main__':
    raise SystemExit(main())
