"""SRA native import command."""
from .archive import parser_for, run_cli
from ..converters.sra2json import SRA2JSONConverter


def _parser():
    return parser_for('SRA')


def main(argv=None):
    return run_cli(_parser(), SRA2JSONConverter, argv)


if __name__ == '__main__':
    raise SystemExit(main())
