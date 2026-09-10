# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Regression for CLI configuration leaking closed capture streams across tests."""
from io import StringIO
import logging
from types import SimpleNamespace
from meta_standards_converter.cli.common import configure_logging
from tests.support.logging_state import preserve_package_logging


def test_cli_logging_state_is_restored_even_after_an_error():
    logger = logging.getLogger("meta_standards_converter")
    original = (list(logger.handlers), logger.level, logger.propagate)
    try:
        with preserve_package_logging():
            configure_logging(SimpleNamespace(quiet=False, verbose=1, log_file=None), stream=StringIO())
            assert logger.propagate is False
            assert logger.handlers != original[0]
            raise RuntimeError("test operation interrupted")
    except RuntimeError as error:
        assert str(error) == "test operation interrupted"
    assert (logger.handlers, logger.level, logger.propagate) == original
