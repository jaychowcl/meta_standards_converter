# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Keep CLI logging configuration scoped to a test operation."""
from contextlib import contextmanager
import logging


@contextmanager
def preserve_package_logging():
    logger = logging.getLogger("meta_standards_converter")
    handlers, level, propagate = list(logger.handlers), logger.level, logger.propagate
    try:
        yield
    finally:
        for handler in list(logger.handlers):
            if handler not in handlers:
                handler.close()
        logger.handlers[:] = handlers
        logger.setLevel(level)
        logger.propagate = propagate
