# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import sys
import unittest

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


class TestProjectScripts(unittest.TestCase):
    def test_all_console_scripts_are_registered(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        for name in (
            "ae2json",
            "geo2ae",
            "geo2json",
            "json2ae",
            "json2h5ad",
            "json2tsv",
            "json2csv",
        ):
            with self.subTest(name=name):
                self.assertEqual(
                    f"meta_standards_converter.cli.{name}:main",
                    pyproject["project"]["scripts"][name],
                )

    def test_test_extra_contains_the_canonical_pytest_stack(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(
            [
                "pytest>=8.2,<9",
                "pytest-subtests>=0.14,<1",
                "anndata>=0.10.8",
                "h5py>=3.10.0",
                "numpy>=1.26.0",
                "pandas>=2.1.0",
                "scanpy>=1.10.0",
                "scipy>=1.11.0",
            ],
            pyproject["project"]["optional-dependencies"]["test"],
        )

    def test_project_license_classifier_matches_gplv3_license_file(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        classifiers = pyproject["project"]["classifiers"]
        self.assertIn(
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            classifiers,
        )
        self.assertNotIn("License :: OSI Approved :: Apache Software License", classifiers)


if __name__ == "__main__":
    unittest.main()
