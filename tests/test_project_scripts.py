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
    def test_ae2json_console_script_is_registered(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(
            "meta_standards_converter.cli.ae2json:main",
            pyproject["project"]["scripts"]["ae2json"],
        )

    def test_json2ae_console_script_is_registered(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(
            "meta_standards_converter.cli.json2ae:main",
            pyproject["project"]["scripts"]["json2ae"],
        )

    def test_json2h5ad_console_script_is_registered(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(
            "meta_standards_converter.cli.json2h5ad:main",
            pyproject["project"]["scripts"]["json2h5ad"],
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
