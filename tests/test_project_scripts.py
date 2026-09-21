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

import tomllib


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


class TestProjectScripts(unittest.TestCase):
    def test_all_console_scripts_are_registered(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        for name in (
            "sra2json",
            "ena2json",
            "miniml-migrate",
            "msc-convert",
            "ae2json",
            "geo2ae",
            "geo2json",
            "json2ae",
            "json2h5ad",
            "json2tsv",
            "json2obs",
        ):
            with self.subTest(name=name):
                self.assertEqual(
                    f"meta_standards_converter.cli.{'convert' if name == 'msc-convert' else name.replace('-', '_')}:main",
                    pyproject["project"]["scripts"][name],
                )

    def test_test_extra_contains_the_canonical_pytest_stack(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(
            [
                "pytest>=9.1.1,<10",
                "anndata>=0.13.4,<1",
                "h5py>=3.16.0,<4",
                "numpy>=2.5.3,<3",
                "pandas>=3.0.6,<4",
                "scanpy>=1.12.4,<2",
                "scipy>=1.18.1,<2",
            ],
            pyproject["project"]["optional-dependencies"]["test"],
        )

    def test_runtime_and_scientific_dependencies_have_tested_major_bounds(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            project = tomllib.load(handle)["project"]

        self.assertEqual(
            ["python-dateutil>=2.9.0.post0,<3", "requests>=2.34.2,<3"],
            project["dependencies"],
        )
        self.assertEqual(
            [
                "anndata>=0.13.4,<1",
                "h5py>=3.16.0,<4",
                "numpy>=2.5.3,<3",
                "pandas>=3.0.6,<4",
                "scanpy>=1.12.4,<2",
                "scipy>=1.18.1,<2",
            ],
            project["optional-dependencies"]["h5ad"],
        )

    def test_modern_build_metadata_and_package_discovery(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as handle:
            config = tomllib.load(handle)
        project = config["project"]
        self.assertEqual(project["version"], "8.0.0")
        self.assertEqual(project["requires-python"], ">=3.12")
        self.assertEqual(project["license"], "GPL-3.0-only")
        self.assertEqual(project["license-files"], ["LICENSE"])
        self.assertFalse(any(c.startswith("License ::") for c in project["classifiers"]))
        self.assertEqual(config["build-system"]["requires"], ["setuptools>=77.0.3"])
        self.assertEqual(config["tool"]["setuptools"]["packages"]["find"]["include"],
                         ["meta_standards_converter*"])

    def test_requirements_installs_project_from_single_dependency_authority(self):
        with open(os.path.join(ROOT, "requirements.txt")) as handle:
            requirements = [line.strip() for line in handle
                            if line.strip() and not line.startswith("#")]
        self.assertEqual(requirements, ["."])
        self.assertFalse(os.path.exists(os.path.join(ROOT, "dependency-provenance")))


if __name__ == "__main__":
    unittest.main()
