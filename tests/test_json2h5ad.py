# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import os
import sys
import tempfile
import unittest


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.json2h5ad import json2h5ad  # noqa: E402


class TestJSON2H5ADConverter(unittest.TestCase):
    def test_existing_json_raises_not_implemented_with_useful_message(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = os.path.join(tmpdir, "GSE1.json")
            with open(json_path, "w", encoding="utf-8") as handle:
                json.dump([{"series": {"accession": [{"value": "GSE1"}]}}], handle)

            with self.assertRaises(NotImplementedError) as error:
                json2h5ad().convert(json_path=json_path, out=tmpdir)

        message = str(error.exception)
        self.assertIn("json2h5ad", message)
        self.assertIn("matrix fetching/parsing", message)
        self.assertIn("AnnData writing", message)

    def test_missing_json_path_raises_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            json2h5ad().convert(json_path="/missing/GSE1.json", out=".")

    def test_none_out_does_not_require_output_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = os.path.join(tmpdir, "GSE1.json")
            with open(json_path, "w", encoding="utf-8") as handle:
                json.dump([], handle)

            with self.assertRaises(NotImplementedError):
                json2h5ad().convert(json_path=json_path, out=None)


if __name__ == "__main__":
    unittest.main()
