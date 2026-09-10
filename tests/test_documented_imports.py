# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Validate public imports in documentation without executing example workflows."""
import ast
import importlib
import re
from pathlib import Path


def test_documented_python_msc_imports_resolve():
    root = Path(__file__).resolve().parents[1]
    checked = []
    for path in [root / "README.md", *sorted((root / "docs").glob("*.md"))]:
        for block in re.finditer(r"```(?:python|py)\s*\n(.*?)```", path.read_text(), re.S):
            # Parse import statements individually; examples may contain shell placeholders.
            for line in block.group(1).splitlines():
                if not line.startswith("from meta_standards_converter"):
                    continue
                tree = ast.parse(line)
                node = tree.body[0]
                module = importlib.import_module(node.module)
                for item in node.names:
                    assert item.name != "*", "Document explicit supported imports"
                    assert getattr(module, item.name) is not None, (path, line)
                checked.append((path, line))
    assert checked, "Documentation must contain executable public Python imports"
