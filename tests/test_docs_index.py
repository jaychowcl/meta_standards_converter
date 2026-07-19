# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import importlib
import re
import subprocess
import unittest
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "docs" / "index.md"
CODEBASE = ROOT / "docs" / "codebase.md"
README = ROOT / "README.md"
PYPROJECT = ROOT / "pyproject.toml"

AUTHOR_TEXT = (
    "Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026\n"
    "https://github.com/jaychowcl\n"
    "https://saezlab.org\n"
    "https://www.ebi.ac.uk/about/teams/functional-genomics/"
)
PYTHON_AUTHOR_HEADER = "\n".join(
    [
        "# =============================================================================",
        "# Authors",
        "#",
        *[f"# {line}" for line in AUTHOR_TEXT.splitlines()],
        "# =============================================================================",
        "",
    ]
)
HASH_AUTHOR_HEADER = PYTHON_AUTHOR_HEADER
HTML_AUTHOR_HEADER = "\n".join(
    [
        "<!--",
        "Authors",
        "",
        *AUTHOR_TEXT.splitlines(),
        "-->",
        "",
    ]
)
README_AUTHORS_LINE = (
    "Created by [jaychowcl](https://github.com/jaychowcl) @ "
    "[Saez-Rodriguez Group](https://saezlab.org) & "
    "[EMBL-EBI Functional Genomics Team](https://www.ebi.ac.uk/about/teams/functional-genomics/) "
    "on May 2026"
)


class DocsIndexTests(unittest.TestCase):
    def test_index_uses_header_references_not_line_ranges(self):
        index_text = INDEX.read_text()

        self.assertNotIn("lines:", index_text)
        self.assertNotIn("line range", index_text.lower())
        self.assertNotIn("sed -n '<start>,<end>p'", index_text)

    def test_every_index_anchor_exists_and_points_to_a_header(self):
        index_text = INDEX.read_text()
        anchors = re.findall(r"^\s+anchor:\s+([a-z0-9-]+)\s*$", index_text, re.MULTILINE)
        codebase_lines = CODEBASE.read_text().splitlines()

        self.assertTrue(anchors)

        for anchor in anchors:
            anchor_line = f'<a id="{anchor}"></a>'
            self.assertIn(anchor_line, codebase_lines)
            line_index = codebase_lines.index(anchor_line)
            self.assertLess(line_index + 1, len(codebase_lines))
            self.assertRegex(codebase_lines[line_index + 1], r"^#{2,6} ", anchor)

    def test_readme_contains_requested_sections_in_order(self):
        readme_text = README.read_text(encoding="utf-8")
        expected_headings = [
            "# meta_standards_converter",
            "## Description",
            "## Installation",
            "### Requirements",
            "## Quickstart",
            "### CLI quickstart",
            "### Python API quickstart",
            "### Docker quickstart",
            "### Rootless Docker Compose quickstart",
            "### Inputs & Outputs",
            "## Guide",
            "### CLI",
            "### Python API",
            "### Docker",
            "### Rootless Docker Compose",
            "### Code flow",
            "## Docs",
            "## Authors",
        ]

        lines = readme_text.splitlines()
        positions = []
        for heading in expected_headings:
            self.assertIn(heading, lines)
            positions.append(lines.index(heading))

        self.assertEqual(sorted(positions), positions)

    def test_readme_quickstarts_link_to_each_interface_guide(self):
        readme_text = README.read_text(encoding="utf-8")

        for label, anchor in (
            ("CLI guide", "cli"),
            ("Python API guide", "python-api"),
            ("Docker guide", "docker"),
            ("Rootless Docker Compose guide", "rootless-docker-compose"),
        ):
            self.assertIn(f"[{label}](#{anchor})", readme_text)

    def test_readme_cli_guide_documents_every_parser_argument(self):
        readme_text = README.read_text(encoding="utf-8")
        modules = {
            command: importlib.import_module(f"meta_standards_converter.cli.{command}")
            for command in ("geo2ae", "geo2json", "json2ae", "ae2json", "json2h5ad")
        }

        for command, module in modules.items():
            match = re.search(
                rf"^#### `{command}`\s*$\n(?P<section>.*?)(?=^#### |^### |\Z)",
                readme_text,
                re.MULTILINE | re.DOTALL,
            )
            self.assertIsNotNone(match, command)
            section = match.group("section")
            for action in module._parser()._actions:
                if action.option_strings:
                    for option in action.option_strings:
                        self.assertIn(f"`{option}`", section, f"{command}: {option}")
                else:
                    self.assertIn(f"`{action.dest}`", section, f"{command}: {action.dest}")

    def test_readme_links_to_docs(self):
        readme_text = README.read_text(encoding="utf-8")

        self.assertIn("[Codebase docs](docs/codebase.md)", readme_text)
        self.assertIn("[Docs index](docs/index.md)", readme_text)

    def test_readme_documents_all_console_scripts(self):
        readme_text = README.read_text(encoding="utf-8")
        with PYPROJECT.open("rb") as handle:
            pyproject = tomllib.load(handle)

        for script_name in pyproject["project"]["scripts"]:
            self.assertIn(f"`{script_name}`", readme_text)

    def test_tracked_commentable_files_have_canonical_author_headers(self):
        tracked_files = subprocess.check_output(
            ["git", "ls-files"],
            cwd=ROOT,
            text=True,
        ).splitlines()
        external_or_noncommentable = {
            "LICENSE",
            "docs/MAGE-TABv1.1_2011_07_28.pdf",
            "docs/MINiML.xsd",
            "tests/GSE328265_family.xml",
        }
        hash_comment_files = {
            ".dockerignore",
            ".gitignore",
            "Dockerfile",
            "compose.yaml",
            "pyproject.toml",
            "requirements.txt",
        }

        for tracked_file in tracked_files:
            if tracked_file in external_or_noncommentable or tracked_file == "README.md":
                continue

            path = ROOT / tracked_file
            text = path.read_text(encoding="utf-8")
            if tracked_file.endswith(".py") or tracked_file in hash_comment_files:
                self.assertTrue(text.startswith(HASH_AUTHOR_HEADER), tracked_file)
            elif tracked_file.endswith(".sh"):
                self.assertTrue(
                    text.startswith("#!/usr/bin/env bash\n" + HASH_AUTHOR_HEADER),
                    tracked_file,
                )
            elif tracked_file.endswith(".md"):
                self.assertTrue(text.startswith(HTML_AUTHOR_HEADER), tracked_file)
            else:
                self.fail(f"Unhandled tracked file for author header policy: {tracked_file}")

    def test_readme_authors_section_uses_linked_author_line(self):
        readme_text = README.read_text(encoding="utf-8")

        self.assertFalse(readme_text.startswith("<!--\nAuthors"))
        self.assertIn(f"## Authors\n\n{README_AUTHORS_LINE}", readme_text)


if __name__ == "__main__":
    unittest.main()
