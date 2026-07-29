# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import importlib
import ast
import re
import subprocess
import unittest
from pathlib import Path

from meta_standards_converter.ae_handlers.ae_constructor import PLATFORM_HANDLER_KEYS

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
    "[EMBL-EBI Functional Genomics Team]"
    "(https://www.ebi.ac.uk/about/teams/functional-genomics/) on May 2026"
)
README_LOGO = (
    '<img width="250" height="250" alt="image" '
    'src="https://github.com/user-attachments/assets/'
    '51b52963-19de-4f67-8977-072b409dae19" />'
)
CANONICAL_CODEBASE_ANCHORS = (
    "architecture",
    "system-context-and-boundaries",
    "architectural-decisions",
    "design-invariants-and-expectations",
    "component-relationships-and-data-flow",
    "entrypoints-and-interfaces",
    "orchestrators-and-core-types",
    "public-api-reference",
    "principal-workflows",
    "extension-and-change-guidance",
)
LEGACY_CODEBASE_ANCHORS = (
    "project-purpose-and-layout",
    "runtime-behavior",
    "end-to-end-geo2ae-flow",
    "end-to-end-json2ae-flow",
    "end-to-end-ae2json-flow",
    "json2h5ad-flow",
    "json2tabular-flow",
    "rootless-json2h5ad-runtime",
    "public-api-and-callable-reference",
    "maintenance-notes",
    "test-plan",
)
CLI_COMMANDS = (
    "geo2ae",
    "geo2json",
    "json2ae",
    "ae2json",
    "json2h5ad",
    "json2tsv",
    "json2csv",
)


class DocsIndexTests(unittest.TestCase):
    def test_codebase_has_substantive_canonical_sections(self):
        codebase_text = CODEBASE.read_text(encoding="utf-8")

        for index, anchor in enumerate(CANONICAL_CODEBASE_ANCHORS):
            start = codebase_text.index(f'<a id="{anchor}"></a>')
            if index + 1 < len(CANONICAL_CODEBASE_ANCHORS):
                end = codebase_text.index(
                    f'<a id="{CANONICAL_CODEBASE_ANCHORS[index + 1]}"></a>',
                    start,
                )
            else:
                end = len(codebase_text)
            section = codebase_text[start:end]
            self.assertRegex(section, rf'<a id="{anchor}"></a>\n## ')
            prose = re.sub(r"<[^>]+>|[#`*|:_-]", " ", section)
            self.assertGreater(len(prose.split()), 20, anchor)

    def test_index_routes_canonical_sections_with_purpose_and_keywords(self):
        index_text = INDEX.read_text(encoding="utf-8")

        for anchor in CANONICAL_CODEBASE_ANCHORS:
            route = re.search(
                rf"- id: {re.escape(anchor)}\n"
                rf"  title: .+\n"
                rf"  anchor: {re.escape(anchor)}\n"
                rf"  purpose: .+\n"
                rf"  keywords: .+",
                index_text,
            )
            self.assertIsNotNone(route, anchor)

    def test_legacy_codebase_anchors_are_preserved(self):
        codebase_text = CODEBASE.read_text(encoding="utf-8")

        for anchor in LEGACY_CODEBASE_ANCHORS:
            self.assertIn(f'<a id="{anchor}"></a>', codebase_text, anchor)

    def test_architectural_decisions_are_evidence_safe(self):
        codebase_text = CODEBASE.read_text(encoding="utf-8")
        start = codebase_text.index('<a id="architectural-decisions"></a>')
        end = codebase_text.index('<a id="design-invariants-and-expectations"></a>')
        section = codebase_text[start:end]
        records = re.findall(
            r"^### AD-\d{3}: .+?(?=^### AD-|\Z)",
            section,
            re.MULTILINE | re.DOTALL,
        )

        self.assertTrue(records)
        for record in records:
            status = re.search(r"\*\*Status:\*\* (Documented|Observed)", record)
            self.assertIsNotNone(status, record)
            self.assertRegex(record, r"\*\*Decision:\*\* .+")
            self.assertRegex(record, r"\*\*Rationale:\*\* .+")
            self.assertRegex(record, r"\*\*Consequences:\*\* .+")
            self.assertRegex(record, r"\*\*Affected components:\*\* .+")
            self.assertRegex(record, r"\*\*Evidence:\*\* .+")
            if status.group(1) == "Observed":
                self.assertIn("**Rationale:** Not documented.", record)

    def test_public_api_reference_covers_formal_exports_and_public_symbols(self):
        codebase_text = CODEBASE.read_text(encoding="utf-8")
        exports = importlib.import_module("meta_standards_converter.converters").__all__

        for symbol in exports:
            self.assertRegex(codebase_text, rf"`(?:[^`]*\.)?{re.escape(symbol)}`")

        source_root = ROOT / "src" / "meta_standards_converter"
        public_symbols = set()
        for source_path in source_root.rglob("*.py"):
            if "__pycache__" in source_path.parts:
                continue
            module = ".".join(source_path.relative_to(ROOT / "src").with_suffix("").parts)
            tree = ast.parse(source_path.read_text(encoding="utf-8"))
            for node in tree.body:
                if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                    if not node.name.startswith("_"):
                        public_symbols.add(f"{module}.{node.name}")

        self.assertTrue(public_symbols)
        for qualified_name in sorted(public_symbols):
            self.assertIn(f"`{qualified_name}`", codebase_text, qualified_name)

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
            "### Configuration",
            "### CLI",
            "### Python API",
            "### Docker",
            "### Rootless Docker Compose",
            "### Code flow",
            "## Docs",
            "## Authors",
        ]

        lines = readme_text.splitlines()
        self.assertTrue(
            readme_text.startswith(f"{README_LOGO}\n\n# meta_standards_converter\n"),
            "README must preserve the historical logo before the H1",
        )
        positions = []
        for heading in expected_headings:
            self.assertIn(heading, lines)
            positions.append(lines.index(heading))

        self.assertEqual(sorted(positions), positions)
        self.assertNotIn("## Configuration", lines)

    def test_readme_configuration_documents_platform_handler_hierarchy(self):
        readme_text = README.read_text(encoding="utf-8")
        match = re.search(
            r"^#### Platform handlers\s*$\n(?P<section>.*?)(?=^#### |^### |\Z)",
            readme_text,
            re.MULTILINE | re.DOTALL,
        )

        self.assertIsNotNone(match)
        section = match.group("section")
        self.assertIn("```mermaid", section)
        self.assertIn("--list-platform-handlers", section)
        for handler_key in PLATFORM_HANDLER_KEYS:
            self.assertIn(handler_key, section)

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
            for command in CLI_COMMANDS
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

    def test_readme_documents_all_seven_conversion_workflows(self):
        readme_text = README.read_text(encoding="utf-8")

        for command in CLI_COMMANDS:
            self.assertIn(f"`{command}`", readme_text)

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
