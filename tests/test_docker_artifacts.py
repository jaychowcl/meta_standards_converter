# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class DockerArtifactsTest(unittest.TestCase):
    def test_dockerfile_installs_local_cli_package(self):
        dockerfile = ROOT / "Dockerfile"

        self.assertTrue(dockerfile.exists())
        content = dockerfile.read_text(encoding="utf-8")

        self.assertIn("FROM python:3.12-slim", content)
        self.assertIn("WORKDIR /app", content)
        self.assertIn("COPY pyproject.toml README.md LICENSE ./", content)
        self.assertIn("COPY src ./src", content)
        self.assertIn("RUN pip install --no-cache-dir .", content)
        self.assertIn('CMD ["geo2ae", "--help"]', content)

    def test_dockerignore_excludes_local_and_generated_files(self):
        dockerignore = ROOT / ".dockerignore"

        self.assertTrue(dockerignore.exists())
        patterns = {
            line.strip()
            for line in dockerignore.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }

        expected_patterns = {
            ".git",
            ".gitignore",
            ".env",
            ".env*",
            "__pycache__/",
            "*.py[cod]",
            "*.egg-info/",
            ".pytest_cache/",
            ".mypy_cache/",
            ".ruff_cache/",
            ".venv/",
            "venv/",
            "output/",
            ".out/",
            ".vscode/",
            ".codex/",
            ".agents/",
            ".dev/",
            "*.log",
        }
        self.assertTrue(expected_patterns.issubset(patterns))


if __name__ == "__main__":
    unittest.main()
