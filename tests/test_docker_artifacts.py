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
        self.assertIn("RUN pip install --no-cache-dir '.[h5ad]'", content)
        self.assertIn('CMD ["geo2ae", "--help"]', content)

    def test_dockerfile_installs_pinned_nfcore_runtime(self):
        content = (ROOT / "Dockerfile").read_text(encoding="utf-8")

        self.assertIn("ARG NEXTFLOW_VERSION=26.04.2", content)
        self.assertIn("ARG DOCKER_CLI_VERSION=29.6.2", content)
        self.assertIn("openjdk-17-jre-headless", content)
        self.assertIn("NEXTFLOW_VERSION", content)
        self.assertIn("nextflow -version", content)
        self.assertIn("docker --version", content)

    def test_compose_uses_only_the_rootless_docker_socket(self):
        compose = ROOT / "compose.yaml"

        self.assertTrue(compose.exists())
        content = compose.read_text(encoding="utf-8")
        self.assertIn("${ROOTLESS_DOCKER_SOCKET:?set ROOTLESS_DOCKER_SOCKET}", content)
        self.assertIn("/run/rootless-docker/docker.sock", content)
        self.assertNotIn("/var/run/docker.sock", content)
        self.assertIn("META_STANDARDS_REQUIRE_ROOTLESS_DOCKER: \"1\"", content)
        self.assertIn("DOCKER_HOST: unix:///run/rootless-docker/docker.sock", content)

    def test_compose_hardens_and_limits_the_runtime_mount(self):
        content = (ROOT / "compose.yaml").read_text(encoding="utf-8")

        self.assertIn("source: ${JSON2H5AD_OUT:?set JSON2H5AD_OUT}", content)
        self.assertIn("target: ${JSON2H5AD_OUT:?set JSON2H5AD_OUT}", content)
        self.assertIn("read_only: true", content)
        self.assertIn("no-new-privileges:true", content)
        self.assertIn("cap_drop:", content)
        self.assertIn("- ALL", content)

    def test_rootless_runner_scripts_are_present_and_parse(self):
        scripts = [
            ROOT / "scripts" / "provision-rootless-json2h5ad.sh",
            ROOT / "scripts" / "json2h5ad-compose.sh",
        ]

        for script in scripts:
            with self.subTest(script=script.name):
                self.assertTrue(script.exists())
                completed = __import__("subprocess").run(
                    ["bash", "-n", str(script)],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(0, completed.returncode, completed.stderr)

        provision = scripts[0].read_text(encoding="utf-8")
        self.assertIn("nfcore-runner", provision)
        self.assertIn("dockerd-rootless-setuptool.sh install", provision)
        self.assertIn("loginctl enable-linger", provision)
        self.assertIn("SecurityOptions", provision)

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
