# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import unittest
import re
import os
import socket
import subprocess
from pathlib import Path

import pytest


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
        self.assertIn('CMD ["msc-convert", "--help"]', content)
        stages = dict(re.findall(r"(?ms)^FROM [^\n]+ AS (\w+)\n(.*?)(?=^FROM |\Z)", content))
        self.assertEqual(list(stages)[-1], "metadata")
        self.assertIn("RUN pip install --no-cache-dir .", stages["base"])
        for excluded in ("h5ad", "openjdk", "nextflow", "gffread", "docker_cli"):
            self.assertNotIn(excluded, stages["base"] + stages["metadata"])
        self.assertIn("FROM base AS full", content)
        self.assertIn("FROM base AS metadata", content)

    def test_dockerfile_installs_pinned_nfcore_runtime(self):
        content = (ROOT / "Dockerfile").read_text(encoding="utf-8")

        self.assertIn("ARG NEXTFLOW_VERSION=26.04.6", content)
        self.assertIn("ARG DOCKER_CLI_VERSION=29.8.1", content)
        self.assertIn("openjdk-21-jre-headless", content)
        self.assertIn("gffread", content)
        self.assertNotIn("openjdk-17-jre-headless", content)
        self.assertIn("182a63c74074e2dc7956ffa3c8cd59de952ed2c44394e21faf5e1736b945444c", content)
        self.assertIn("sha256sum --check --strict", content)
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

        self.assertIn('user: "0:0"', content)
        self.assertNotIn('user: "${RUNNER_UID', content)
        self.assertIn(
            "JSON2H5AD_OUT: ${JSON2H5AD_OUT:?set JSON2H5AD_OUT}",
            content,
        )
        self.assertIn("source: ${JSON2H5AD_OUT:?set JSON2H5AD_OUT}", content)
        self.assertIn("target: ${JSON2H5AD_OUT:?set JSON2H5AD_OUT}", content)
        self.assertIn("target: full", content)
        self.assertIn('command: ["msc-convert", "--help"]', content)
        self.assertNotIn("NEXTFLOW_VERSION:", content)
        self.assertNotIn("DOCKER_CLI_VERSION:", content)
        self.assertIn("read_only: true", content)
        self.assertIn("no-new-privileges:true", content)
        self.assertIn("cap_drop:", content)
        self.assertIn("- ALL", content)
        self.assertIn("NXF_OPTS: -Djava.io.tmpdir=/nextflow-tmp", content)
        self.assertIn("/tmp:size=2g,mode=1777", content)
        self.assertIn("/nextflow-tmp:size=2g,mode=1777,exec", content)

    @pytest.mark.fake_process
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
        self.assertIn(
            '"u:${project_owner_uid}:rwx,u:${runner_name}:rwx"',
            provision,
        )
        self.assertIn(
            '"d:u:${project_owner_uid}:rwx,d:u:${runner_name}:rwx,'
            'd:u::rwx,d:g::---,d:o::---"',
            provision,
        )
        self.assertIn("verify_output_acl", provision)
        self.assertIn("getfacl", provision)

        runner = scripts[1].read_text(encoding="utf-8")
        self.assertIn("verify_output_acl", runner)
        self.assertIn("Mapped output ownership", runner)
        self.assertNotIn("chown", runner)

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


@pytest.mark.fake_process
def test_compose_runner_uses_only_reviewed_fake_path_and_socket_seam(
    tmp_path, monkeypatch
):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()

    def executable(name, content):
        path = fake_bin / name
        path.write_text(f"#!/usr/bin/bash\n{content}\n", encoding="utf-8")
        path.chmod(0o755)

    for name in ("bash", "cut", "dirname", "mkdir", "readlink"):
        os.symlink(f"/usr/bin/{name}", fake_bin / name)
    executable(
        "id",
        'case "${1:-}" in -un) echo nfcore-runner;; -u) echo 12345;; *) exit 2;; esac',
    )
    executable(
        "getent",
        'echo "nfcore-runner:x:12345:12345::/tmp/fake-runner:/bin/false"',
    )
    executable("getfacl", "exit 0")
    executable("find", "exit 0")
    executable("stat", 'echo "12345:12345"')
    executable(
        "docker",
        'echo "$*" >>"${FAKE_DOCKER_LOG}"\n'
        'if [[ "${1:-}" == "info" ]]; then echo \'["name=rootless"]\'; fi',
    )

    socket_path = tmp_path / "rootless.sock"
    rootless = socket.socket(socket.AF_UNIX)
    try:
        rootless.bind(os.fspath(socket_path))
    except PermissionError:
        rootless.close()
        rootless = None
        socket_path = next(
            path
            for path in (Path("/run/uuidd/request"), Path("/run/docker.sock"))
            if path.is_socket()
        )
    output = tmp_path / "output"
    docker_log = tmp_path / "docker.log"
    monkeypatch.setenv("PATH", os.fspath(fake_bin))
    monkeypatch.setenv("ROOTLESS_DOCKER_SOCKET", os.fspath(socket_path))
    monkeypatch.setenv("JSON2H5AD_OUT", os.fspath(output))
    monkeypatch.setenv("FAKE_DOCKER_LOG", os.fspath(docker_log))

    try:
        completed = subprocess.run(
            [os.fspath(ROOT / "scripts" / "json2h5ad-compose.sh"), "run", "--rm"],
            text=True,
            capture_output=True,
            check=False,
        )
    finally:
        if rootless is not None:
            rootless.close()

    assert completed.returncode == 0, completed.stderr
    calls = docker_log.read_text(encoding="utf-8").splitlines()
    assert calls[0] == "info --format {{json .SecurityOptions}}"
    assert calls[1].endswith("compose.yaml run --rm")


if __name__ == "__main__":
    unittest.main()
