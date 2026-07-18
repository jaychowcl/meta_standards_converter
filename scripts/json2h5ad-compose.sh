#!/usr/bin/env bash
# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
set -euo pipefail

runner_name="nfcore-runner"
if [[ "$(id -un)" != "${runner_name}" ]]; then
    echo "Run this command as ${runner_name}, for example:" >&2
    echo "  sudo -u ${runner_name} -H $0 $*" >&2
    exit 1
fi

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
runner_uid="$(id -u)"
runner_gid="$(id -g)"
runner_home="$(getent passwd "${runner_name}" | cut -d: -f6)"
runtime_dir="/run/user/${runner_uid}"
rootless_socket="${runtime_dir}/docker.sock"
output_root="${JSON2H5AD_OUT:-${project_root}/.out/json2h5ad}"

if [[ ! -S "${rootless_socket}" ]]; then
    echo "Rootless Docker socket is unavailable: ${rootless_socket}" >&2
    exit 1
fi

mkdir -p "${output_root}/home" "${output_root}/.nextflow"

export DOCKER_HOST="unix://${rootless_socket}"
export HOME="${runner_home}"
export JSON2H5AD_OUT="$(readlink -f "${output_root}")"
export ROOTLESS_DOCKER_SOCKET="${rootless_socket}"
export RUNNER_GID="${runner_gid}"
export RUNNER_UID="${runner_uid}"
export XDG_RUNTIME_DIR="${runtime_dir}"

security_options="$(docker info --format '{{json .SecurityOptions}}')"
if [[ "${security_options,,}" != *rootless* ]]; then
    echo "Refusing to run against a non-rootless Docker daemon." >&2
    exit 1
fi

exec docker compose --file "${project_root}/compose.yaml" "$@"
