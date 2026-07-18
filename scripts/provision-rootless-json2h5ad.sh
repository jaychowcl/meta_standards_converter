#!/usr/bin/env bash
set -euo pipefail

runner_name="nfcore-runner"
project_root="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
output_root="${2:-${project_root}/.out/json2h5ad}"
subid_count=65536

if [[ "$(id -u)" -ne 0 ]]; then
    echo "This provisioning command must run as root." >&2
    exit 1
fi
if [[ ! -f "${project_root}/compose.yaml" ]]; then
    echo "Project root does not contain compose.yaml: ${project_root}" >&2
    exit 1
fi

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
    acl \
    dbus-user-session \
    fuse-overlayfs \
    slirp4netns \
    uidmap

if ! id "${runner_name}" >/dev/null 2>&1; then
    useradd --create-home --shell /bin/bash "${runner_name}"
    passwd --lock "${runner_name}"
fi

range_available() {
    local start="$1"
    awk -F: -v candidate_start="${start}" -v candidate_count="${subid_count}" '
        BEGIN {
            candidate_end = candidate_start + candidate_count - 1
            available = 1
        }
        NF >= 3 {
            existing_start = $2
            existing_end = $2 + $3 - 1
            if (candidate_start <= existing_end && existing_start <= candidate_end) {
                available = 0
            }
        }
        END { exit available ? 0 : 1 }
    ' /etc/subuid /etc/subgid
}

has_subuid=false
has_subgid=false
if grep --quiet "^${runner_name}:" /etc/subuid; then
    has_subuid=true
fi
if grep --quiet "^${runner_name}:" /etc/subgid; then
    has_subgid=true
fi
if [[ "${has_subuid}" != "${has_subgid}" ]]; then
    echo "${runner_name} must have matching entries in /etc/subuid and /etc/subgid." >&2
    exit 1
fi
if [[ "${has_subuid}" == "false" ]]; then
    candidate=100000
    while ! range_available "${candidate}"; do
        candidate=$((candidate + subid_count))
    done
    usermod \
        --add-subuids "${candidate}-$((candidate + subid_count - 1))" \
        --add-subgids "${candidate}-$((candidate + subid_count - 1))" \
        "${runner_name}"
fi

runner_uid="$(id -u "${runner_name}")"
runner_gid="$(id -g "${runner_name}")"
runner_home="$(getent passwd "${runner_name}" | cut -d: -f6)"
runtime_dir="/run/user/${runner_uid}"
project_owner_uid="$(stat --format '%u' "${project_root}")"

mkdir -p "${output_root}"
path="$(dirname "${project_root}")"
while [[ "${path}" != "/" ]]; do
    setfacl --modify "u:${runner_name}:--x" "${path}"
    path="$(dirname "${path}")"
done
path="$(dirname "${output_root}")"
while [[ "${path}" != "${project_root}" && "${path}" != "/" ]]; do
    setfacl --modify "u:${runner_name}:--x" "${path}"
    path="$(dirname "${path}")"
done
setfacl --modify "u:${runner_name}:r-x" "${project_root}"
setfacl --modify "u:${runner_name}:r--" \
    "${project_root}/.dockerignore" \
    "${project_root}/Dockerfile" \
    "${project_root}/LICENSE" \
    "${project_root}/README.md" \
    "${project_root}/compose.yaml" \
    "${project_root}/pyproject.toml"
setfacl --recursive --modify "u:${runner_name}:r-X" "${project_root}/src"
setfacl --recursive --modify \
    "u:${project_owner_uid}:rwx" \
    "u:${runner_name}:rwx" \
    "${output_root}"
setfacl --modify \
    "d:u:${project_owner_uid}:rwx" \
    "d:u:${runner_name}:rwx" \
    "d:u::rwx" \
    "d:g::---" \
    "d:o::---" \
    "${output_root}"

loginctl enable-linger "${runner_name}"
systemctl start "user@${runner_uid}.service"
install --directory --owner "${runner_name}" --group "${runner_gid}" --mode 0700 "${runtime_dir}"

runner_env=(
    "HOME=${runner_home}"
    "XDG_RUNTIME_DIR=${runtime_dir}"
    "DBUS_SESSION_BUS_ADDRESS=unix:path=${runtime_dir}/bus"
)
runuser --user "${runner_name}" -- env "${runner_env[@]}" \
    dockerd-rootless-setuptool.sh install --force
runuser --user "${runner_name}" -- env "${runner_env[@]}" \
    systemctl --user enable --now docker.service

security_options="$(
    runuser --user "${runner_name}" -- env \
        "DOCKER_HOST=unix://${runtime_dir}/docker.sock" \
        docker info --format '{{json .SecurityOptions}}'
)"
if [[ "${security_options,,}" != *rootless* ]]; then
    echo "Rootless Docker verification failed: ${security_options}" >&2
    exit 1
fi

echo "Rootless json2h5ad runner is ready."
echo "Runner UID:GID: ${runner_uid}:${runner_gid}"
echo "Output: ${output_root}"
echo "SecurityOptions: ${security_options}"
echo "Run: sudo -u ${runner_name} -H ${project_root}/scripts/json2h5ad-compose.sh build converter"
