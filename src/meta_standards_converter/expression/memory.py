# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
import logging
import os
from pathlib import Path


logger = logging.getLogger(__name__)

def _available_memory_bytes() -> int:
    """Return the conservative host/cgroup memory currently available."""

    candidates: list[int] = []
    try:
        pages = int(os.sysconf("SC_AVPHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        if pages > 0 and page_size > 0:
            candidates.append(pages * page_size)
    except (AttributeError, OSError, TypeError, ValueError):
        pass

    for limit_path, current_path in (
        (Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")),
        (
            Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
            Path("/sys/fs/cgroup/memory/memory.usage_in_bytes"),
        ),
    ):
        try:
            limit_text = limit_path.read_text(encoding="ascii").strip()
            if limit_text == "max":
                continue
            remaining = int(limit_text) - int(
                current_path.read_text(encoding="ascii").strip()
            )
            if remaining > 0:
                candidates.append(remaining)
        except (OSError, TypeError, ValueError):
            continue

    if not candidates:
        raise RuntimeError("Unable to determine available memory safely.")
    return min(candidates)


def _path_size_bytes(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    if path.is_dir():
        total = 0
        for member in path.rglob("*"):
            if member.is_file():
                total += member.stat().st_size
        return total
    return 0


def _estimate_asset_memory_bytes(path: str, asset: "Asset") -> int:
    """Estimate peak resident bytes before loading an expression asset."""

    local = Path(path)
    on_disk = max(1, _path_size_bytes(local))
    name = local.name.casefold()
    for compression_suffix in (".gz", ".bz2", ".xz", ".zip"):
        if name.endswith(compression_suffix):
            name = name[: -len(compression_suffix)]
            break
    suffix = Path(name).suffix
    multiplier = 6
    if str(path).casefold().endswith((".gz", ".bz2", ".xz", ".zip")):
        multiplier = 12
    estimate = on_disk * multiplier

    if asset.kind == "h5ad" or suffix == ".h5":
        try:
            import h5py

            logical_bytes = 0
            with h5py.File(local, "r") as handle:
                def account(_name, item):
                    nonlocal logical_bytes
                    if isinstance(item, h5py.Dataset):
                        logical_bytes += int(item.size) * int(item.dtype.itemsize)

                handle.visititems(account)
            estimate = max(estimate, logical_bytes * 3)
        except (ImportError, OSError, TypeError, ValueError):
            pass
    return max(1, int(estimate))
