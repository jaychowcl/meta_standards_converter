# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit route settings; unknown options are never forwarded silently."""

from collections.abc import Mapping

TARGETS = frozenset({"json", "magetab", "tsv", "csv", "h5ad", "obs"})
EXPRESSION = frozenset(
    {
        "explicit_assets",
        "asset_manifest",
        "asset_specs",
        "force_reprocess",
        "matrix_orientation",
        "pipeline",
        "genome",
        "fasta",
        "gtf",
        "gff",
        "accept_inferred_reference",
        "profile",
        "revision",
        "params_file",
        "nextflow_config",
        "work_dir",
        "resume",
        "force_memory",
        "processed_checkpoint_dir",
        "allow_invalid",
        "allow_unverified_combination",
        "replacement_profile",
    }
)
OUTPUT = {
    "json": {"aggregate"},
    "magetab": {"enrich", "platform_handler", "replacement_profile"},
    "tsv": {"aggregate", "allow_invalid", "replacement_profile"},
    "csv": {"aggregate", "allow_invalid", "replacement_profile"},
    "h5ad": set(EXPRESSION),
    "obs": set(EXPRESSION) | {"include_var", "include_uns"},
}
INPUT = {
    "sra_records": set(),
    "ena_records": set(),
    "geo_accession": {"enrich", "related_series", "remove_empty"},
    "sra_accession": {"enrich_from_geo_ae", "include_peer", "evidence_dir"},
    "ena_accession": {"enrich_from_geo_ae", "include_peer", "evidence_dir"},
    "ae_accession": {"sdrf_sources"},
    "magetab": {"sdrf_sources"},
    "geo_xml": {"remove_empty"},
    "geo_archive": {"remove_empty"},
    "sra_xml": set(),
    "ena_xml": set(),
    "json": set(),
    "miniml": set(),
    "atlas": set(),
    "curator": set(),
    "h5ad": set(),
    "anndata": set(),
    "matrix": {"orientation"},
    "fastq": set(),
}
RUNTIME_DEFAULTS = {
    "overwrite": False,
    "fail_fast": False,
    "recursive": False,
    "allow_processing": False,
    "resource_profile": "standard",
    "resource_overrides": None,
    "insdc_default": "ena",
    "allowed_hosts": (),
}
BOOLS = {
    "overwrite",
    "fail_fast",
    "recursive",
    "allow_processing",
    "enrich",
    "related_series",
    "remove_empty",
    "enrich_from_geo_ae",
    "include_peer",
    "aggregate",
    "allow_invalid",
    "force_reprocess",
    "accept_inferred_reference",
    "resume",
    "force_memory",
    "allow_unverified_combination",
    "include_var",
    "include_uns",
}


def settings(value, allowed, label):
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    unknown = set(value) - set(allowed)
    if unknown:
        raise ValueError(f'Unsupported {label}: {", ".join(sorted(map(str, unknown)))}')
    result = dict(value)
    for key, item in result.items():
        if key in BOOLS and not isinstance(item, bool):
            raise TypeError(f"{key} must be a boolean")
    for key in ("orientation", "matrix_orientation"):
        if key in result and result[key] not in {
            "auto",
            "genes-by-observations",
            "observations-by-genes",
        }:
            raise ValueError(f"Unsupported {key}")
    if "pipeline" in result and result["pipeline"] not in {
        "auto",
        "rnaseq",
        "scrnaseq",
    }:
        raise ValueError("Unsupported pipeline")
    if result.get("gtf") and result.get("gff"):
        raise ValueError("gtf and gff are mutually exclusive")
    if result.get("force_memory") and not result.get("resume"):
        raise ValueError("force_memory requires resume=True")
    if "insdc_default" in result and result["insdc_default"] not in {"ena", "sra"}:
        raise ValueError("insdc_default must be ena or sra")
    return result
