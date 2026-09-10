# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
import hashlib
import logging
import os
import shutil
import subprocess
from pathlib import Path
from urllib.parse import urlparse
from meta_standards_converter.expression.planning import SourcePlanner

logger = logging.getLogger(__name__)

class ReferenceResolver:
    """Resolve explicit or safely confirmed nf-core reference parameters."""

    GENOME_BY_TAXID = {"9606": "GRCh38", "10090": "GRCm39"}
    GENOME_BY_NAME = {"homo sapiens": "GRCh38", "mus musculus": "GRCm39"}

    def resolve(
        self,
        packages: list[dict],
        genome: str | None = None,
        fasta: str | None = None,
        gtf: str | None = None,
        gff: str | None = None,
        accept_inferred: bool = False,
    ) -> dict:
        if gtf and gff:
            raise ValueError("Reference annotation options gtf and gff are mutually exclusive.")
        annotation = gtf or gff
        if genome:
            return {
                "genome": genome,
                **({"gtf": gtf} if gtf else {}),
                **({"gff": gff} if gff else {}),
            }
        if fasta or annotation:
            if annotation and not fasta:
                raise ValueError("A custom annotation requires either genome or fasta.")
            if fasta and not annotation:
                raise ValueError("A custom fasta reference requires an annotation in gtf or gff format.")
            return {
                "fasta": fasta,
                **({"gtf": gtf} if gtf else {}),
                **({"gff": gff} if gff else {}),
            }

        taxids = set()
        names = set()
        planner = SourcePlanner()
        for package in packages:
            for sample in planner._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                for channel in planner._as_list(sample.get("channel")):
                    if not isinstance(channel, dict):
                        continue
                    for organism in planner._as_list(channel.get("organism")):
                        if isinstance(organism, dict):
                            if organism.get("taxid"):
                                taxids.add(str(organism["taxid"]))
                            if organism.get("value") or organism.get("name"):
                                names.add(str(organism.get("value") or organism.get("name")).lower())
                        elif organism:
                            names.add(str(organism).lower())
        candidates = {
            self.GENOME_BY_TAXID[value]
            for value in taxids
            if value in self.GENOME_BY_TAXID
        } | {
            self.GENOME_BY_NAME[value]
            for value in names
            if value in self.GENOME_BY_NAME
        }
        if len(candidates) != 1:
            raise ValueError(
                "Unable to infer one supported genome reference; provide --genome or --fasta and --gtf."
            )
        inferred = next(iter(candidates))
        if not accept_inferred:
            raise ValueError(
                f"Inferred reference {inferred}; rerun with --accept-inferred-reference "
                "or provide an explicit reference."
            )
        return {"genome": inferred, "inferred": True}


class AnnotationConverter:
    """Validate local annotations and normalize GFF3 input to GTF."""

    def __init__(self, command_runner=None, which=None):
        self.command_runner = command_runner or subprocess.run
        self.which = which or shutil.which

    def prepare(self, reference: dict, reference_dir: Path) -> dict:
        prepared = dict(reference)
        if prepared.get("fasta"):
            prepared["fasta"] = str(self._local_file(prepared["fasta"], "FASTA"))

        annotation = prepared.get("gtf") or prepared.get("gff")
        if not annotation:
            return prepared
        source = self._local_file(annotation, "annotation")
        annotation_format = self._annotation_format(source)
        digest = self._sha256(source)
        metadata = {
            "annotation_source": str(source),
            "annotation_format": annotation_format,
            "annotation_sha256": digest,
        }

        if annotation_format == "gtf":
            prepared["gtf"] = str(source)
            prepared["effective_annotation"] = str(source)
            return {**prepared, **metadata}

        reference_dir.mkdir(parents=True, exist_ok=True)
        destination = reference_dir / f"{digest}.gtf"
        temporary = Path(str(destination) + ".tmp")
        if not destination.exists() or destination.stat().st_size == 0:
            if not self.which("gffread"):
                raise RuntimeError("GFF3 annotations require gffread on PATH.")
            command = ["gffread", str(source), "-T", "-o", str(temporary)]
            completed = self.command_runner(
                command,
                text=True,
                capture_output=True,
                check=False,
            )
            if completed.returncode or not temporary.exists() or temporary.stat().st_size == 0:
                temporary.unlink(missing_ok=True)
                detail = (completed.stderr or completed.stdout or "conversion produced no GTF").strip()
                raise RuntimeError(f"gffread annotation conversion failed: {detail}")
            os.replace(temporary, destination)

        prepared.pop("gff", None)
        prepared["gtf"] = str(destination)
        prepared["effective_annotation"] = str(destination)
        return {**prepared, **metadata}

    def _local_file(self, value: str, label: str) -> Path:
        parsed = urlparse(str(value))
        if parsed.scheme not in ("", "file"):
            raise ValueError(f"User-supplied {label} must be a local file: {value}")
        path = Path(parsed.path if parsed.scheme == "file" else value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"User-supplied {label} file not found: {path}")
        return path

    def _annotation_format(self, path: Path) -> str:
        name = path.name.lower()
        if name.endswith((".gtf", ".gtf.gz")):
            return "gtf"
        if name.endswith((".gff", ".gff.gz", ".gff3", ".gff3.gz")):
            return "gff3"
        raise ValueError(f"Unsupported annotation format: {path}")

    def _sha256(self, path: Path) -> str:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
