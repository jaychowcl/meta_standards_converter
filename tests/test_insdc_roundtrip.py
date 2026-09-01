# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from pathlib import Path

import pytest

from meta_standards_converter.converters.ae2json import ae2json
from meta_standards_converter.converters.ena2json import ena2json
from meta_standards_converter.converters.json2ae import json2ae
from meta_standards_converter.converters.sra2json import sra2json
from meta_standards_converter.insdc_handlers.study_models import ProviderDocument, StudyFetchResult
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser


ROOT = Path(__file__).resolve().parents[1]


class Fetcher:
    def __init__(self, result):
        self.result = result

    def fetch(self, accession):
        return [self.result]


def _document(path: Path, kind: str) -> ProviderDocument:
    return ProviderDocument(
        kind,
        path.name,
        None,
        "application/json" if path.suffix == ".json" else "application/xml",
        path.read_bytes(),
    )


@pytest.mark.parametrize("provider", ["sra", "ena"])
def test_insdc_miniml_magetab_miniml_roundtrip_without_render_enrichment(
    tmp_path, provider
) -> None:
    if provider == "sra":
        fixture = ROOT / "docs" / "sra" / "fixtures" / "SRX017289"
        result = StudyFetchResult(
            "sra",
            "SRX017289",
            "SRP002056",
            (
                _document(fixture / "sra-efetch.xml", "sra_efetch"),
                _document(fixture / "pubmed-esummary.xml", "pubmed_esummary"),
            ),
        )
        converter = sra2json(fetcher=Fetcher(result), parser=INSDCStudyParser())
    else:
        fixture = ROOT / "docs" / "ena" / "fixtures" / "SRX017289"
        result = StudyFetchResult(
            "ena",
            "SRX017289",
            "SRP002056",
            (
                *tuple(
                _document(fixture / name, kind)
                for kind, name in (
                    ("ena_study", "study.xml"),
                    ("ena_sample", "sample.xml"),
                    ("ena_experiment", "experiment.xml"),
                    ("ena_run", "run.xml"),
                    ("ena_file_report", "file-report.json"),
                )
                ),
                _document(
                    ROOT / "docs" / "sra" / "fixtures" / "SRX017289" / "pubmed-esummary.xml",
                    "pubmed_esummary",
                ),
            ),
        )
        converter = ena2json(fetcher=Fetcher(result), parser=INSDCStudyParser())

    json_dir = tmp_path / "json"
    [package] = converter.convert("SRX017289", out=json_dir)
    source = next(json_dir.glob("*.json"))
    ae_dir = tmp_path / "ae"
    magetabs = json2ae().convert(str(source), out=str(ae_dir), enrich=False)
    idf = ae_dir / "SRP002056.idf.txt"
    [roundtripped] = ae2json().convert(str(idf))

    assert magetabs
    assert idf.exists()
    assert (ae_dir / "SRP002056.sdrf.txt").exists()
    assert roundtripped.series.iid == "SRP002056"
    assert roundtripped.series.title == package.series.title
