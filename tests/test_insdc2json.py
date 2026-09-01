# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from meta_standards_converter.converters.ena2json import ena2json
from meta_standards_converter.converters.insdc2json import (
    StudyConversionFailures,
    merge_origin_package,
)
from meta_standards_converter.converters.sra2json import sra2json
from meta_standards_converter.insdc_handlers.study_models import ProviderDocument, StudyFetchResult
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser
from meta_standards_converter.miniml import MINiMLCodec


ROOT = Path(__file__).resolve().parents[1]


class Fetcher:
    def __init__(self, results):
        self.results = results

    def fetch(self, accession):
        return self.results


def fixture_result(provider="sra"):
    if provider == "sra":
        path = ROOT / "docs" / "sra" / "fixtures" / "SRX017289" / "sra-efetch.xml"
        return StudyFetchResult(
            provider="sra",
            requested_accession="SRX017289",
            study_accession="SRP002056",
            documents=(ProviderDocument("sra_efetch", path.name, "https://example.test/sra", "application/xml", path.read_bytes()),),
        )
    fixture = ROOT / "docs" / "ena" / "fixtures" / "SRX017289"
    return StudyFetchResult(
        provider="ena",
        requested_accession="SRX017289",
        study_accession="SRP002056",
        documents=tuple(
            ProviderDocument(kind, name, f"https://example.test/{name}", "application/xml" if name.endswith(".xml") else "application/json", (fixture / name).read_bytes())
            for kind, name in (
                ("ena_study", "study.xml"),
                ("ena_sample", "sample.xml"),
                ("ena_experiment", "experiment.xml"),
                ("ena_run", "run.xml"),
                ("ena_file_report", "file-report.json"),
            )
        ),
    )


@pytest.mark.parametrize(
    ("converter_type", "provider", "suffix"),
    [(sra2json, "sra", ".sra.json"), (ena2json, "ena", ".ena.json")],
)
def test_converter_returns_typed_packages_and_writes_one_list_per_study(tmp_path, converter_type, provider, suffix) -> None:
    converter = converter_type(fetcher=Fetcher([fixture_result(provider)]), parser=INSDCStudyParser())

    packages = converter.convert("SRX017289", out=tmp_path)

    assert [package.series.iid for package in packages] == ["SRP002056"]
    destination = tmp_path / f"SRP002056{suffix}"
    assert destination.exists()
    payload = json.loads(destination.read_text())
    assert isinstance(payload, list) and payload[0]["series"]["iid"] == "SRP002056"


def test_converter_rejects_mutually_exclusive_enrichment() -> None:
    converter = sra2json(fetcher=Fetcher([fixture_result()]), parser=INSDCStudyParser())
    with pytest.raises(ValueError, match="mutually exclusive"):
        converter.convert("SRX017289", enrich_geo=True, enrich_ae=True)


def test_requested_enrichment_requires_compatible_link() -> None:
    path = ROOT / "docs" / "sra" / "fixtures" / "SRX7812918" / "sra-efetch.xml"
    result = StudyFetchResult(
        provider="sra",
        requested_accession="SRX7812918",
        study_accession="SRP250911",
        documents=(ProviderDocument("sra_efetch", path.name, None, "application/xml", path.read_bytes()),),
    )
    converter = sra2json(fetcher=Fetcher([result]), parser=INSDCStudyParser())

    with pytest.raises(ValueError, match="has no compatible GEO link"):
        converter.convert("SRX7812918", enrich_geo=True)


def test_multi_study_results_are_stably_sorted() -> None:
    first = fixture_result()
    second = StudyFetchResult(
        provider=first.provider,
        requested_accession=first.requested_accession,
        study_accession="ARP000001",
        documents=first.documents,
    )
    class Parser:
        def parse(self, result, *, enrichment="none"):
            package = INSDCStudyParser().parse(first, enrichment=enrichment)
            data = json.loads(json.dumps(package.to_mapping()))
            data["series"]["iid"] = result.study_accession
            data["series"]["accession"] = [{"value": result.study_accession}]
            return MINiMLCodec().decode(data).package

    converter = sra2json(fetcher=Fetcher([first, second]), parser=Parser())

    packages = converter.convert("PRJNA1")

    assert [package.series.iid for package in packages] == ["ARP000001", "SRP002056"]


def test_duplicate_resolved_studies_are_converted_once() -> None:
    result = fixture_result()
    parser = Mock(wraps=INSDCStudyParser())
    converter = sra2json(fetcher=Fetcher([result, result]), parser=parser)

    packages = converter.convert("SRX017289")

    assert len(packages) == 1
    assert parser.parse.call_count == 1


def test_project_continues_after_one_resolved_study_fails(tmp_path) -> None:
    good = fixture_result()
    bad = StudyFetchResult(
        good.provider,
        good.requested_accession,
        "ARP000001",
        good.documents,
    )

    class Parser:
        def parse(self, result, *, enrichment="none"):
            if result.study_accession == "ARP000001":
                raise ValueError("broken study")
            return INSDCStudyParser().parse(result, enrichment=enrichment)

    converter = sra2json(fetcher=Fetcher([good, bad]), parser=Parser())

    with pytest.raises(StudyConversionFailures, match="ARP000001"):
        converter.convert("PRJNA1", out=tmp_path)

    assert (tmp_path / "SRP002056.sra.json").exists()


def test_explicit_geo_enrichment_uses_origin_suffix_and_precedence(tmp_path) -> None:
    result = fixture_result()
    base = INSDCStudyParser().parse(result)
    origin_data = MINiMLCodec().encode(base)
    origin_data["series"]["title"] = "GEO submitted title"
    origin = MINiMLCodec().decode(origin_data).package
    geo = Mock()
    geo.convert.return_value = [origin]
    converter = sra2json(
        fetcher=Fetcher([result]),
        parser=INSDCStudyParser(),
        geo_converter=geo,
    )

    packages = converter.convert("SRX017289", enrich_geo=True, out=tmp_path)

    assert packages[0].series.title == "GEO submitted title"
    assert (tmp_path / "SRP002056.sra.geo.json").exists()
    geo.convert.assert_called_once_with(gse="GSE18729", enrich=False)


def test_ambiguous_origin_sample_alignment_is_preserved_without_guessing() -> None:
    package = INSDCStudyParser().parse(fixture_result())
    base_data = MINiMLCodec().encode(package)
    duplicate = json.loads(json.dumps(base_data["sample"][0]))
    duplicate["iid"] = "SRS999999"
    duplicate["title"] = "second INSDC sample"
    base_data["sample"].append(duplicate)
    base = MINiMLCodec().decode(base_data).package

    origin_data = MINiMLCodec().encode(package)
    origin_sample = origin_data["sample"][0]
    origin_sample["iid"] = "GSM465245"
    origin_sample["title"] = "GEO biology title"
    origin_sample["accession"] = [{"value": "GSM465245", "database": "GEO"}]
    origin_sample.pop("sra_run", None)
    origin = MINiMLCodec().decode(origin_data).package

    merged = MINiMLCodec().encode(merge_origin_package(base, origin, source="geo"))

    assert [sample["title"] for sample in merged["sample"]] == [
        "small RNAs_wild type (fog-2)_adult males_part1",
        "second INSDC sample",
        "GEO biology title",
    ]
    codes = [
        item["code"]
        for item in merged["extensions"]["insdc"]["enrichment"]["diagnostics"]
    ]
    assert codes.count("ambiguous_sample_alignment") == 2
    assert "unmatched_origin_sample" in codes
