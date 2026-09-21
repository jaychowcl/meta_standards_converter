# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""MAGE-TAB field meaning survives the typed MINiML boundary."""
from copy import deepcopy
from collections import Counter
import csv
import json
from pathlib import Path

import pytest

from meta_standards_converter.converters import Converter
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import AEWebFetcher
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.model import Reference

FIXTURES = Path(__file__).parents[1] / "fixtures"


def table(path):
    with path.open(newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def test_array_roundtrip_is_strict_and_preserves_namespace_and_units():
    path = FIXTURES / "fidelity/E-MTAB-16847.idf.txt"
    package = AEParser().parse(AEWebFetcher().resolve(str(path)))
    MINiMLCodec().decode(package, strict=True)
    assert package.platforms[0].technology == "other"
    result = Converter().convert(package, out_type="magetab", enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    rows = next(row[1] for row in result.items[0].payload[0] if row[0] == "SDRF File")
    assert rows == table(path.with_name("E-MTAB-16847.sdrf.txt"))
    ref = next(node.array_design_ref for node in package.series.assay_paths[0].steps
               if getattr(node, "array_design_ref", None))
    assert ref.term_source_ref == "ArrayExpress"


def test_reference_qualifiers_are_typed_and_backward_compatible():
    assert Reference.from_mapping({"ref": "A-GEOD-1"}).to_mapping() == {"ref": "A-GEOD-1"}
    value = {"ref": "A-GEOD-1", "term_source_ref": "ArrayExpress", "term_accession_number": "A-GEOD-1"}
    reference = Reference.from_mapping(value)
    assert reference.term_source_ref == "ArrayExpress"
    assert reference.term_accession_number == "A-GEOD-1"
    assert reference.to_mapping() == value


def test_geo_array_declares_the_namespace_of_its_actual_identifier():
    result = Converter().convert(FIXTURES / "studies/GSE100/inputs/geo.xml",
        out_type="magetab", enrichment="off", options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    idf = result.items[0].payload[0]
    sdrf = next(row[1] for row in idf if row[0] == "SDRF File")
    col = sdrf[0].index("Array Design REF")
    assert sdrf[0][col + 1] == "Term Source REF"
    assert {(row[col], row[col+1]) for row in sdrf[1:]} == {("GPL221", "GEO")}
    assert "GEO" in next(row[1:] for row in idf if row[0] == "Term Source Name")


@pytest.mark.parametrize("reference,namespace", [("GPL221", "GEO"), ("A-GEOD-21185", "ArrayExpress"), ("custom-design", None)])
def test_platform_namespace_import_is_not_assumed(reference, namespace):
    parser = AEParser()
    platforms, sample = {}, {}
    parser._map_platform(["Array Design REF", "Technology Type"], [reference, "array assay"], sample, platforms)
    assert platforms[reference]["accession"][0].get("database") == namespace
    assert platforms[reference]["technology"] == "other"


def valid_tables():
    return [["Protocol Name", "P1"], ["Experimental Factor Name", "dose"],
            ["Term Source Name", "EFO"], ["SDRF File", [
                ["Source Name", "Protocol REF", "Assay Name", "Factor Value[dose]", "Unit", "Term Source REF"],
                ["s1", "P1", "a1", "0", "nanomolar", "EFO"]]]]


@pytest.mark.parametrize("defect", ["width", "protocol", "factor", "source", "qualifier"])
def test_output_validation_rejects_broken_references_and_qualifiers(defect):
    from meta_standards_converter.magetab.validation import validate_magetab, MAGETabValidationError
    rows = valid_tables()
    validate_magetab(rows)
    sdrf = rows[-1][1]
    if defect == "width": sdrf[1].pop()
    if defect == "protocol": sdrf[1][1] = "missing"
    if defect == "factor": rows[1][1] = "other factor"
    if defect == "source": sdrf[1][-1] = "undeclared"
    if defect == "qualifier": sdrf[0][-2] = "Comment[unrelated]"
    with pytest.raises(MAGETabValidationError):
        validate_magetab(rows)


def test_writer_rejects_invalid_bundle_before_creating_files(tmp_path):
    from meta_standards_converter.magetab.writer import MAGETabWriter
    rows = valid_tables()
    rows[-1][1][1][1] = "missing"
    with pytest.raises(ValueError, match="protocol"):
        MAGETabWriter().write(rows, str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()


PANEL = ["E-MTAB-13662", "E-MTAB-14560", "E-MTAB-14566", "E-MTAB-16192",
         "E-MTAB-16253", "E-MTAB-13131", "E-MTAB-16856", "E-MTAB-14416",
         "E-MTAB-16847", "E-MTAB-14960"]


@pytest.mark.parametrize("accession", PANEL)
def test_recorded_archive_panel_preserves_ordered_fields_and_rows(accession):
    path = FIXTURES / "fidelity" / (accession + ".idf.txt")
    result = Converter().convert(path, out_type="magetab", enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    actual = next(r[1] for r in result.items[0].payload[0] if r[0] == "SDRF File")
    expected = table(path.with_name(accession + ".sdrf.txt"))
    # Empty optional file nodes have no scientific identity to preserve.
    omit = {i for i, h in enumerate(expected[0]) if h == "Derived Array Data File"
            and all(not row[i].strip() for row in expected[1:])}
    expected = [[v for i, v in enumerate(row) if i not in omit] for row in expected]
    normalize = lambda row: tuple("" if v is None else str(v).strip() for v in row)
    def columns(headers):
        counts = Counter()
        keys = []
        for h in headers:
            counts[h] += 1
            keys.append((h, counts[h]))
        return keys
    akeys, ekeys = columns(actual[0]), columns(expected[0])
    assert set(akeys) == set(ekeys)
    positions = [akeys.index(key) for key in ekeys]
    assert Counter(normalize([r[i] for i in positions]) for r in actual[1:]) == Counter(map(normalize, expected[1:]))
    # Attribute ordering within a node may normalize; graph node/edge order may not.
    graph = lambda headers: [h for h in headers if h.endswith(" Name") or h.endswith(" File") or h == "Protocol REF"]
    assert graph(actual[0]) == graph(expected[0])


@pytest.mark.parametrize("target", ["json", "magetab", "tsv", "csv"])
def test_recorded_erp000263_package_converts_without_enrichment(target):
    value = json.loads((FIXTURES / "fidelity/ERP000263.ena.json").read_text())
    result = Converter().convert(value, out_type=target, enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()


def test_node_comments_do_not_migrate_onto_a_characteristic_or_factor():
    from meta_standards_converter.magetab.semantics import _miniml_path_columns
    from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
    import io
    steps = [{"kind": "source", "name": "sample", "comments": [{"name": "ENA_SAMPLE", "value": "ERS1"}],
              "characteristics": [{"name": "organism", "value": "Homo sapiens"}]},
             {"kind": "scan", "name": "run", "comments": [{"name": "ENA_RUN", "value": "ERR1"}],
              "factor_values": [{"name": "dose", "value": "0"}]}]
    cells = _miniml_path_columns(steps)
    sdrf = io.StringIO()
    csv.writer(sdrf, delimiter="\t").writerows([[k for k,v in cells], [v for k,v in cells]])
    source = MAGETabInput(TextResource("test.idf.txt", "Investigation Accession\tE-MTAB-1\nExperimental Factor Name\tdose\n", "test"),
                         (TextResource("test.sdrf.txt", sdrf.getvalue(), "test"),), "test", "path")
    parsed = AEParser().parse(source).to_mapping()["series"]["assay_paths"][0]["steps"]
    for expected, actual in zip(steps, parsed):
        assert expected["comments"] == actual["comments"]


def test_native_sra_broker_matches_available_library_facts_without_enrichment():
    value = json.loads((FIXTURES / 'fidelity/ERP185509.sra.json').read_text())[0]
    assert value['series']['iid'] == 'ERP185509'
    result = Converter().convert(value, out_type='magetab', enrichment='off',
                                 options={'expand_studies': False})
    assert result.status == 'complete', result.to_dict()
    native = next(r[1] for r in result.items[0].payload[0] if r[0] == 'SDRF File')
    broker = table(FIXTURES / 'fidelity/E-MTAB-16253.sdrf.txt')
    fields = ['Comment[ENA_RUN]', 'Comment[LIBRARY_LAYOUT]', 'Comment[LIBRARY_STRATEGY]',
              'Comment[LIBRARY_SOURCE]', 'Comment[LIBRARY_SELECTION]']
    def facts(rows):
        indices = [rows[0].index(h) for h in fields]
        return {tuple(row[i] for i in indices) for row in rows[1:]}
    assert len(facts(native)) == 8
    assert facts(native) == facts(broker)
    expected_runs = {row[0] for row in facts(broker)}
    assert {r['run'] for s in value['sample'] for r in s.get('sra_run', [])} == expected_runs
    # These factors are authored in AE but are not explicit per-sample facts in
    # the captured native package. Shared treatment prose cannot assign doses.
    assert 'Factor Value[compound]' in broker[0] and 'Factor Value[dose]' in broker[0]
    assert 'Factor Value[dose]' not in native[0]
    assert all(record['status'] == 'skipped' for record in result.items[0].preparation)


def test_repeated_unknown_idf_comments_and_scoped_dates_survive_export():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
    source = MAGETabInput(TextResource('test.idf.txt',
        'Investigation Accession\tE-MTAB-1\nPublic Release Date\t2025-01-01\n'
        'Comment [Custom note]\t0\tfirst\nComment [Custom note]\tsecond\n'
        'Comment[ArrayExpressReleaseDate]\t2025-02-01\t2025-03-01\n', 'test'),
        (TextResource('test.sdrf.txt', 'Source Name\nsample\n', 'test'),), 'test', 'path')
    package = AEParser().parse(source)
    rows = AEConstructor().miniml2magetab(package)
    assert next(r[1:] for r in rows if r[0] == 'Comment[Custom note]') == ['0', 'first', 'second']
    assert next(r[1:] for r in rows if r[0] == 'Comment[ArrayExpressReleaseDate]') == ['2025-02-01', '2025-03-01']
    assert next(r[1:] for r in rows if r[0] == 'Public Release Date') == ['2025-01-01']
