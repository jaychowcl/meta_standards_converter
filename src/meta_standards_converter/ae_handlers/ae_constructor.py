# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Constructor class for ae MAGETAB idf and sdrf
'''
from meta_standards_converter.ae_handlers.ae_idf_handlers import IDFConstructor
from meta_standards_converter.ae_handlers.ae_model import overlay_core, render_model
from meta_standards_converter.ae_handlers.ae_roundtrip import restore_extensions, semantic_sha256, unchanged_magetab
from meta_standards_converter.helpers.json_helper import JSONHandler
from meta_standards_converter.harmonizers.harmonizers import Harmonizer

import copy
import csv
import os
import re


class ProtocolRegistry:
    LABEL_BY_KIND = {
        "manufacture": "Manufacture-Protocol",
        "treatment": "Treatment-Protocol",
        "growth": "Growth-Protocol",
        "extraction": "Extract-Protocol",
        "extract": "Extract-Protocol",
        "library construction": "Library-Construction-Protocol",
        "labeling": "Label-Protocol",
        "label": "Label-Protocol",
        "hybridization": "Hybridization-Protocol",
        "scan": "Scan-Protocol",
        "data processing": "Data-Processing",
        "sample collection": "Sample-Collection-Protocol",
        "nucleic acid sequencing": "Nucleic-Acid-Sequencing-Protocol",
    }

    def __init__(self, series_accession: str):
        self.series_accession = series_accession
        self.by_key = {}

    def get_ref(self, kind: str, text: str | None, label: str | None = None) -> str | None:
        text = self.clean(text)
        if not text:
            return None

        label = label or self.LABEL_BY_KIND.get(kind, kind)
        key = (kind, text)
        if key not in self.by_key:
            self.by_key[key] = {
                "ref": f"P-{self.series_accession}-{len(self.by_key) + 1}",
                "kind": kind,
                "label": label,
                "text": text,
            }
        return self.by_key[key]["ref"]

    def ensure_required(self, kind: str, label: str | None = None) -> str:
        label = label or self.LABEL_BY_KIND.get(kind, kind)
        required_type = Harmonizer().geoprotocols2efo(protocol_type=label)[0]
        for record in self.records():
            record_type = Harmonizer().geoprotocols2efo(protocol_type=record["label"])[0]
            if record_type == required_type:
                return record["ref"]

        key = (kind, "")
        if key not in self.by_key:
            self.by_key[key] = {
                "ref": f"P-{self.series_accession}-{len(self.by_key) + 1}",
                "kind": kind,
                "label": label,
                "text": "",
                "required": True,
            }
        return self.by_key[key]["ref"]

    def records(self) -> list[dict]:
        return list(self.by_key.values())

    def clean(self, value):
        if value is None:
            return None
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split())


from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor, normalized_extension


class AEConstructor:
    def __init__(self, idf_constructor=None, sdrf_constructor=None):
        self.idf_constructor = idf_constructor or IDFConstructor()
        self.sdrf_constructor = sdrf_constructor or SDRFConstructor()

    def miniml2magetab(self, data: dict)-> list:
        """
        converts miniml json to magetab idf. Walks through sections of idf to extract from miniml
        """
        preserved = unchanged_magetab(data)
        if preserved is not None:
            return preserved
        mage_tab = data.get("mage_tab") if isinstance(data, dict) else None
        model = mage_tab.get("model") if isinstance(mage_tab, dict) else None
        modeled = render_model(model) if isinstance(model, dict) else None
        roundtrip = mage_tab.get("roundtrip") if isinstance(mage_tab, dict) else None
        core_changed = (
            isinstance(roundtrip, dict)
            and roundtrip.get("semantic_sha256") != semantic_sha256(data)
        )
        if modeled is not None and not core_changed:
            return modeled
        protocol_registry = ProtocolRegistry(series_accession=self._series_accession(data=data))
        technology_type = self._detect_ae_technology(data=data)
        sdrf = self.sdrf_constructor._miniml2sdrf(
            data=data,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        )
        idf = self.idf_constructor.miniml2idf(
            data=data,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        )
        sdrf_index = self._sdrf_row_index(rows=idf)
        if sdrf_index is None:
            raise ValueError("IDF does not contain an SDRF File row.")
        idf[sdrf_index] = ["SDRF File", sdrf, *idf[sdrf_index][2:]]
        if modeled is not None:
            return overlay_core(modeled, idf)
        return restore_extensions(data, idf)

    def _detect_ae_technology(self, data: dict) -> str:
        handler = JSONHandler()
        platform_tech = " ".join(str(x).lower() for x in handler._from_path(data, "platform.*.technology") if x)
        platform_text = " ".join(str(x).lower() for x in handler._from_path(data, "platform.*.title") if x)
        library_source = " ".join(str(x).lower() for x in handler._from_path(data, "sample.*.library_source") if x)
        library_strategy = " ".join(str(x).lower() for x in handler._from_path(data, "sample.*.library_strategy") if x)
        library_selection = " ".join(str(x).lower() for x in handler._from_path(data, "sample.*.library_selection") if x)
        sample_type = " ".join(str(x).lower() for x in handler._from_path(data, "sample.*.type") if x)
        text_paths = [
            "series.title",
            "series.summary",
            "series.overall_design",
            "series.type.*",
            "sample.*.description",
            "sample.*.data_processing",
            "sample.*.channel.*.extract_protocol",
            "sample.*.channel.*.growth_protocol",
            "sample.*.channel.*.treatment_protocol",
            "sample.*.channel.*.molecule",
            "sample.*.channel.*.characteristics.*.tag",
            "sample.*.channel.*.characteristics.*.value",
            "sample.*.supplementary_data.*.value",
            "sample.*.raw_data.*.value",
            "series.supplementary_data.*.value",
        ]
        text = " ".join(
            str(x).lower()
            for path in text_paths
            for x in handler._from_path(data, path)
            if x
        )
        relations = [x for x in handler._from_path(data, "sample.*.relation.*") if isinstance(x, dict)]
        has_sra = any((relation.get("type") or "").lower() == "sra" for relation in relations)
        has_array = "array" in platform_tech or self._has_array_files(data=data)

        if "high-throughput sequencing" in platform_tech or has_sra or library_strategy or sample_type == "sra":
            if "single cell" in library_source or "single-cell" in text or "single cell" in text or "10x" in text:
                if "visium" in text or "spatial" in text:
                    return "spatial_sequencing"
                if "10x" not in text and "droplet" not in text and "chromium" not in text:
                    return "plate_single_cell_sequencing"
                if self._has_tenx_version(text=text, version="3"):
                    return "tenx_v3_droplet_single_cell_sequencing"
                if self._has_tenx_version(text=text, version="2"):
                    return "tenx_v2_droplet_single_cell_sequencing"
                return "droplet_single_cell_sequencing"
            return "bulk_sequencing"

        if has_array:
            return "array"

        return "generic"

    def _has_tenx_version(self, text: str, version: str) -> bool:
        if "10x" not in text and "chromium" not in text:
            return False
        return re.search(rf"(?<![a-z0-9])v{version}(?![a-z0-9])", text) is not None

    def _has_array_files(self, data: dict) -> bool:
        handler = JSONHandler()
        values = []
        for path in (
            "platform.*.supplementary_data.*.value",
            "sample.*.supplementary_data.*.value",
            "sample.*.raw_data.*.value",
            "series.supplementary_data.*.value",
        ):
            values.extend(x for x in handler._from_path(data, path) if x)

        array_extensions = (".cel", ".gpr", ".idat", ".chp", ".txt", ".tif", ".tiff", ".exp", ".rpt", ".cab")
        return any(
            normalized_extension(value) in array_extensions
            for value in values
        )

    def _series_accession(self, data: dict):
        series = data.get("series") if isinstance(data, dict) else None
        series_values = series if isinstance(series, list) else [series]
        for series_item in series_values:
            if not isinstance(series_item, dict):
                continue
            accessions = series_item.get("accession")
            accession_values = accessions if isinstance(accessions, list) else [accessions]
            for accession in accession_values:
                if isinstance(accession, dict) and accession.get("value"):
                    return accession.get("value")
        return "GEO"
    
    def magetab2file(self, magetab:list, out:str = None) -> str:
        '''
        Write magetab to idf and sdrf
        '''
        out = out or "."
        os.makedirs(out, exist_ok=True)

        rows = self._normalize_magetab_rows(magetab=copy.deepcopy(magetab))
        sdrf_index = self._sdrf_row_index(rows=rows)
        if sdrf_index is None:
            raise ValueError("MAGETAB does not contain an SDRF File row.")

        sdrf_row = rows[sdrf_index]
        if len(sdrf_row) < 2:
            raise ValueError("SDRF File row does not contain an SDRF payload.")

        sdrf = sdrf_row[1]
        if not self._is_table(sdrf):
            raise ValueError("SDRF payload must be a non-empty list of row lists.")

        ae_accession = self._magetab_accession(rows=rows)
        idf_filename = f"{ae_accession}.idf.txt"
        sdrf_filename = f"{ae_accession}.sdrf.txt"
        idf_path = os.path.join(out, idf_filename)
        sdrf_path = os.path.join(out, sdrf_filename)

        rows[sdrf_index] = ["SDRF File", sdrf_filename, *sdrf_row[2:]]
        idf_rows = rows
        if not idf_rows:
            raise ValueError("MAGETAB does not contain usable IDF rows.")

        self._write_tsv(path=idf_path, rows=idf_rows)
        self._write_tsv(path=sdrf_path, rows=sdrf)

        return idf_path

    def _normalize_magetab_rows(self, magetab: list) -> list:
        rows = []
        index = 0
        while index < len(magetab):
            item = magetab[index]

            if (
                isinstance(item, str)
                and item.strip().lower() == "sdrf file"
                and index + 1 < len(magetab)
            ):
                rows.append(["SDRF File", magetab[index + 1]])
                index += 2
                continue

            if isinstance(item, str):
                if "," in item:
                    key, value = item.split(",", 1)
                    rows.append([key.strip(), value.strip()])
                elif item.strip():
                    rows.append([item.strip()])
                index += 1
                continue

            if isinstance(item, (list, tuple)):
                rows.append(list(item))
                index += 1
                continue

            rows.append([item])
            index += 1

        return rows

    def _strip_quotes_from_table(self, rows: list) -> list:
        return [self._strip_quotes(row) for row in rows]

    def _strip_quotes(self, value):
        if isinstance(value, list):
            return [self._strip_quotes(item) for item in value]
        if isinstance(value, tuple):
            return tuple(self._strip_quotes(item) for item in value)
        if isinstance(value, str):
            return value.replace('"', "").replace("'", "")
        return value

    def _sdrf_row_index(self, rows: list):
        for index, row in enumerate(rows):
            if row and str(row[0]).strip().lower() == "sdrf file":
                return index
        return None

    def _magetab_accession(self, rows: list) -> str:
        labels = [
            "comment[arrayexpressaccession]",
            "investigation accession",
            "comment[secondaryaccession]",
        ]
        rows_by_label = {
            str(row[0]).strip().lower(): row
            for row in rows
            if row
        }

        for label in labels:
            row = rows_by_label.get(label)
            if not row:
                continue
            for value in row[1:]:
                if value is not None and str(value).strip():
                    return self._safe_filename_token(value=str(value).strip())

        return "AE"

    def _safe_filename_token(self, value: str) -> str:
        token = value.replace(os.sep, "_")
        if os.altsep:
            token = token.replace(os.altsep, "_")
        return token or "AE"

    def _is_table(self, value) -> bool:
        return (
            isinstance(value, list)
            and bool(value)
            and all(isinstance(row, (list, tuple)) for row in value)
        )

    def _write_tsv(self, path: str, rows: list) -> None:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for row in rows:
                writer.writerow(["" if value is None else value for value in row])
