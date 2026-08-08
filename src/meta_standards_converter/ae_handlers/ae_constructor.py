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
from meta_standards_converter.ae_handlers.ae_common import (
    ProtocolRegistry,
    _has_tenx_version,
    detect_ae_technology,
    has_array_files,
    normalized_extension,
)

import copy
import csv
import os


PLATFORM_HANDLER_KEYS = (
    "plate_single_cell_sequencing",
    "droplet_single_cell_sequencing",
    "tenx_v2_droplet_single_cell_sequencing",
    "tenx_v3_droplet_single_cell_sequencing",
    "single_cell_sequencing",
    "spatial_sequencing",
    "bulk_sequencing",
    "sequencing",
    "array",
    "generic",
)


def validate_platform_handler(value: str) -> str:
    if value not in PLATFORM_HANDLER_KEYS:
        raise ValueError(
            f"Unsupported platform handler: {value}. "
            f"Choose one of: {', '.join(PLATFORM_HANDLER_KEYS)}"
        )
    return value


from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage


class AEConstructor:
    def __init__(self, idf_constructor=None, sdrf_constructor=None):
        self.idf_constructor = idf_constructor or IDFConstructor()
        self.sdrf_constructor = sdrf_constructor or SDRFConstructor()

    def miniml2magetab(self, data: MINiMLPackage, platform_handler: str | None = None) -> list:
        """
        converts miniml json to magetab idf. Walks through sections of idf to extract from miniml
        """
        data = MINiMLCodec().encode(MINiMLCodec().decode(data).package)
        forced = platform_handler is not None
        if forced:
            technology_type = validate_platform_handler(platform_handler)
            modeled = None
        else:
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
            technology_type = self._detect_ae_technology(data=data)
        protocol_registry = ProtocolRegistry(series_accession=self._series_accession(data=data))
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
        return detect_ae_technology(data)

    def _has_tenx_version(self, text: str, version: str) -> bool:
        return _has_tenx_version(text, version)

    def _has_array_files(self, data: dict) -> bool:
        return has_array_files(data)

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
