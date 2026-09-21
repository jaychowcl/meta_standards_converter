# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Write constructed MAGE-TAB tables to IDF/SDRF files."""
import copy
import csv
import os

class MAGETabWriter:
    def write(self, magetab:list, out:str = None) -> str:
        '''
        Write magetab to idf and sdrf
        '''
        out = out or "."
        rows = self._normalize_magetab_rows(magetab=copy.deepcopy(magetab))
        from .validation import validate_magetab
        validate_magetab(rows)
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

        os.makedirs(out, exist_ok=True)
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
