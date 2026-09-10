# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from meta_standards_converter.magetab.sdrf.model import SDRFAttr
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _SingleCellSequencingSDRFHandler

class _SpatialSequencingSDRFHandler(_SingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        # Spatial presets are outside the scoped single-cell chemistry change.
        attrs = []
        construction = self.library_construction(sample)
        if construction:
            attrs.append(SDRFAttr(label="Comment[library construction]", value=construction))
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode size]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode size]": "12",
            "Comment[sample barcode read]": "index1",
        }
        for label, value in values.items():
            attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def library_construction(self, sample: dict):
        text = self.study_text(sample=sample)
        if "visium" in text:
            return "10x Visium"
        return None

    def extra_assay_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = super().extra_assay_attrs(sample=sample, run=run)
        filename = self.clean(run.get("submitted_file_name") if run else None) or ""
        lower_filename = filename.lower()
        if "_i1_" in lower_filename or "index" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="sample_barcode"),
                SDRFAttr(label="Comment[read_index]", value="index1"),
            ])
        elif "_r1_" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="cell_barcode"),
                SDRFAttr(label="Comment[read_index]", value="read1"),
            ])
        elif "_r2_" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="cdna"),
                SDRFAttr(label="Comment[read_index]", value="read2"),
            ])
        return attrs
