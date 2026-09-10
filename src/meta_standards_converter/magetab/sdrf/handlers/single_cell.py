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
from meta_standards_converter.magetab.sdrf.handlers.sequencing import _SequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.sequencing import _BulkSequencingSDRFHandler

class _SingleCellSequencingSDRFHandler(_SequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = []
        values = {
            "Comment[library construction]": self.library_construction(sample=sample),
        }
        read_lengths = run.get("read_lengths") if run else None
        if read_lengths:
            values["Comment[cdna read size]"] = read_lengths[-1]
        for label, value in values.items():
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def extra_assay_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        biosample = self.biosample_accession(sample=sample)
        return [SDRFAttr(label="Comment[technical replicate group]", value=biosample)] if biosample else []

    def library_construction(self, sample: dict):
        return None

    def study_text(self, sample: dict) -> str:
        values = [
            sample.get("title"),
            sample.get("description"),
            sample.get("data_processing"),
        ]
        for channel in self._as_list(sample.get("channel")):
            if isinstance(channel, dict):
                values.extend([
                    channel.get("extract_protocol"),
                    channel.get("growth_protocol"),
                    channel.get("treatment_protocol"),
                ])
        for series in self._as_list(self.data.get("series")):
            if not isinstance(series, dict):
                continue
            values.extend([
                series.get("title"),
                series.get("summary"),
                series.get("overall_design"),
            ])
        return " ".join(str(value).lower() for value in values if value)


class _DropletSingleCellSequencingSDRFHandler(_SingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = super().extra_library_attrs(sample=sample, channel=channel, run=run)
        text = self.study_text(sample=sample)
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode size]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode size]": "12",
            "Comment[sample barcode read]": "index1",
            "Comment[single cell isolation]": self.single_cell_isolation(sample=sample),
        }
        for label, value in values.items():
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def library_construction(self, sample: dict):
        text = self.study_text(sample=sample)
        if "v3" in text or "3' v3" in text or "3 v3" in text:
            return "10xV3"
        if "10x" in text or "chromium" in text or "droplet" in text:
            return "10x technology"
        return None

    def single_cell_isolation(self, sample: dict):
        text = self.study_text(sample=sample)
        if "10x" in text or "chromium" in text:
            return "10x technology"
        if "droplet" in text:
            return "droplet"
        return None


class _TenXDropletSingleCellSequencingSDRFHandler(_DropletSingleCellSequencingSDRFHandler):
    def tenx_library_attrs(
        self,
        library_construction: str,
        cdna_read_size: str,
        umi_barcode_size: str,
    ) -> list[SDRFAttr]:
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cdna read offset]": "0",
            "Comment[cdna read size]": cdna_read_size,
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode size]": "16",
            "Comment[end bias]": "3 prime tag",
            "Comment[input molecule]": "polyA RNA",
            "Comment[library construction]": library_construction,
            "Comment[primer]": "oligo-dT",
            "Comment[LIBRARY_STRAND]": "not applicable",
            "Comment[sample barcode offset]": "0",
            "Comment[sample barcode read]": "index1",
            "Comment[sample barcode size]": "8",
            "Comment[single cell isolation]": "10x technology",
            "Comment[spike in]": "",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode size]": umi_barcode_size,
        }
        return [
            SDRFAttr(label=label, value=self.clean(value), required=value == "")
            for label, value in values.items()
        ]


class _TenXV2DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        return self.tenx_library_attrs(
            library_construction="10xV2",
            cdna_read_size="98",
            umi_barcode_size="10",
        )


class _TenXV3DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        return self.tenx_library_attrs(
            library_construction="10xV3",
            cdna_read_size="91",
            umi_barcode_size="12",
        )


class _PlateSingleCellSequencingSDRFHandler(_BulkSequencingSDRFHandler):
    def extra_source_attrs(self, sample: dict, channel: dict) -> list[SDRFAttr]:
        barcode = self.clean(sample.get("barcode"))
        description = self.clean(sample.get("description"))
        attrs = []
        if barcode:
            attrs.append(SDRFAttr(label="Comment[index]", value=barcode))
        if description:
            attrs.append(SDRFAttr(label="Description", value=description))
        return attrs
