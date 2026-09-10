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
        from meta_standards_converter.magetab.chemistry import resolve_chemistry
        result = resolve_chemistry(sample, channel=channel, run=run, series=self.data.get("series"))
        for diagnostic in result.diagnostics:
            message = f"Sample {self.sample_accession(sample=sample)} chemistry {diagnostic.code}: {', '.join(diagnostic.paths)}"
            if message not in self.audit.warnings:
                self.audit.warnings.append(message)
        return [SDRFAttr(label=f"Comment[{name}]", value=value) for name, value in result.attributes]

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
    """Droplet rendering; chemistry attributes require source evidence."""


class _TenXDropletSingleCellSequencingSDRFHandler(_DropletSingleCellSequencingSDRFHandler):
    """Retained rendering selection without a scientific preset."""


class _TenXV2DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    """Legacy selection name; a version alone does not identify chemistry."""


class _TenXV3DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    """Legacy selection name; a version alone does not identify chemistry."""


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
