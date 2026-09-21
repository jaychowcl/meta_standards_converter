"""One strict MINiML policy at every metadata acceptance boundary."""
from meta_standards_converter.miniml import MINiMLCodec


def validate_metadata(metadata):
    if metadata is not None:
        codec = MINiMLCodec()
        for group in metadata.groups:
            for package in group.packages:
                codec.decode(package, strict=True)
