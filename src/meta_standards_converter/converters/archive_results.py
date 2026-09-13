"""Native import outcomes and independent atomic file publication."""
from dataclasses import dataclass, field
from pathlib import Path
import json
import os
import re
from uuid import uuid4
from meta_standards_converter.miniml import MINiMLPackage


@dataclass
class StudyImportOutcome:
    study: str
    primary: str
    status: str
    package: MINiMLPackage | None = None
    output: str | None = None
    issues: list[str] = field(default_factory=list)

    def to_mapping(self):
        return {'study': self.study, 'primary': self.primary, 'status': self.status,
                'output': self.output, 'issues': list(self.issues)}


@dataclass
class ArchiveImportResult:
    accession: str
    studies: list[StudyImportOutcome] = field(default_factory=list)

    @property
    def packages(self):
        return [s.package for s in self.studies if s.package is not None]

    @property
    def ok(self):
        return all(s.status == 'complete' for s in self.studies)

    def to_mapping(self):
        return {'accession': self.accession, 'studies': [s.to_mapping() for s in self.studies]}


def publish_json(path, value, overwrite=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid4().hex + '.tmp')
    try:
        with temporary.open('x', encoding='utf-8') as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)  # atomic no-replace, including competing imports
    finally:
        temporary.unlink(missing_ok=True)


def output_name(seed, seeds):
    value = seed.primary
    if sum(other.primary == seed.primary for other in seeds) > 1:
        value += '__' + seed.study
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
        raise ValueError('Unsafe study filename')
    return value + '.json'
