# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Native import outcomes and independent atomic file publication."""
from meta_standards_converter.conversion_logging import measured, count, progress
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


def _load_native(source, parser, seed):
    """Retain a strictly valid core before optional provider discovery can stop."""
    from inspect import signature
    from ..metadata.retrieval_budget import _active, _phase, _budget_stage, _BudgetExhausted, _coverage
    from ..miniml import MINiMLCodec
    from ..miniml.codec import MINiMLCompatibilityError
    budget, candidate, captured = _active.get(), None, None
    def core(records):
        nonlocal candidate, captured
        captured = records
        if not records.xml and not any(records.indexed.values()): return
        try: candidate = MINiMLCodec().decode(parser.parse(records), strict=True).package
        except MINiMLCompatibilityError: pass
        _phase.set('native_supplementary')  # Restored by the enclosing stage scope.
    try:
        with _budget_stage('native_core'):
            options = {'_on_core': core} if budget and '_on_core' in signature(source.fetch).parameters else {}
            records = source.fetch(seed, **options)
        if not records.xml and not any(records.indexed.values()):
            raise ValueError('No native metadata records retrieved')
        return parser.parse(records), records.issues
    except _BudgetExhausted as exc:
        if candidate is None: raise
        try:
            candidate = MINiMLCodec().decode(parser.parse(captured), strict=True).package
        except Exception:
            pass  # Supplementary records cannot invalidate the accepted core.
        return _coverage(candidate, budget), [*captured.issues, str(exc) + '; valid native candidate retained']


@measured("publication")
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
        count("artifacts_committed")
        progress(1, 1, "artifacts")
    finally:
        temporary.unlink(missing_ok=True)


def output_name(seed, seeds):
    value = seed.primary
    if sum(other.primary == seed.primary for other in seeds) > 1:
        value += '__' + seed.study
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
        raise ValueError('Unsafe study filename')
    return value + '.json'
