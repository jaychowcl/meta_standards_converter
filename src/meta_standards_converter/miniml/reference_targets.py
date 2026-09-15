"""Field-contract parsing shared by reference import and residual pruning."""
import re

_LIST_FIELDS = {'ENA-STUDY','ENA-SAMPLE','ENA-EXPERIMENT','ENA-RUN','ENA-ANALYSIS','ENA-SUBMISSION'}


def parse_reference_targets(database, literal, verified_ranges=()):
    if re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:', literal):
        return [literal]
    for record in verified_ranges:
        if record.get('database') == database and record.get('literal') == literal and record.get('accessions'):
            return list(record['accessions'])
    parts = literal.split(',')
    if database in _LIST_FIELDS and all(re.fullmatch(r'[SED]R[PSXRZA]\d+', p.strip()) for p in parts):
        return [p.strip() for p in parts]
    return [literal]
