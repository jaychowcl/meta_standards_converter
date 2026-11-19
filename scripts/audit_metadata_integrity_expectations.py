# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Generate separate review candidates; never replace historical golden files."""
import argparse
from datetime import date
import json
import logging
from pathlib import Path
import subprocess
import tempfile

from pytest import MonkeyPatch
from scripts.audit_authorship_candidate import sha256
from tests.support.contracts import GEO, artifacts, install_replay
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy

ROOT = Path(__file__).resolve().parents[1]


def generate(destination):
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError('Use a fresh reviewed-expectation directory')
    destination.mkdir(parents=True)
    records = []
    def save(path, value, original, reason):
        target = destination / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n' if isinstance(value, (dict, list)) else value)
        records.append({'original': str(original.relative_to(ROOT)), 'original_sha256': sha256(original),
                        'derivative': str(path), 'sha256': sha256(target), 'reason': reason})
    with tempfile.TemporaryDirectory(prefix='msc-expectation-') as temporary, MonkeyPatch.context() as patches:
        workspace = Path(temporary)
        patches.chdir(workspace)
        install_replay(patches)
        for command, source, target in [('geo2json', 'GSE328265', 'json'), ('geo2ae', 'GSE328265', 'magetab'),
                                        ('json2ae', GEO / 'inputs/miniml.json', 'magetab')]:
            policy = PreparationPolicy(enrichment='off' if command == 'json2ae' else 'standard',
                                       expand_studies=False, linked_metadata=False, peer_provider=False)
            options = {'platform_handler': 'single_cell_sequencing'} if target == 'magetab' else {}
            converted(Converter(), source, out_type=target, preparation=policy,
                      output_options=options, outdir=workspace / command)
            for name, value in artifacts(workspace / command, workspace).items():
                save(Path('GSE328265') / command / name, value, GEO / 'expected' / command / name,
                     'Canonical graph and source-file provenance; unknown file links retained; explicit read roles only.')
        # This routing fixture contains growth and treatment. Change the ordered
        # application sequence while retaining every definition and other cell.
        original = ROOT / 'tests/fixtures/chemistry/routing/GSM9254695.sdrf.txt'
        text = original.read_text()
        assert text.count('P-GSE900012-1\tP-GSE900012-2') == 1
        save(Path('GSM9254695.sdrf.txt'), text.replace('P-GSE900012-1\tP-GSE900012-2', 'P-GSE900012-2\tP-GSE900012-1'),
             original, 'GEO-generated growth precedes treatment; explicit definitions and all other fields unchanged.')
        # Protocol identifiers are opaque graph-local keys. This independent
        # expected transformation renames the sole generated blank definition
        # and its two applications together; no scientific value changes.
        original = ROOT / 'tests/fixtures/sra_publication_concurrency.json'
        value = json.loads(original.read_text())
        assert value['magetab'][30][1] == 'P-SRP250911-2'
        value['magetab'][30][1] = 'P-SRP250911-1'
        for row in value['magetab'][40][1][1:]:
            assert row[14] == 'P-SRP250911-2'
            row[14] = 'P-SRP250911-1'
        save(Path(original.name), value, original, 'One blank protocol reference and two applications consistently renamed after redundant graph construction was removed.')
    manifest = {'production_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'historical_fixtures_unchanged': True, 'files': records,
                'review_requirement': 'Generated GSE328265 candidates require source/field review and independent regressions before acceptance.'}
    (destination / 'review.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    generate(args.destination)
