# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""One isolated, warmed, network-free baseline or candidate measurement."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import time


def measure(source, output, manifest):
    from tests.support.recorded_replay import ReplayTransport
    expected_origin = Path(importlib.util.find_spec('meta_standards_converter').origin).resolve()
    root = expected_origin.parents[2]
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True)
    replay = ReplayTransport(manifest)
    with replay.offline():
        import meta_standards_converter
        from meta_standards_converter.converters import Converter, PreparationPolicy
        from meta_standards_converter.conversion_logging import operation
        origin = Path(meta_standards_converter.__file__).resolve()
        assert origin == expected_origin
        policy = PreparationPolicy(enrichment='off', expand_studies=False)
        before = hashlib.sha256(source.read_bytes()).hexdigest()
        class Emitter:
            def __init__(self): self.events = []
            def emit(self, event_type, **values): self.events.append({'event_type': event_type, **values})
        warm = Emitter()
        with operation('audit_timing_warmup', heartbeat=False, emitter=warm):
            result = Converter().convert(source, out_type='magetab', preparation=policy)
        assert result.items[0].execution == 'succeeded'
        events = Emitter()
        wall, cpu = time.perf_counter(), time.process_time()
        with operation('audit_timing', heartbeat=False, emitter=events):
            result = Converter().convert(source, out_type='magetab', preparation=policy)
        wall, cpu = time.perf_counter() - wall, time.process_time() - cpu
        assert result.items[0].execution == 'succeeded'
        assert hashlib.sha256(source.read_bytes()).hexdigest() == before
        tables = result.items[0].payload
        def primitive(value):
            if hasattr(value, 'to_mapping'): return value.to_mapping()
            if isinstance(value, (tuple, list)): return [primitive(v) for v in value]
            if isinstance(value, dict): return {k: primitive(v) for k, v in value.items()}
            return value
        payload = primitive(tables)
        counts = [{'idf_rows': len(table), 'sdrf_rows': len(next(r[1] for r in table if r[0] == 'SDRF File')) - 1}
                  for table in payload]
        record = {'source': str(source), 'source_sha256': before, 'revision': revision,
                  'msc_origin': str(origin), 'wall_seconds': wall, 'cpu_seconds': cpu,
                  'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                  'peak_rss_scope': 'process lifetime including one warm-up',
                  'output_sha256': hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(),
                  'counts': counts, 'status': result.status,
                  'items': [{'execution': i.execution, 'completeness': i.completeness, 'validation': i.validation} for i in result.items],
                  'stage_events': [e for e in events.events if e['event_type'] == 'stage_finished'],
                  'live_requests': 0}
    replay.assert_complete()
    with output.open('x') as stream: json.dump(record, stream, indent=2)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'output', 'manifest'): parser.add_argument('--' + name, type=Path, required=True)
    measure(**vars(parser.parse_args()))
