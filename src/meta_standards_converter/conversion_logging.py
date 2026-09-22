# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operation-local telemetry shared by unified and standalone conversions.

This module never installs handlers or changes application logging levels.
"""
from collections import Counter
from contextlib import contextmanager
from contextvars import ContextVar, copy_context
from functools import wraps
import logging
import threading
import time
from uuid import uuid4

from .operational_events import OperationalEventEmitter, redact
from .runtime_contracts import SafeErrorEnvelope

_run = ContextVar("msc_log_run", default=None)
_labels = ContextVar("msc_log_labels", default={})
_span = ContextVar("msc_log_span", default=None)


def current():
    return _run.get()


def context_callable(fn):
    """Capture the submitting operation when dispatching a worker."""
    context = copy_context()
    return lambda *args, **kwargs: context.run(fn, *args, **kwargs)


@contextmanager
def labels(**values):
    token = _labels.set(_labels.get() | {k: v for k, v in values.items() if v is not None})
    try:
        yield
    finally:
        _labels.reset(token)


def bind(**attributes):
    """Bind additional identity inside the enclosing operation/item scope."""
    _labels.set(_labels.get() | {k: v for k, v in attributes.items() if v is not None})


def count(name, amount=1):
    run = current()
    if run is not None:
        with run.lock:
            run.counters[name] += amount


def emit(event_type, *, level="INFO", status=None, duration_seconds=None, **attributes):
    run = current()
    if run is None:
        return None
    span = _span.get()
    return run.emitter.emit(
        event_type, component=run.component, level=level, run_id=run.id,
        stage=attributes.pop("stage", None) or (span.name if span else None),
        status=status, duration_seconds=duration_seconds,
        attributes=_labels.get() | attributes,
    )


class Operation:
    def __init__(self, component, clock, emitter=None):
        self.component = component
        self.clock = clock
        self.id = uuid4().hex
        self.emitter = emitter or OperationalEventEmitter(repository="meta_standards_converter")
        self.counters = Counter()
        self.lock = threading.RLock()
        self.spans = []
        self.started = clock()
        self.last_tick = self.started
        self.failure_stage = None
        self.failure_error = None
        self.status = "complete"
        self.statistics = {}
        self.results_seen = set()

    def tick(self):
        with self.lock:
            now = self.clock()
            if now - self.last_tick < 30 or not self.spans:
                return
            self.last_tick = now
            span = self.spans[-1]
            self.emitter.emit("progress", component=self.component, run_id=self.id,
                stage=span.name, status="running", duration_seconds=now - span.started,
                attributes=span.labels | span.statistics | dict(self.counters))


@contextmanager
def operation(component, *, clock=time.monotonic, heartbeat=True, emitter=None):
    existing = current()
    if existing is not None:
        yield existing
        return
    run = Operation(component, clock, emitter)
    token = _run.set(run)
    label_token = _labels.set(dict(_labels.get()))
    stop = threading.Event()
    worker = None
    def pulse():
        while not stop.wait(1):
            run.tick()
    try:
        emit("conversion_started", status="running")
        if heartbeat:
            worker = threading.Thread(target=context_callable(pulse), daemon=True, name="msc-progress")
            worker.start()
        yield run
    except BaseException as error:
        run.status = "interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed"
        safe = SafeErrorEnvelope.from_exception(error, stage=run.failure_stage,
                                                correlation_id=run.id)
        run.statistics["error"] = safe.to_dict()
        try:
            error._msc_logged_envelope = safe
        except (AttributeError, TypeError):
            pass
        raise
    finally:
        stop.set()
        if worker is not None:
            worker.join(timeout=2)
        try:
            level = "ERROR" if run.status in {"failed", "interrupted"} else "WARNING" if run.status == "partial" else "INFO"
            emit("conversion_finished", level=level, status=run.status,
                 stage=run.failure_stage if run.status in {"failed", "interrupted"} else None, duration_seconds=clock() - run.started,
                 **(dict(run.counters) | run.statistics))
        finally:
            _labels.reset(label_token)
            _run.reset(token)


class Stage:
    def __init__(self, name, statistics):
        self.name = name
        self.statistics = statistics
        self.labels = dict(_labels.get())
        self.status = "completed"
        self.started = 0

    def update(self, **statistics):
        run = current()
        if run is not None:
            with run.lock:
                self.statistics.update(statistics)


@contextmanager
def stage(name, **statistics):
    span = Stage(name, statistics)
    run = current()
    if run is None:
        yield span
        return
    parent = _span.get()
    token = _span.set(span)
    phases = {"discovery", "input_validation", "load", "validation", "expansion",
              "preparation", "conversion", "projection", "expression", "publication"}
    level = "INFO" if parent is None or (name in phases and name != parent.name) else "DEBUG"
    span.started = run.clock()
    with run.lock:
        run.spans.append(span)
    emit("stage_started", level=level, status="running", **statistics)
    try:
        yield span
    except BaseException as error:
        span.status = "failed"
        if run.failure_error is not error:
            run.failure_stage = name
            run.failure_error = error
        raise
    finally:
        try:
            emit("stage_finished", level="DEBUG" if span.status == "failed" else level,
                 status=span.status, duration_seconds=run.clock() - span.started,
                 **span.statistics)
        finally:
            with run.lock:
                run.spans.remove(span)
            _span.reset(token)


def measured(name):
    """Measure a stage without changing its signature or return contract."""
    def decorate(fn):
        @wraps(fn)
        def wrapped(*args, **kwargs):
            with stage(name) as span:
                result = fn(*args, **kwargs)
                span.update(**statistics(result))
                return result
        return wrapped
    return decorate


def statistics(value):
    """Read bounded public result metadata; never inspect scientific matrices."""
    stats = {}
    metadata = getattr(value, "metadata", None)
    if metadata is not None and hasattr(metadata, "groups"):
        value = metadata
    if hasattr(value, "samples") and hasattr(value, "to_mapping"):
        value = [value]
    groups = getattr(value, "groups", None)
    packages = getattr(value, "packages", None)
    if isinstance(groups, (tuple, list)):
        stats["datasets"] = len(groups)
        packages = [p for g in groups for p in g.packages]
    if isinstance(value, (tuple, list)) and value and hasattr(value[0], "samples"):
        packages = value
    if isinstance(packages, (tuple, list)):
        stats["packages"] = len(packages)
        stats["samples"] = sum(len(getattr(p, "samples", ())) for p in packages)
        stats["runs"] = sum(len(sample.to_mapping().get("sra_run", ()))
                            for p in packages for sample in getattr(p, "samples", ()))
    for attr, key in (("row_count", "rows"), ("n_obs", "observations"), ("n_vars", "features")):
        number = getattr(value, attr, None)
        if isinstance(number, int):
            stats[key] = number
    artifacts = getattr(value, "artifacts", None)
    if isinstance(artifacts, dict):
        stats["artifacts"] = len(artifacts)
    items = getattr(value, "items", None)
    if isinstance(items, list):
        stats["inputs"] = len(items)
        stats.update({state: sum(i.status == state for i in items)
                      for state in ("complete", "partial", "failed", "skipped")})
        stats["artifacts"] = sum(len(i.artifacts) for i in items)
    studies = getattr(value, "studies", None)
    if isinstance(studies, list):
        stats["studies"] = len(studies)
        stats.update({state: sum(i.status == state for i in studies)
                      for state in ("complete", "partial", "failed")})
    if isinstance(value, (list, tuple)):
        tables = value if value and isinstance(value[0], (list, tuple)) and value[0] and isinstance(value[0][0], (list, tuple)) else [value]
        row_counts = [len(row[1]) - 1 for table in tables for row in table
                      if isinstance(row, (list, tuple)) and len(row) > 1
                      and row[0] == "SDRF File" and isinstance(row[1], list)]
        if row_counts:
            stats.update(magetab_packages=len(row_counts), sdrf_rows=sum(row_counts))
    sample_files = getattr(value, "sample_h5ads", None)
    if isinstance(sample_files, dict):
        stats["sample_outputs"] = len(sample_files)
    obs = getattr(value, "obs", None)
    if hasattr(obs, "shape") and isinstance(obs.shape, tuple):
        stats["observations"] = obs.shape[0]
    columns = getattr(value, "columns", None)
    if isinstance(columns, (tuple, list)):
        stats["columns"] = len(columns)
    return stats


def progress(completed, total=None, unit="units", **attributes):
    span = _span.get()
    if span is not None:
        span.update(completed=completed, total=total, unit=unit, **attributes)
    emit("progress", level="DEBUG", completed=completed, total=total, unit=unit, **attributes)


def diagnostic(code, message, *, severity="warning", stage_name=None, **attributes):
    """Retain each diagnostic occurrence; sinks control presentation."""
    run = current()
    if run is None:
        return
    count(severity + "s")
    emit("diagnostic", level=severity.upper(), stage=stage_name, code=code,
         message=redact(message), **attributes)


def result_diagnostics(result):
    run = current()
    if run is not None:
        with run.lock:
            if id(result) in run.results_seen:
                return
            run.results_seen.add(id(result))
    items = getattr(result, "items", None)
    if isinstance(items, list):
        for item in items:
            with labels(item_id=item.id):
                result_diagnostics(item)
    diagnostics = getattr(result, "diagnostics", None)
    if isinstance(diagnostics, (tuple, list)):
        for d in diagnostics:
            severity = "debug" if d.code in {"duplicate_input", "fail_fast", "enrichment_skipped"} else "error" if d.code in {"sample_skipped", "source_failed"} else d.severity
            diagnostic(d.code, d.message, severity=severity, stage_name=d.stage)
    for attr, severity in (("warnings", "warning"), ("errors", "error"), ("failures", "error"), ("issues", "warning")):
        values = getattr(result, attr, None)
        if isinstance(values, (tuple, list)):
            for value in values:
                diagnostic(attr, str(value), severity="error" if attr == "issues" and getattr(result, "status", None) == "failed" else severity)
    for attr in ("studies", "conversions"):
        values = getattr(result, attr, None)
        if isinstance(values, dict):
            values = list(values.values())
        if isinstance(values, (tuple, list)):
            for value in values:
                with labels(dataset_id=getattr(value, "study", None) or getattr(value, "study_accession", None)):
                    result_diagnostics(value)


def outcome(result):
    status = getattr(result, "status", None)
    if isinstance(status, str) and status in {"complete", "partial", "failed", "skipped"}:
        return status
    studies = getattr(result, "studies", None)
    if isinstance(studies, list):
        if studies and all(s.status == "failed" for s in studies):
            return "failed"
        return "complete" if all(s.status == "complete" for s in studies) else "partial"
    if getattr(result, "partial", False) is True:
        return "partial"
    return "complete"


def logged_conversion(fn):
    @wraps(fn)
    def wrapped(*args, **kwargs):
        if current() is not None:
            emit("converter_delegated", level="DEBUG", converter=fn.__qualname__)
            with labels():
                return fn(*args, **kwargs)
        with operation(fn.__qualname__) as run:
            source = next((kwargs[k] for k in ("input", "source", "json_path", "gse", "accession") if k in kwargs), args[1] if len(args)>1 else None)
            if source is not None:
                label = source if isinstance(source, str) and not source.lstrip().startswith(("{", "[", "<")) else type(source).__name__
                bind(source=label)
            result = fn(*args, **kwargs)
            result_diagnostics(result)
            run.status = outcome(result)
            run.statistics.update(statistics(result))
            return result
    return wrapped


@contextmanager
def item_scope(item, *, total, index):
    run = current()
    if run is None:
        yield
        return
    run.failure_stage = None
    with labels(item_id=item.id):
        emit("item_started", status="running", index=index, inputs=total, source=item.source)
        try:
            yield
        finally:
            result_diagnostics(item)
            count("inputs_finished")
            emit("item_finished", status=item.status, index=index, inputs=total,
                 finished=index, execution=item.execution, completeness=item.completeness,
                 validation=item.validation, datasets=len(item.dataset_ids), artifacts=len(item.artifacts))
