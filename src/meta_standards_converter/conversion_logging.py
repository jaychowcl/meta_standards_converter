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
        self.status = "complete"
        self.statistics = {}
        self.diagnostics_seen = set()

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
        raise
    finally:
        stop.set()
        if worker is not None:
            worker.join(timeout=2)
        try:
            level = "ERROR" if run.status in {"failed", "interrupted"} else "WARNING" if run.status == "partial" else "INFO"
            emit("conversion_finished", level=level, status=run.status,
                 stage=run.failure_stage, duration_seconds=clock() - run.started,
                 **(dict(run.counters) | run.statistics))
        finally:
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
    level = "DEBUG" if parent is not None else "INFO"
    span.started = run.clock()
    with run.lock:
        run.spans.append(span)
    emit("stage_started", level=level, status="running", **statistics)
    try:
        yield span
    except BaseException:
        span.status = "failed"
        run.failure_stage = name
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
    groups = getattr(value, "groups", None)
    packages = getattr(value, "packages", None)
    if isinstance(groups, (tuple, list)):
        stats["datasets"] = len(groups)
        packages = [p for g in groups for p in g.packages]
    if isinstance(value, (tuple, list)) and value and hasattr(value[0], "samples"):
        packages = value
    if isinstance(packages, (tuple, list)):
        stats["packages"] = len(packages)
        stats["samples"] = sum(len(p.samples) for p in packages)
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
    return stats


def diagnostic(code, message, *, severity="warning", stage_name=None, **attributes):
    """Record a diagnostic once per logical item; sinks control presentation."""
    run = current()
    if run is None:
        return
    key = (_labels.get().get("item_id"), stage_name, code, message)
    with run.lock:
        if key in run.diagnostics_seen:
            return
        run.diagnostics_seen.add(key)
    count(severity + "s")
    emit("diagnostic", level=severity.upper(), stage=stage_name, code=code,
         message=redact(message), **attributes)


def result_diagnostics(result):
    items = getattr(result, "items", None)
    if isinstance(items, list):
        for item in items:
            with labels(item_id=item.id):
                result_diagnostics(item)
    diagnostics = getattr(result, "diagnostics", None)
    if isinstance(diagnostics, (tuple, list)):
        for d in diagnostics:
            diagnostic(d.code, d.message, severity=d.severity, stage_name=d.stage)
    for attr, severity in (("warnings", "warning"), ("errors", "error"), ("failures", "error"), ("issues", "warning")):
        values = getattr(result, attr, None)
        if isinstance(values, (tuple, list)):
            for value in values:
                diagnostic(attr, str(value), severity=severity)
    for attr in ("studies", "conversions"):
        values = getattr(result, attr, None)
        if isinstance(values, dict):
            values = list(values.values())
        if isinstance(values, (tuple, list)):
            for value in values:
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
            return fn(*args, **kwargs)
        with operation(fn.__qualname__) as run:
            result = fn(*args, **kwargs)
            result_diagnostics(result)
            run.status = outcome(result)
            run.statistics.update(statistics(result))
            return result
    return wrapped
