# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""CLI-owned log sinks and bounded, stderr-only progress rendering."""
from collections import Counter
from contextvars import ContextVar
from datetime import datetime, timezone
from functools import wraps
import json
import logging
import shutil
import sys
import threading
import time

from meta_standards_converter.operational_events import redact
from meta_standards_converter.conversion_logging import current, _labels, _span

_session = ContextVar("msc_cli_logging", default=None)
_cli_lock = threading.RLock()


def event_for(record):
    if hasattr(record, "msc_event"):
        return record.msc_event
    run = current()
    span = _span.get()
    return dict(schema_version="1.0", timestamp=datetime.fromtimestamp(record.created, timezone.utc).isoformat(),
        repository="meta_standards_converter", component=record.name, event_type="log",
        level=record.levelname, run_id=run.id if run else None, stage=span.name if span else None,
        status=None, duration_seconds=None,
        attributes=redact(dict(_labels.get()) | {"message": record.getMessage()}))


class EventFormatter(logging.Formatter):
    def __init__(self, *, jsonl=False, timestamps=False):
        super().__init__()
        self.jsonl = jsonl
        self.timestamps = timestamps

    def format(self, record):
        event = event_for(record)
        if self.jsonl:
            return json.dumps(event, ensure_ascii=True, separators=(",", ":"))
        attributes = event["attributes"]
        prefix = (event["timestamp"] + " ") if self.timestamps else ""
        prefix += event["level"] + " " + event["component"] + ": "
        if event["event_type"] == "log":
            return prefix + attributes["message"]
        bits = [event["event_type"]]
        for key in ("run_id", "stage", "status"):
            if event.get(key) is not None:
                bits.append(f"{key}={event[key]}")
        for key, value in attributes.items():
            bits.append(f"{key}={json.dumps(value, ensure_ascii=True) if isinstance(value, (dict, list)) else value}")
        if event["duration_seconds"] is not None:
            bits.append(f"elapsed_s={event['duration_seconds']:.3f}")
        return prefix + " ".join(bits)


class ProgressRenderer:
    def __init__(self, stream, *, mode="auto", quiet=False, clock=time.monotonic,
                 width=None, animate=True, total=None):
        self.stream = stream
        self.clock = clock
        self.width = width or (lambda: shutil.get_terminal_size((80, 24)).columns)
        self.mode = "off" if quiet else mode
        self.live = self.mode == "auto" and bool(getattr(stream, "isatty", lambda: False)())
        self.animate = animate
        self.total = total
        self.finished = 0
        self.states = Counter()
        self.last = None
        self.last_draw = float("-inf")
        self.line_length = 0
        self.closed = False
        self.lock = threading.RLock()
        self.stop = threading.Event()
        self.worker = None
        self.items_in_run = 0
        self.activity_started = clock()

    def update(self, event):
        kind = event["event_type"]
        attrs = event["attributes"]
        if kind not in {"conversion_started", "conversion_finished", "batch_planned", "item_started", "item_finished", "stage_started", "stage_finished", "progress", "sample_started", "sample_finished", "study_discovered", "preparation_outcome"}:
            return
        with self.lock:
            if kind == "conversion_started":
                self.closed = False
                self.items_in_run = 0
            if kind == "batch_planned":
                self.total = attrs.get("inputs")
            if kind == "item_finished":
                self.finished += 1
                self.items_in_run += 1
                self.states[event["status"]] += 1
            if kind == "conversion_finished" and not self.items_in_run:
                self.finished += 1
                self.states[event["status"]] += 1
            if kind in {"stage_started", "item_started"}:
                self.activity_started = self.clock()
            self.last = event
            if self.mode == "plain" and event["level"] != "DEBUG" and kind in {"conversion_started", "conversion_finished", "stage_started", "stage_finished", "item_finished", "progress"}:
                marker = {"running":"RUN", "completed":"OK", "complete":"OK", "partial":"PARTIAL",
                          "failed":"FAIL", "skipped":"SKIP"}.get(event["status"], "RUN")
                description = self._description(event)
                if kind in {"item_finished", "conversion_finished"}:
                    finished = f"{self.finished}/{self.total}" if self.total is not None else str(self.finished)
                    counts = " ".join(f"{key}={self.states[key]}" for key in ("complete", "partial", "failed", "skipped"))
                    description += f" | {finished} inputs finished | {counts}"
                self.stream.write(f"[{marker}] {description}\n")
                self.stream.flush()
            if self.live and self.animate and self.worker is None and kind != "conversion_finished":
                self.stop.clear()
                self.worker = threading.Thread(target=self._animate, name="msc-terminal-progress", daemon=True)
                self.worker.start()
        if kind == "conversion_finished":
            self.close()
        else:
            self.draw()

    def _description(self, event):
        attrs = event["attributes"]
        bits = [str(attrs.get("dataset_id") or attrs.get("item_id") or attrs.get("source") or event["component"]),
                str(event.get("stage") or event["event_type"])]
        if attrs.get("total") is not None:
            bits.append(f"{attrs.get('completed', 0)}/{attrs['total']} {attrs.get('unit', 'units')}")
        elif "completed" in attrs:
            bits.append(f"{attrs['completed']} processed")
        if event.get("status") == "running":
            bits.append(f"{self.clock() - self.activity_started:.1f}s")
        elif event.get("duration_seconds") is not None:
            bits.append(f"{event['duration_seconds']:.1f}s")
        return " | ".join(bits)

    def _animate(self):
        while not self.stop.wait(.25):
            self.draw()

    def draw(self):
        with self.lock:
            now = self.clock()
            if not self.live or self.closed or self.last is None or now - self.last_draw < .25:
                return
            self.last_draw = now
            if self.total is not None and self.total > 0:
                done = min(10, int(10 * self.finished / self.total))
                prefix = f"[{'#'*done}{'-'*(10-done)}] {self.finished}/{self.total} inputs finished"
            else:
                prefix = f"[{ '|/-\\'[int(now*4)%4] }] {self.finished} inputs finished"
            counts = " ".join(f"{key}={self.states[key]}" for key in ("complete", "partial", "failed", "skipped") if self.states[key])
            line = prefix + (" " + counts if counts else "") + " | " + self._description(self.last)
            line = line[:max(1, self.width()-1)]
            self._clear()
            self.stream.write("\r" + line)
            self.line_length = len(line)
            self.stream.flush()

    def _clear(self):
        if self.line_length:
            self.stream.write("\r" + " " * self.line_length + "\r")
            self.line_length = 0

    def message(self, text):
        with self.lock:
            self._clear()
            self.stream.write(text + "\n")
            self.stream.flush()
        self.draw()

    def close(self):
        self.stop.set()
        worker = self.worker
        if worker is not None and worker is not threading.current_thread():
            worker.join(timeout=2)
        with self.lock:
            self._clear()
            self.stream.flush()
            self.worker = None
            self.closed = True


class ConsoleHandler(logging.Handler):
    def __init__(self, stream, *, threshold, mode="auto", quiet=False, total=None):
        super().__init__(logging.DEBUG)
        self.threshold = threshold
        self.progress = ProgressRenderer(stream, mode=mode, quiet=quiet, total=total)
        self.setFormatter(EventFormatter())
        self.repeats = {}
        self.results = Counter()
        self._msc_owned = True

    def emit(self, record):
        event = event_for(record)
        self.progress.update(event)
        if event["event_type"] == "conversion_finished":
            self.results[event["status"]] += 1
            self.flush_repeats(event["run_id"])
        if record.levelno < self.threshold:
            return
        if self.progress.mode == "plain" and record.levelno >= logging.INFO and event["event_type"] in {"conversion_started", "conversion_finished", "stage_started", "stage_finished", "progress"}:
            return
        if record.levelno == logging.WARNING and event["event_type"] == "diagnostic":
            attrs = event["attributes"]
            key = (event["run_id"], attrs.get("item_id"), event["stage"], attrs.get("code"))
            self.repeats[key] = self.repeats.get(key, 0) + 1
            if self.repeats[key] > 1:
                return
        self.progress.message(self.format(record))

    def flush_repeats(self, run_id=None):
        for key, count in list(self.repeats.items()):
            if run_id is None or key[0] == run_id:
                if count > 1:
                    self.progress.message(f"WARNING item={key[1]} stage={key[2]} code={key[3]}: {count} occurrences (details in log/report)")
                del self.repeats[key]

    def close(self):
        self.flush_repeats()
        self.progress.close()
        super().close()


def register_configuration(args, logger, handlers, previous):
    session = _session.get()
    if session is not None:
        session.update(args=args, logger=logger, handlers=handlers, previous=previous)


def cli_entrypoint(fn):
    """Release CLI-owned sinks even when argparse, conversion or output fails."""
    @wraps(fn)
    def wrapped(*args, **kwargs):
        _cli_lock.acquire()
        session = {}
        token = _session.set(session)
        code = None
        try:
            code = fn(*args, **kwargs)
            configured = session.get("args")
            command = fn.__module__.rsplit(".", 1)[-1]
            if configured is not None and not configured.quiet and command in {"geo2ae", "geo2json", "json2ae", "ae2json", "sra2json", "ena2json"}:
                console = next((h for h in session["handlers"] if isinstance(h, ConsoleHandler)), None)
                states = dict(console.results) if console else {}
                human_summary("complete" if code == 0 else "partial/failed", counts=states)
            return code
        except Exception as error:
            from meta_standards_converter.diagnostics import RequestError
            from .common import record_safe_cli_error
            if isinstance(error, RequestError):
                sys.stderr.write('error: ' + str(error) + '\n')
                raise SystemExit(2) from None
            record_safe_cli_error(logging.getLogger(fn.__module__), error,
                                  location=fn.__module__.rsplit('.', 1)[-1], stage='command')
            return 1
        finally:
            try:
                for handler in session.get("handlers", ()):
                    session["logger"].removeHandler(handler)
                    handler.close()
                if "previous" in session:
                    session["logger"].setLevel(session["previous"][0])
                    session["logger"].propagate = session["previous"][1]
            finally:
                _session.reset(token)
                _cli_lock.release()
    return wrapped


def human_summary(status, *, counts=None, artifacts=()):
    """Render compact human results to stderr, leaving stdout for data."""
    print("Conversion: " + status, file=sys.stderr)
    if counts:
        print(" ".join(f"{key}={value}" for key, value in counts.items()), file=sys.stderr)
    for identity, name, path in artifacts:
        print(f"  {identity} {name}: {path}", file=sys.stderr)
