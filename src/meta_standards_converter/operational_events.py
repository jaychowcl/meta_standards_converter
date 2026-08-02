# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operational event v1 logging, JSONL, redaction, and metrics export."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import re
import threading
from typing import Any, Mapping
from urllib.parse import urlsplit, urlunsplit

EVENT_SCHEMA_VERSION = "1.0"
_SENSITIVE_KEY = re.compile(r"(?:api[_-]?key|authorization|credential|password|secret|token)", re.I)
_LABEL = re.compile(r"[^a-zA-Z0-9_:.-]+")


def redact(value: Any, *, key: str = "") -> Any:
    if _SENSITIVE_KEY.search(key):
        return "[REDACTED]"
    if isinstance(value, Mapping):
        return {str(k): redact(v, key=str(k)) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [redact(item) for item in value]
    if isinstance(value, str) and key.lower() in {"url", "uri", "endpoint"}:
        parsed = urlsplit(value)
        if parsed.scheme and parsed.netloc:
            return urlunsplit((parsed.scheme, parsed.netloc, parsed.path, "", ""))
    return value


class OperationalEventEmitter:
    """Emit the cross-repository schema-1.0 operational envelope."""

    def __init__(self, *, repository: str, jsonl_path: str | Path | None = None, logger: logging.Logger | None = None) -> None:
        self.repository = repository
        self.jsonl_path = Path(jsonl_path) if jsonl_path is not None else None
        self.logger = logger or logging.getLogger(f"{repository}.operations")
        self._metrics = defaultdict(lambda: {"count": 0, "duration_seconds": 0.0})
        self._lock = threading.Lock()

    def emit(self, event_type: str, *, component: str, level: str = "INFO", run_id: str | None = None, stage: str | None = None, status: str | None = None, duration_seconds: float | None = None, attributes: Mapping[str, Any] | None = None) -> dict[str, Any]:
        if not event_type.strip() or not component.strip():
            raise ValueError("event_type and component must not be blank")
        event = {"schema_version": EVENT_SCHEMA_VERSION, "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"), "repository": self.repository, "component": component, "event_type": event_type, "level": level.upper(), "run_id": run_id, "stage": stage, "status": status, "duration_seconds": duration_seconds, "attributes": redact(dict(attributes or {}))}
        payload = json.dumps(event, sort_keys=True, separators=(",", ":"))
        with self._lock:
            metric = self._metrics[(event_type, status or "unknown")]
            metric["count"] += 1
            metric["duration_seconds"] += float(duration_seconds or 0)
            if self.jsonl_path is not None:
                self.jsonl_path.parent.mkdir(parents=True, exist_ok=True)
                with self.jsonl_path.open("a", encoding="utf-8") as stream:
                    stream.write(payload + "\n")
        self.logger.log(getattr(logging, level.upper(), logging.INFO), payload)
        return event

    def metrics(self) -> dict[str, Any]:
        with self._lock:
            events = {f"{event}|{status}": dict(values) for (event, status), values in sorted(self._metrics.items())}
        return {"schema_version": EVENT_SCHEMA_VERSION, "repository": self.repository, "events": events}

    def export_json(self, path: str | Path) -> dict[str, Any]:
        value = self.metrics()
        self._write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")
        return value

    def prometheus_text(self) -> str:
        lines = ["# HELP operational_events_total Operational event count.", "# TYPE operational_events_total counter", "# HELP operational_event_duration_seconds_sum Aggregate event duration.", "# TYPE operational_event_duration_seconds_sum counter"]
        with self._lock:
            items = sorted(self._metrics.items())
        for (event, status), values in items:
            labels = f'repository="{_LABEL.sub("_", self.repository)}",event_type="{_LABEL.sub("_", event)}",status="{_LABEL.sub("_", status)}"'
            lines.append(f"operational_events_total{{{labels}}} {values['count']}")
            lines.append(f"operational_event_duration_seconds_sum{{{labels}}} {values['duration_seconds']:.9g}")
        return "\n".join(lines) + "\n"

    def export_prometheus(self, path: str | Path) -> str:
        value = self.prometheus_text()
        self._write(path, value)
        return value

    @staticmethod
    def _write(path: str | Path, value: str) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        stage = destination.with_name(f".{destination.name}.{os.getpid()}.stage")
        stage.write_text(value, encoding="utf-8")
        os.replace(stage, destination)


__all__ = ["EVENT_SCHEMA_VERSION", "OperationalEventEmitter", "redact"]
