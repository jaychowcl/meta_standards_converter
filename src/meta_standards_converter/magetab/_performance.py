"""Private export measurements; scientific payloads never carry timing state."""
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import wraps
import sys
import time

from ..conversion_logging import stage


@dataclass
class _Frame:
    span: object
    child_wall: float = 0.0
    child_cpu: float = 0.0


_active = ContextVar("msc_magetab_phase", default=None)


def _peak_rss():
    try:
        import resource
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value if sys.platform == "darwin" else value * 1024)
    except (ImportError, AttributeError, OSError):
        return None


@contextmanager
def phase(name, **counts):
    """Report inclusive work and work excluding nested export phases.

    CPU and RSS describe this process. RSS is a lifetime high-water mark, not
    memory allocated by this phase; unrelated existing logging spans are not
    subtracted from the explicitly named exclusive-instrumented measurements.
    """
    parent = _active.get()
    with stage(name, **counts) as span:
        frame = _Frame(span)
        token = _active.set(frame)
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        initial_peak = _peak_rss()
        try:
            yield span
        finally:
            wall = max(0.0, time.perf_counter() - start_wall)
            cpu = max(0.0, time.process_time() - start_cpu)
            span.update(
                wall_seconds=wall,
                process_cpu_seconds=cpu,
                exclusive_instrumented_wall_seconds=max(0.0, wall - frame.child_wall),
                exclusive_instrumented_cpu_seconds=max(0.0, cpu - frame.child_cpu),
                peak_rss_bytes=_peak_rss(),
                initial_peak_rss_bytes=initial_peak,
                peak_rss_scope="process_lifetime_high_water",
            )
            _active.reset(token)
            if parent is not None:
                parent.child_wall += wall
                parent.child_cpu += cpu


def measured_phase(name):
    def decorate(function):
        @wraps(function)
        def measured(*args, **kwargs):
            with phase(name):
                return function(*args, **kwargs)
        return measured
    return decorate


def work(**counts):
    """Add cheap semantic work counts to the current export phase."""
    frame = _active.get()
    if frame is not None:
        frame.span.update(**{
            key: frame.span.statistics.get(key, 0) + value
            for key, value in counts.items()
        })
