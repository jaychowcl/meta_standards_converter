# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conversion-owned connection pools, without cross-request cookie state."""
from contextlib import contextmanager
from contextvars import ContextVar
from http.cookiejar import DefaultCookiePolicy
import threading
import time
from urllib3.connectionpool import HTTPConnectionPool, HTTPSConnectionPool
from urllib3.exceptions import EmptyPoolError
from urllib3.util import Timeout
from ..runtime_contracts import get_resource_profile

import requests

_active_pool = ContextVar('msc_http_pool', default=None)
_connect_wait = ContextVar('msc_pool_connect_wait', default=None)


class _AcquisitionStopped(EmptyPoolError):
    def __init__(self, pool, error):
        super().__init__(pool, 'Pool acquisition stopped before a connection was obtained')
        self.error = error


class _BoundedAcquisition:
    def urlopen(self, *args, **kwargs):
        try:
            return super().urlopen(*args, **kwargs)
        except _AcquisitionStopped as stopped:
            # urllib3 must see EmptyPoolError to avoid returning a phantom slot.
            # Restore cancellation/deadline/timeout only after its cleanup.
            raise stopped.error from None

    def _get_conn(self, timeout=None):
        from ..metadata.retrieval_budget import _active, _check_cancelled
        limit = _connect_wait.get()
        if limit is None:
            limit = get_resource_profile('standard').connect_timeout_seconds
        if timeout is not None:
            limit = min(limit, timeout)
        deadline = time.monotonic() + limit
        while True:
            try:
                _check_cancelled()
                budget = _active.get()
                remaining = deadline - time.monotonic()
                if budget is not None:
                    budget_remaining = budget.remaining()
                    if budget_remaining is not None:
                        remaining = min(remaining, budget_remaining)
                if remaining <= 0:
                    raise requests.exceptions.ConnectTimeout('Connection-pool acquisition deadline reached')
            except BaseException as error:
                raise _AcquisitionStopped(self, error) from None
            try:
                return super()._get_conn(timeout=min(.1, remaining))
            except EmptyPoolError:
                continue  # Waiting only; the outer requester still owns retries.


class _HTTPPool(_BoundedAcquisition, HTTPConnectionPool):
    pass


class _HTTPSPool(_BoundedAcquisition, HTTPSConnectionPool):
    pass


class _BoundedAdapter(requests.adapters.HTTPAdapter):
    def init_poolmanager(self, *args, **kwargs):
        super().init_poolmanager(*args, **kwargs)
        self.poolmanager.pool_classes_by_scheme = {'http': _HTTPPool, 'https': _HTTPSPool}

    def proxy_manager_for(self, *args, **kwargs):
        manager = super().proxy_manager_for(*args, **kwargs)
        # Preserve specialized proxy pools (such as SOCKS); standard proxies use
        # the same two urllib3 classes as PoolManager.
        if manager.pool_classes_by_scheme == {'http': HTTPConnectionPool, 'https': HTTPSConnectionPool}:
            manager.pool_classes_by_scheme = {'http': _HTTPPool, 'https': _HTTPSPool}
        return manager



class _NoStoredCookies(DefaultCookiePolicy):
    def set_ok(self, cookie, request):
        return False


class _Pool:
    def __init__(self):
        self.local = threading.local()
        self.lock = threading.Lock()
        self.sessions = []

    def get(self, url, **kwargs):
        session = getattr(self.local, 'session', None)
        if session is None:
            session = requests.Session()
            session.cookies.set_policy(_NoStoredCookies())
            # Each worker owns one session. The outer requester exclusively
            # owns retries, redirects, pacing and the request allowance.
            for scheme in ('http://', 'https://'):
                session.mount(scheme, _BoundedAdapter(
                    pool_connections=8, pool_maxsize=1, max_retries=0, pool_block=True))
            self.local.session = session
            with self.lock:
                self.sessions.append(session)
        timeout = kwargs.get('timeout')
        connect = timeout.connect_timeout if isinstance(timeout, Timeout) else (
            timeout[0] if isinstance(timeout, (list, tuple)) else timeout)
        token = _connect_wait.set(connect)
        try:
            return session.get(url, **kwargs)
        finally:
            _connect_wait.reset(token)

    def close(self):
        for session in self.sessions:
            session.close()
        self.sessions.clear()


@contextmanager
def transport_scope():
    """Nested converter calls share pools; separate conversions do not."""
    if _active_pool.get() is not None:
        yield
        return
    pool = _Pool()
    token = _active_pool.set(pool)
    try:
        yield
    finally:
        try:
            pool.close()
        finally:
            _active_pool.reset(token)


def get(url, **kwargs):
    pool = _active_pool.get()
    return pool.get(url, **kwargs) if pool is not None else requests.get(url, **kwargs)
