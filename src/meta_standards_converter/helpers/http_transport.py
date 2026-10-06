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

import requests

_active_pool = ContextVar('msc_http_pool', default=None)


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
                session.mount(scheme, requests.adapters.HTTPAdapter(
                    pool_connections=8, pool_maxsize=1, max_retries=0, pool_block=True))
            self.local.session = session
            with self.lock:
                self.sessions.append(session)
        return session.get(url, **kwargs)

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
