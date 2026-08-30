"""Shared HTTP session for the external data collectors.

A single pooled ``requests.Session`` is used by every collector so that the
thread pool in ``app.py`` reuses connections instead of opening a new socket
per gene, and so that transient upstream failures are retried once or twice
before the collector gives up.
"""
from __future__ import annotations

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

USER_AGENT = "nw_overlap_app/1.0 (+https://github.com/pepsea/disease_gene_network01)"


def build_session(
    total_retries: int = 2,
    backoff_factor: float = 0.5,
    pool_maxsize: int = 20,
) -> requests.Session:
    """Return a session that retries idempotent failures with backoff."""
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    retry = Retry(
        total=total_retries,
        connect=total_retries,
        read=total_retries,
        status=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "POST"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=pool_maxsize, pool_maxsize=pool_maxsize)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


SESSION = build_session()
