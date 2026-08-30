"""Test doubles for the ``requests`` layer."""
from __future__ import annotations

import json as _json
from typing import Any


class FakeResponse:
    def __init__(self, payload: Any = None, status: int = 200, text: str | None = None):
        self._payload = payload
        self.status_code = status
        self._text = text

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        if self._text is not None:
            return _json.loads(self._text)  # raises ValueError on bad JSON
        return self._payload


class FakeSession:
    """Records calls and replays queued or callable responses."""

    def __init__(self, get=None, post=None):
        self._get = get
        self._post = post
        self.get_calls: list[dict] = []
        self.post_calls: list[dict] = []

    def _resolve(self, handler, url, kwargs):
        if handler is None:
            raise AssertionError(f"unexpected request to {url}")
        if callable(handler):
            return handler(url, **kwargs)
        if isinstance(handler, Exception):
            raise handler
        return handler

    def get(self, url, **kwargs):
        self.get_calls.append({"url": url, **kwargs})
        return self._resolve(self._get, url, kwargs)

    def post(self, url, **kwargs):
        self.post_calls.append({"url": url, **kwargs})
        return self._resolve(self._post, url, kwargs)
