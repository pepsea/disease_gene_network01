"""On-disk cache shared by the PPI collectors.

The original project cached SIGNOR's bulk TSV and per-gene STRING results
under ``ppi_cache/``. The same scheme is kept here, but the location is
configurable so that a container can mount it as a volume and keep the cache
across restarts.
"""
from __future__ import annotations

import os
from pathlib import Path

CACHE_DIR = Path(os.environ.get("PPI_CACHE_DIR", Path(__file__).parent.parent / "ppi_cache"))


def cache_path(*parts: str) -> Path:
    """Return a path inside the cache directory, creating its parent."""
    path = CACHE_DIR.joinpath(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def cache_enabled() -> bool:
    """Disk caching can be turned off with ``PPI_CACHE_DISABLED=1``."""
    return os.environ.get("PPI_CACHE_DISABLED", "").strip().lower() not in {"1", "true", "yes", "on"}
