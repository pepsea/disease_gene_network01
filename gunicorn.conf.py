"""Gunicorn settings for running the app as a long-lived service.

The workload is IO-bound (every request fans out to Open Targets, SIGNOR,
STRING and BioGRID), so threads rather than processes carry the concurrency.
"""
import os


def _int_env(name: str, default: int) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


bind = f"0.0.0.0:{_int_env('PORT', 5005)}"

worker_class = "gthread"
# Each worker keeps its own SIGNOR index, so more workers means more copies of
# that bulk download. Scale threads first; add workers only for CPU headroom.
workers = _int_env("GUNICORN_WORKERS", 2)
threads = _int_env("GUNICORN_THREADS", 8)

# A cold SIGNOR download plus PPI lookups for a full gene batch can run long,
# so the request timeout sits well above the collectors' own read timeouts.
timeout = _int_env("GUNICORN_TIMEOUT", 300)
graceful_timeout = 30
keepalive = 5

# Recycle workers periodically so a leak in a long-running container cannot
# accumulate indefinitely.
max_requests = _int_env("GUNICORN_MAX_REQUESTS", 1000)
max_requests_jitter = 100

accesslog = "-"
errorlog = "-"
loglevel = os.environ.get("LOG_LEVEL", "info")
