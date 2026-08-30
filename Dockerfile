# syntax=docker/dockerfile:1

FROM python:3.11-slim

# Faster, quieter, log-friendly Python in a container.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=5005

WORKDIR /app

# Dependencies first so application edits do not invalidate the install layer.
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app.py gunicorn.conf.py ./
COPY collectors/ ./collectors/
COPY templates/ ./templates/
COPY nw_overlap.py ./

# Run as an unprivileged user.
RUN useradd --create-home --uid 10001 appuser && chown -R appuser:appuser /app
USER appuser

EXPOSE 5005

# python:3.11-slim ships no curl or wget, so the probe uses the stdlib.
HEALTHCHECK --interval=30s --timeout=10s --start-period=20s --retries=3 \
  CMD python -c "import os,sys,urllib.request; \
url='http://127.0.0.1:%s/healthz' % os.environ.get('PORT','5005'); \
sys.exit(0 if urllib.request.urlopen(url, timeout=5).status == 200 else 1)"

CMD ["gunicorn", "-c", "gunicorn.conf.py", "app:app"]
