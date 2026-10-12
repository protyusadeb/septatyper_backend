# SeptaTyper API: FastAPI + BLAST+ (the prediction engine shells out to blastn/makeblastdb)
FROM python:3.12-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONPATH=/app

RUN apt-get update \
 && apt-get install -y --no-install-recommends ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

# non-root user; /data/submissions is the persistent upload volume
RUN useradd --system --uid 10001 --home-dir /app app \
 && mkdir -p /app /data/submissions \
 && chown -R app:app /app /data

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY --chown=app:app septatyper ./septatyper
COPY --chown=app:app scripts ./scripts
COPY --chown=app:app --chmod=755 docker/entrypoint.sh /entrypoint.sh

COPY --chown=app:app setup_references.py
RUN python setup_references.py

USER app
EXPOSE 8000
ENTRYPOINT ["/entrypoint.sh"]
# One worker on purpose: the rate limiter and the BLAST concurrency cap live in process memory.
CMD ["uvicorn", "septatyper.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1", "--proxy-headers", "--forwarded-allow-ips", "*"]
