#!/bin/sh
# Builds DATABASE_URL from POSTGRES_* (URL-encoding the password) and refuses to start the
# server with a weak JWT secret. Applies to `docker compose run` commands too.
set -eu

if [ -z "${DATABASE_URL:-}" ]; then
  : "${POSTGRES_USER:?POSTGRES_USER or DATABASE_URL must be set}"
  : "${POSTGRES_PASSWORD:?POSTGRES_PASSWORD or DATABASE_URL must be set}"
  : "${POSTGRES_DB:?POSTGRES_DB or DATABASE_URL must be set}"
  PW="$(python3 -c 'import os, urllib.parse as u; print(u.quote(os.environ["POSTGRES_PASSWORD"], safe=""))')"
  export DATABASE_URL="postgresql+psycopg2://${POSTGRES_USER}:${PW}@${POSTGRES_HOST:-db}:5432/${POSTGRES_DB}"
fi

if [ "${1:-}" = "uvicorn" ]; then
  if [ "${#SEPTATYPER_JWT_SECRET}" -lt 32 ]; then
    echo "SEPTATYPER_JWT_SECRET must be at least 32 characters. Generate one with: openssl rand -hex 32" >&2
    exit 1
  fi
fi

exec "$@"
