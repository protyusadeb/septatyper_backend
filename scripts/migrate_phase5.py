"""Upgrade a Phase 2 database for Module C.

    python -m scripts.migrate_phase5            (uses $DATABASE_URL)

Adds the new nullable columns on `submissions` and the helper indexes. Safe to re-run.
For later schema changes, move to Alembic.
"""
import sys

from sqlalchemy import inspect, text

from septatyper.database import get_engine
from septatyper.models import Base, Submission


def main() -> int:
    engine = get_engine()
    Base.metadata.create_all(engine)  # creates any missing tables; never alters existing ones
    insp = inspect(engine)
    existing = {c["name"] for c in insp.get_columns("submissions")}
    added = []
    with engine.begin() as conn:
        for col in Submission.__table__.columns:
            if col.name in existing:
                continue
            ddl = col.type.compile(dialect=engine.dialect)
            conn.execute(text(f'ALTER TABLE submissions ADD COLUMN "{col.name}" {ddl}'))
            added.append(col.name)
        for col in ("file_sha256", "bv_brc_id_supplied", "analysis_state"):
            conn.execute(text(f'CREATE INDEX IF NOT EXISTS ix_submissions_{col} ON submissions ("{col}")'))
    print("Added columns:", ", ".join(added) if added else "none (already up to date)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
