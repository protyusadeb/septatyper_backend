"""Engine and session handling.

Set DATABASE_URL, e.g.
    postgresql+psycopg2://septatyper:secret@localhost:5432/septatyper
"""
from __future__ import annotations

import os
from functools import lru_cache
from typing import Iterator

from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker

DEFAULT_URL = "postgresql+psycopg2://septatyper:septatyper@localhost:5432/septatyper"

SessionLocal = sessionmaker(autoflush=False, expire_on_commit=False, future=True)


def get_database_url() -> str:
    return os.getenv("DATABASE_URL", DEFAULT_URL)


@lru_cache(maxsize=4)
def get_engine(url: str | None = None) -> Engine:
    target = url or get_database_url()
    kwargs = {"connect_args": {"check_same_thread": False}} if target.startswith("sqlite") else {}
    engine = create_engine(target, pool_pre_ping=True, future=True, **kwargs)
    SessionLocal.configure(bind=engine)
    return engine


def init_db(url: str | None = None) -> Engine:
    """Create all tables (use Alembic for later migrations)."""
    from .models import Base

    engine = get_engine(url)
    Base.metadata.create_all(engine)
    return engine


def get_session() -> Iterator[Session]:
    """FastAPI dependency: `db: Session = Depends(get_session)`."""
    get_engine()
    session = SessionLocal()
    try:
        yield session
    finally:
        session.close()
