"""Process-wide prediction engine and concurrency cap, shared by /predict and the submission worker."""
import os
import threading
from functools import lru_cache
from typing import Optional

from .prediction_engine import SeptaTyperEngine

_slots = threading.BoundedSemaphore(int(os.getenv("SEPTATYPER_MAX_CONCURRENT", "2")))


class EngineBusy(RuntimeError):
    pass


@lru_cache(maxsize=1)
def get_prediction_engine() -> SeptaTyperEngine:
    return SeptaTyperEngine(
        reference_dir=os.getenv("SEPTATYPER_REF_DIR", "references"),
        blastn=os.getenv("BLASTN_PATH", "blastn"),
        makeblastdb=os.getenv("MAKEBLASTDB_PATH", "makeblastdb"),
        threads=int(os.getenv("SEPTATYPER_THREADS", "2")),
    )


def run_engine(text: str, species: Optional[str], wait: float = 20.0) -> dict:
    """Run one prediction while holding a concurrency slot. Raises EngineBusy if none frees up."""
    engine = get_prediction_engine()
    if not _slots.acquire(timeout=wait):
        raise EngineBusy("The server is busy")
    try:
        return engine.predict(text, species=species)
    finally:
        _slots.release()
