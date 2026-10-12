"""BV-BRC (formerly PATRIC) API client, shared by the ingestion script and the submission portal."""
from __future__ import annotations

import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .models import BvbrcPrediction

log = logging.getLogger("septatyper.bvbrc")

class BVBRCClient:
    """Polite batched client for the BV-BRC (formerly PATRIC) Solr API."""

    def __init__(self, base: str = "https://www.bv-brc.org/api", sleep: float = 0.5,
                 batch_size: int = 50, timeout: int = 90, retries: int = 5):
        self.base, self.sleep, self.batch_size, self.timeout = base.rstrip("/"), sleep, batch_size, timeout
        self.session = requests.Session()
        retry = Retry(total=retries, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504),
                      allowed_methods=frozenset(["GET", "POST"]), respect_retry_after_header=True)
        self.session.mount("https://", HTTPAdapter(max_retries=retry))
        self.session.headers.update({"User-Agent": "SeptaTyper-ingest/1.0 (academic database build)"})

    def _rql(self, collection: str, rql: str) -> List[Dict[str, Any]]:
        resp = self.session.post(
            f"{self.base}/{collection}/",
            data=rql.encode("utf-8"),
            headers={"Content-Type": "application/rqlquery+x-www-form-urlencoded",
                     "Accept": "application/json"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        return resp.json()

    @staticmethod
    def _id_list(ids: Iterable[str]) -> str:
        # IDs MUST be quoted: unquoted "562.7710" is parsed as the number 562.771.
        return ",".join(f'"{i}"' for i in ids)

    def fetch_batch(self, ids: List[str]) -> Dict[str, Dict[str, Any]]:
        idl = self._id_list(ids)
        genomes = self._rql("genome", f"in(genome_id,({idl}))&select(genome_id,genome_name)&limit({len(ids)})")
        time.sleep(self.sleep)
        limit = max(500, len(ids) * 20)
        amr = self._rql(
            "genome_amr",
            f'and(in(genome_id,({idl})),eq(antibiotic,"ciprofloxacin"))'
            "&select(genome_id,antibiotic,resistant_phenotype,evidence,computational_method,"
            "computational_method_version,laboratory_typing_method)"
            f"&limit({limit})",
        )
        if len(amr) >= limit:
            log.warning("genome_amr hit the %d-record limit; some records may be truncated", limit)
        out: Dict[str, Dict[str, Any]] = {i: {"genome_name": None, "amr": []} for i in ids}
        for g in genomes:
            if g.get("genome_id") in out:
                out[g["genome_id"]]["genome_name"] = g.get("genome_name")
        for a in amr:
            if a.get("genome_id") in out:
                out[a["genome_id"]]["amr"].append(a)
        time.sleep(self.sleep)
        return out

    def fetch_all(self, ids: List[str], cache_path: Optional[Path]) -> Dict[str, Dict[str, Any]]:
        cache: Dict[str, Dict[str, Any]] = {}
        if cache_path and cache_path.exists():
            cache = json.loads(cache_path.read_text())
            log.info("Loaded %d cached BV-BRC records from %s", len(cache), cache_path)
        todo = [i for i in ids if i not in cache]
        for n in range(0, len(todo), self.batch_size):
            chunk = todo[n:n + self.batch_size]
            try:
                result = self.fetch_batch(chunk)
            except requests.RequestException as exc:
                log.error("BV-BRC batch %d-%d failed: %s (will be retried on the next run)", n, n + len(chunk), exc)
                continue
            stamp = datetime.now(timezone.utc).isoformat()
            for gid, rec in result.items():
                rec["fetched_at"] = stamp
                cache[gid] = rec
            if cache_path:
                cache_path.write_text(json.dumps(cache))
            log.info("BV-BRC: %d / %d fetched", min(n + len(chunk), len(todo)), len(todo))
        return cache


def summarise_computational(amr_records: List[Dict[str, Any]]) -> Tuple[BvbrcPrediction, Optional[str], Optional[str]]:
    """Collapse BV-BRC ciprofloxacin records into one computational call.

    Only records with evidence == 'Computational Method' count. No such record
    -> No_call. Conflicting records -> No_call. Intermediate and similar
    values -> No_call (raw value is kept in bv_brc_prediction_raw).
    """
    comp = [r for r in amr_records if str(r.get("evidence", "")).strip().lower() == "computational method"]
    if not comp:
        return BvbrcPrediction.NO_CALL, None, None
    raw_values = sorted({str(r.get("resistant_phenotype", "")).strip() for r in comp})
    methods = sorted({f"{r.get('computational_method') or 'unknown'} {r.get('computational_method_version') or ''}".strip() for r in comp})
    mapped = {v.lower() for v in raw_values}
    method_str = "; ".join(methods)
    if mapped == {"resistant"}:
        return BvbrcPrediction.RESISTANT, raw_values[0], method_str
    if mapped == {"susceptible"}:
        return BvbrcPrediction.SUSCEPTIBLE, raw_values[0], method_str
    return BvbrcPrediction.NO_CALL, " | ".join(raw_values), method_str


