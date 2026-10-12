"""Shared FASTA upload reading (size caps, gzip, encoding)."""
import gzip
import io

from fastapi import HTTPException, UploadFile


def read_fasta_upload(upload: UploadFile, max_bytes: int) -> str:
    data = upload.file.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise HTTPException(413, f"File exceeds the {max_bytes // (1024 * 1024)} MB limit")
    if data[:2] == b"\x1f\x8b":  # gzip; cap the decompressed size to defuse zip bombs
        try:
            data = gzip.GzipFile(fileobj=io.BytesIO(data)).read(max_bytes + 1)
        except (OSError, EOFError) as exc:
            raise HTTPException(400, "Invalid gzip file") from exc
        if len(data) > max_bytes:
            raise HTTPException(413, "Decompressed file exceeds the size limit")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(400, "File is not UTF-8 text; upload a plain or gzipped FASTA") from exc


def looks_like_nucleotide_fasta(text: str) -> bool:
    """Cheap pre-check so obviously wrong files are rejected before any background work."""
    t = text.lstrip()
    if not t.startswith(">"):
        return False
    body = "".join(line.strip() for line in t.splitlines() if line and not line.startswith(">"))[:20000]
    if not body:
        return False
    ok = sum(1 for ch in body.upper() if ch in "ACGTRYSWKMBDHVN")
    return ok / len(body) > 0.98
