"""Minimal JWT auth for the curator area. Passwords use scrypt from the standard library."""
import base64
import hashlib
import hmac
import os
import secrets
import time
from typing import Callable, Optional

import jwt
from fastapi import Depends, HTTPException
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sqlalchemy import select
from sqlalchemy.orm import Session

from .database import get_session
from .models import User, UserRole

TOKEN_TTL_SECONDS = int(os.getenv("SEPTATYPER_TOKEN_TTL", str(8 * 3600)))
_bearer = HTTPBearer(auto_error=False)

_N, _R, _P = 2 ** 14, 8, 1


def hash_password(password: str) -> str:
    salt = secrets.token_bytes(16)
    dk = hashlib.scrypt(password.encode(), salt=salt, n=_N, r=_R, p=_P, dklen=32)
    return "scrypt${}${}${}${}${}".format(_N, _R, _P, base64.b64encode(salt).decode(), base64.b64encode(dk).decode())


def verify_password(password: str, stored: Optional[str]) -> bool:
    """Always does the same amount of work, even for unknown users (stored is None)."""
    if not stored or not stored.startswith("scrypt$"):
        hashlib.scrypt(password.encode(), salt=b"0" * 16, n=_N, r=_R, p=_P, dklen=32)
        return False
    _, n, r, p, salt_b64, dk_b64 = stored.split("$")
    dk = hashlib.scrypt(password.encode(), salt=base64.b64decode(salt_b64), n=int(n), r=int(r), p=int(p), dklen=32)
    return hmac.compare_digest(dk, base64.b64decode(dk_b64))


def _secret() -> str:
    secret = os.getenv("SEPTATYPER_JWT_SECRET")
    if not secret or len(secret) < 32:
        raise HTTPException(503, "Authentication is not configured (set SEPTATYPER_JWT_SECRET to 32+ random characters)")
    return secret


def create_token(user: User) -> str:
    now = int(time.time())
    return jwt.encode(
        {"sub": str(user.user_id), "role": user.role.value, "iat": now, "exp": now + TOKEN_TTL_SECONDS},
        _secret(), algorithm="HS256",
    )


def current_user(
    creds: Optional[HTTPAuthorizationCredentials] = Depends(_bearer),
    db: Session = Depends(get_session),
) -> User:
    if creds is None:
        raise HTTPException(401, "Authentication required", headers={"WWW-Authenticate": "Bearer"})
    try:
        payload = jwt.decode(creds.credentials, _secret(), algorithms=["HS256"])
    except jwt.PyJWTError as exc:
        raise HTTPException(401, "Invalid or expired token", headers={"WWW-Authenticate": "Bearer"}) from exc
    user = db.get(User, int(payload["sub"]))
    if user is None:  # role is re-read from the database, never trusted from the token
        raise HTTPException(401, "Unknown user")
    return user


def require_roles(*roles: UserRole) -> Callable[..., User]:
    def dep(user: User = Depends(current_user)) -> User:
        if user.role not in roles:
            raise HTTPException(403, "You do not have permission to do this")
        return user
    return dep


require_curator = require_roles(UserRole.CURATOR, UserRole.ADMIN)


def get_user_by_email(db: Session, email: str) -> Optional[User]:
    return db.scalar(select(User).where(User.email == email.strip().lower()))
