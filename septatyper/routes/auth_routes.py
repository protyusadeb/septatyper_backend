import os

from fastapi import APIRouter, Depends, HTTPException, Request
from pydantic import BaseModel
from sqlalchemy.orm import Session

from ..auth import (TOKEN_TTL_SECONDS, create_token, current_user, get_user_by_email, verify_password)
from ..database import get_session
from ..models import User
from ..ratelimit import RateLimiter, enforce

router = APIRouter(prefix="/auth", tags=["auth"])
_login_limiter = RateLimiter(int(os.getenv("SEPTATYPER_LOGIN_LIMIT", "10")), 60)
_TRUST_PROXY = os.getenv("SEPTATYPER_TRUST_PROXY") == "1"


class LoginRequest(BaseModel):
    email: str
    password: str


class TokenOut(BaseModel):
    access_token: str
    token_type: str = "bearer"
    role: str
    expires_in: int


class MeOut(BaseModel):
    email: str
    role: str


@router.post("/login", response_model=TokenOut)
def login(body: LoginRequest, request: Request, db: Session = Depends(get_session)):
    enforce(_login_limiter, request, _TRUST_PROXY)
    user = get_user_by_email(db, body.email)
    ok = verify_password(body.password, user.hashed_password if user else None)
    if not user or not ok:
        raise HTTPException(401, "Incorrect email or password")
    return TokenOut(access_token=create_token(user), role=user.role.value, expires_in=TOKEN_TTL_SECONDS)


@router.get("/me", response_model=MeOut)
def me(user: User = Depends(current_user)):
    return MeOut(email=user.email, role=user.role.value)
