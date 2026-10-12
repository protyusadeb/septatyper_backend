"""FastAPI entry point:  uvicorn septatyper.main:app --reload"""
import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .routes import auth_routes, isolates, predict, stats, submissions

app = FastAPI(
    title="SeptaTyper API",
    version="1.0.0",
    description="Browse and query ciprofloxacin resistance septatypes in Enterobacteriaceae.",
)

origins = [o.strip() for o in os.getenv("FRONTEND_ORIGINS", "http://localhost:3000").split(",") if o.strip()]
app.add_middleware(CORSMiddleware, allow_origins=origins, allow_methods=["GET", "POST"], allow_headers=["*"])

app.include_router(isolates.router, prefix="/api/v1")
app.include_router(stats.router, prefix="/api/v1")
app.include_router(predict.router, prefix="/api/v1")
app.include_router(auth_routes.router, prefix="/api/v1")
app.include_router(submissions.router, prefix="/api/v1")


@app.get("/api/v1/health", tags=["meta"])
def health():
    return {"status": "ok"}
