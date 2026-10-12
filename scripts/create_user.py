"""Create or update a curator/admin account.

    python -m scripts.create_user curator@lab.org curator
    python -m scripts.create_user admin@lab.org admin
"""
import getpass
import sys

from sqlalchemy import select

from septatyper.auth import hash_password
from septatyper.database import SessionLocal, get_engine, init_db
from septatyper.models import User, UserRole


def main(argv) -> int:
    if len(argv) != 3 or argv[2] not in {r.value for r in UserRole}:
        print("usage: python -m scripts.create_user EMAIL {public|submitter|curator|admin}")
        return 2
    email, role = argv[1].strip().lower(), UserRole(argv[2])
    pw = getpass.getpass("Password (min 12 characters): ")
    if len(pw) < 12 or pw != getpass.getpass("Repeat password: "):
        print("Passwords must match and be at least 12 characters")
        return 1
    init_db()
    get_engine()
    with SessionLocal() as db:
        user = db.scalar(select(User).where(User.email == email))
        if user is None:
            user = User(email=email)
            db.add(user)
        user.hashed_password, user.role = hash_password(pw), role
        db.commit()
    print(f"Saved {email} as {role.value}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
