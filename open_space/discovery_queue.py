"""A local durable queue for independent archive measurements and follow-up decisions."""

from contextlib import closing
import json
from pathlib import Path
import random
import sqlite3
import threading
import time
import uuid

from .dataio import canonical, digest
from .ledger import transaction

LEASE = 180


def connect(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=40, isolation_level=None)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA busy_timeout=40000")
    return db


def initialize(path, project):
    with closing(connect(path)) as db:
        db.execute("PRAGMA journal_mode=WAL")
        db.executescript("""
        CREATE TABLE IF NOT EXISTS campaign (id INTEGER PRIMARY KEY CHECK(id=1), project TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS task(
          id TEXT PRIMARY KEY, kind TEXT NOT NULL, subject TEXT NOT NULL, payload TEXT NOT NULL,
          state TEXT NOT NULL, owner TEXT, expires REAL, created REAL NOT NULL,
          attempts INTEGER NOT NULL DEFAULT 0, calls INTEGER NOT NULL DEFAULT 0,
          result TEXT, error TEXT);
        CREATE TABLE IF NOT EXISTS dependency (
            task_id TEXT NOT NULL, required_id TEXT NOT NULL, PRIMARY KEY(task_id,required_id)
        );
        CREATE TABLE IF NOT EXISTS lineage (
            parent_id TEXT NOT NULL, child_id TEXT NOT NULL, PRIMARY KEY(parent_id,child_id)
        );
        CREATE INDEX IF NOT EXISTS lineage_child ON lineage(child_id);
        CREATE INDEX IF NOT EXISTS task_order ON task(state,created);
        """)
        with transaction(db):
            db.execute(
                "INSERT OR IGNORE INTO campaign VALUES (1,?)", (canonical(project),)
            )
            if config(db) != project:
                raise ValueError(
                    "This directory already contains a different frozen campaign"
                )


def config(db):
    row = db.execute("SELECT project FROM campaign WHERE id=1").fetchone()
    if not row:
        raise ValueError("Initialize a discovery campaign first")
    return json.loads(row[0])


def enqueue(db, kind, subject, payload):
    key = digest(canonical([kind, subject, payload]).encode())
    db.execute(
        """INSERT OR IGNORE INTO task(id,kind,subject,payload,state,created)
                  VALUES (?,?,?,?,'queued',?)""",
        (key, kind, subject, canonical(payload), time.time()),
    )
    if kind == "review":
        db.executemany(
            "INSERT OR IGNORE INTO dependency VALUES(?,?)",
            [(key, parent) for parent in payload["children"]],
        )
    return key


def claim(db, max_attempts):
    now = time.time()
    with transaction(db):
        # A crash on the last attempt must not leave a job running forever.
        db.execute(
            """UPDATE task SET state='failed',error='final worker lease expired',owner=NULL,expires=NULL
                      WHERE state='running' AND expires<=? AND attempts>=?""",
            (now, max_attempts),
        )
        db.execute(
            """UPDATE task SET state='failed',error='worker attempt ceiling reached'
                      WHERE state='queued' AND attempts>=?""",
            (max_attempts,),
        )
        row = db.execute(
            """SELECT * FROM task WHERE (state='queued' OR (state='running' AND expires<=?))
                            AND attempts<? AND NOT EXISTS (
                                SELECT 1 FROM dependency d JOIN task t ON d.required_id=t.id
                                WHERE d.task_id=task.id AND t.state NOT IN ('done','failed'))
                            ORDER BY created,id LIMIT 1""",
            (now, max_attempts),
        ).fetchone()
        if row is None:
            return None
        owner = uuid.uuid4().hex
        db.execute(
            "UPDATE task SET owner=?,state='running',expires=?,attempts=attempts+1,error=NULL WHERE id=?",
            (owner, now + LEASE, row["id"]),
        )
    return {
        **dict(row),
        "payload": json.loads(row["payload"]),
        "owner": owner,
        "attempts": row["attempts"] + 1,
    }


class LostLease(RuntimeError):
    pass


def check(db, job):
    row = db.execute(
        "SELECT owner,state,expires FROM task WHERE id=?", (job["id"],)
    ).fetchone()
    if (
        row["state"] != "running"
        or row["owner"] != job["owner"]
        or row["expires"] <= time.time()
    ):
        raise LostLease("A different worker owns this measurement")


class Heartbeat:
    def __init__(self, path, job):
        self.path, self.job = path, job
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        with closing(connect(self.path)) as db:
            while not self.stop.wait(LEASE / 4):
                cursor = db.execute(
                    """UPDATE task SET expires=? WHERE id=? AND state='running'
                                       AND owner=? AND expires>?""",
                    (
                        time.time() + LEASE,
                        self.job["id"],
                        self.job["owner"],
                        time.time(),
                    ),
                )
                if not cursor.rowcount:
                    return

    def __enter__(self):
        self.thread.start()

    def __exit__(self, *_):
        self.stop.set()
        self.thread.join()


def reserve_model(db, job, cap):
    with transaction(db):
        check(db, job)
        if db.execute("SELECT COALESCE(SUM(calls),0) FROM task").fetchone()[0] >= cap:
            return False
        db.execute("UPDATE task SET calls=calls+1 WHERE id=?", (job["id"],))
    return True


def complete(db, job, result, children):
    encoded = canonical(result)
    with transaction(db):
        check(db, job)
        for kind, subject, payload in children:
            child = enqueue(db, kind, subject, payload)
            db.execute("INSERT OR IGNORE INTO lineage VALUES(?,?)", (job["id"], child))
            if job["kind"] in ("pixels", "variable_pixels"):
                # A reviewer scheduled before image localization must also see any
                # subsequently requested measurement at the displaced position.
                db.execute(
                    """INSERT OR IGNORE INTO dependency(task_id,required_id)
                              SELECT task_id,? FROM dependency WHERE required_id=?""",
                    (child, job["id"]),
                )
        db.execute(
            """UPDATE task SET result=?,state='done',owner=NULL,expires=NULL,error=NULL WHERE id=?""",
            (encoded, job["id"]),
        )


def fail(db, job, exc, max_attempts, transient=False):
    with transaction(db):
        check(db, job)
        state = "queued" if transient and job["attempts"] < max_attempts else "failed"
        db.execute(
            "UPDATE task SET state=?,owner=NULL,expires=NULL,error=? WHERE id=?",
            (state, f"{type(exc).__name__}: {str(exc)[:400]}", job["id"]),
        )
    if state == "queued":
        time.sleep(min(10, 2 ** job["attempts"]) + random.random())


def retry(db, task_id):
    """Rerun failed measurements and invalidate critics derived from the old failure."""
    with transaction(db):
        changed = db.execute(
            "UPDATE task SET state='queued',attempts=0,error=NULL WHERE state='failed' AND id=?",
            (task_id,),
        )
        if not changed.rowcount:
            raise ValueError("Expected a failed task ID")
        # Removing ownership here prevents a stale running reviewer from committing.
        db.execute(
            """UPDATE task SET state='queued',attempts=0,result=NULL,error=NULL,owner=NULL,expires=NULL
                      WHERE id IN (SELECT task_id FROM dependency WHERE required_id=?)""",
            (task_id,),
        )
    return {"requeued": task_id}


def stats(db):
    return {
        "tasks": {
            r["state"]: r["n"]
            for r in db.execute("SELECT state,count(*) n FROM task GROUP BY state")
        },
        "by_kind": [
            dict(r)
            for r in db.execute(
                "SELECT kind,state,count(*) AS count FROM task GROUP BY kind,state ORDER BY kind,state"
            )
        ],
        "model_calls_reserved": db.execute(
            "SELECT COALESCE(SUM(calls),0) FROM task"
        ).fetchone()[0],
    }
