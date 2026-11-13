"""Durable observation jobs. All write transitions are short SQLite transactions."""

from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
import time
import uuid
import zlib

from .dataio import canonical, digest, load_observation


LEASE_SECONDS = 600


def fingerprint():
    """Identify the code actually running, including local, uncommitted edits."""
    root = Path(__file__).parent
    modules = {p.name: digest(p.read_bytes()) for p in sorted(root.glob("*.py"))}
    return digest(canonical(modules).encode())


def connect(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=30, isolation_level=None)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA journal_mode=WAL")
    db.executescript("""
        CREATE TABLE IF NOT EXISTS observations (
            id TEXT PRIMARY KEY, manifest TEXT NOT NULL, rows BLOB NOT NULL
        );
        CREATE TABLE IF NOT EXISTS jobs (
            id TEXT PRIMARY KEY, observation TEXT NOT NULL REFERENCES observations(id),
            code TEXT NOT NULL, config TEXT NOT NULL, state TEXT NOT NULL,
            owner TEXT, lease REAL, created REAL NOT NULL, attempts INTEGER NOT NULL DEFAULT 0,
            calls INTEGER NOT NULL DEFAULT 0, error TEXT
        );
        CREATE TABLE IF NOT EXISTS stages (
            job TEXT NOT NULL REFERENCES jobs(id), stage TEXT NOT NULL, value TEXT NOT NULL,
            PRIMARY KEY(job,stage)
        );
        CREATE INDEX IF NOT EXISTS work_queue ON jobs(code,state,created);
    """)
    return db


@contextmanager
def transaction(db):
    db.execute("BEGIN IMMEDIATE")
    try:
        yield
    except BaseException:
        db.rollback()
        raise
    else:
        db.commit()


def enqueue(db, manifest_path, config, code):
    manifest, rows = load_observation(manifest_path)
    metadata = canonical(manifest)
    observation = digest(metadata.encode())
    settings = canonical(config)
    job_id = digest(canonical([observation, code, config]).encode())
    snapshot = zlib.compress(canonical(rows).encode())
    with transaction(db):
        db.execute("INSERT OR IGNORE INTO observations VALUES(?,?,?)", (observation, metadata, snapshot))
        db.execute("""INSERT OR IGNORE INTO jobs(id,observation,code,config,state,created)
                      VALUES(?,?,?,?,'queued',?)""", (job_id, observation, code, settings, time.time()))
    return job_id


def claim(db, code, only_job=None):
    now = time.time()
    with transaction(db):
        job = db.execute("""SELECT * FROM jobs WHERE code=? AND (? IS NULL OR id=?) AND
                (state='queued' OR (state='running' AND lease < ?))
                ORDER BY created,id LIMIT 1""", (code, only_job, only_job, now)).fetchone()
        if job is None:
            return None
        owner = uuid.uuid4().hex
        db.execute("""UPDATE jobs SET state='running', owner=?, lease=?, attempts=attempts+1,
                      error=NULL WHERE id=?""", (owner, now + LEASE_SECONDS, job["id"]))
    return {**dict(job), "owner": owner}


class LostLease(RuntimeError):
    """Another worker owns this job or its lease has expired."""


def _check_owner(db, job):
    row = db.execute("SELECT owner,lease,state FROM jobs WHERE id=?", (job["id"],)).fetchone()
    if row["state"] != "running" or row["owner"] != job["owner"] or row["lease"] <= time.time():
        raise LostLease("Worker lease lost")
    db.execute("UPDATE jobs SET lease=? WHERE id=?", (time.time() + LEASE_SECONDS, job["id"]))


def read_stages(db, job_id):
    return {r["stage"]: json.loads(r["value"]) for r in
            db.execute("SELECT stage,value FROM stages WHERE job=?", (job_id,))}


def checkpoint(db, job, stage, value):
    encoded = canonical(value)
    with transaction(db):
        _check_owner(db, job)
        db.execute("INSERT INTO stages VALUES(?,?,?)", (job["id"], stage, encoded))


def reserve_call(db, job):
    with transaction(db):
        _check_owner(db, job)
        db.execute("UPDATE jobs SET calls=calls+1 WHERE id=?", (job["id"],))


def finish(db, job, state, error=None):
    with transaction(db):
        _check_owner(db, job)
        db.execute("UPDATE jobs SET state=?,owner=NULL,lease=NULL,error=? WHERE id=?",
                   (state, error, job["id"]))


def retry(db, job_id, code):
    with transaction(db):
        cursor = db.execute("""UPDATE jobs SET state='queued',error=NULL WHERE id=?
                               AND code=? AND state='failed'""", (job_id, code))
        if not cursor.rowcount:
            raise ValueError("Retry requires a failed job from this code revision")


def source(db, job):
    row = db.execute("SELECT manifest,rows FROM observations WHERE id=?", (job["observation"],)).fetchone()
    return json.loads(row["manifest"]), json.loads(zlib.decompress(row["rows"]))


def summary(db):
    return [dict(row) for row in db.execute("""
        SELECT jobs.id, json_extract(observations.manifest,'$.observation_id') AS observation,
               jobs.state,jobs.attempts,jobs.calls,jobs.error,jobs.code
        FROM jobs JOIN observations ON observations.id=jobs.observation ORDER BY jobs.created,jobs.id
    """)]


def evidence(db, job_id):
    job = db.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
    if job is None:
        raise ValueError("Unknown job ID")
    manifest = db.execute("SELECT manifest FROM observations WHERE id=?", (job["observation"],)).fetchone()[0]
    return {"schema_version": 1, "job_id": job_id, "state": job["state"], "execution_error": job["error"],
            "software_sha256": job["code"], "observation": json.loads(manifest),
            "model_config": json.loads(job["config"]), "provider_calls_reserved": job["calls"],
            "stages": read_stages(db, job_id), "claim_status": "screening_only"}
