from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Iterable

from .models import GlossaryEntry, StyleProfile, utc_now


SCHEMA = """
CREATE TABLE IF NOT EXISTS translation_memory (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_text TEXT NOT NULL,
    target_text TEXT NOT NULL,
    source_lang TEXT NOT NULL,
    target_lang TEXT NOT NULL,
    domain TEXT NOT NULL DEFAULT 'General',
    status TEXT NOT NULL DEFAULT 'approved',
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE(source_text, source_lang, target_lang, domain)
);

CREATE TABLE IF NOT EXISTS glossary (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_term TEXT NOT NULL,
    target_term TEXT NOT NULL,
    source_lang TEXT NOT NULL,
    target_lang TEXT NOT NULL,
    domain TEXT NOT NULL DEFAULT 'General',
    case_sensitive INTEGER NOT NULL DEFAULT 0,
    notes TEXT NOT NULL DEFAULT '',
    created_at TEXT NOT NULL,
    UNIQUE(source_term, source_lang, target_lang, domain)
);

CREATE TABLE IF NOT EXISTS style_profiles (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE,
    tone TEXT NOT NULL,
    audience TEXT NOT NULL,
    rules TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS audit_log (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    event_type TEXT NOT NULL,
    source_text TEXT NOT NULL,
    target_text TEXT NOT NULL,
    actor TEXT NOT NULL DEFAULT 'local_user',
    metadata TEXT NOT NULL DEFAULT '',
    created_at TEXT NOT NULL
);
"""


class StudioStore:
    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.initialize()

    def connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def initialize(self) -> None:
        with self.connect() as conn:
            conn.executescript(SCHEMA)
            self.seed_profiles(conn)

    def seed_profiles(self, conn: sqlite3.Connection) -> None:
        defaults = [
            StyleProfile("Enterprise Formal", "Formal", "Enterprise", "Use precise terminology, avoid contractions, preserve legal and compliance nuance."),
            StyleProfile("Technical Product", "Technical", "Product teams", "Preserve product names, UI labels, placeholders, code, and measurement units."),
            StyleProfile("Customer Friendly", "Conversational", "Customers", "Use natural language while keeping approved terminology and meaning intact."),
        ]
        for profile in defaults:
            conn.execute(
                "INSERT OR IGNORE INTO style_profiles(name, tone, audience, rules) VALUES (?, ?, ?, ?)",
                (profile.name, profile.tone, profile.audience, profile.rules),
            )

    def get_memories(self, source_lang: str, target_lang: str, domain: str = "General") -> list[dict]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT * FROM translation_memory
                WHERE source_lang = ? AND target_lang = ? AND (domain = ? OR domain = 'General')
                ORDER BY updated_at DESC
                """,
                (source_lang, target_lang, domain),
            ).fetchall()
        return [dict(row) for row in rows]

    def upsert_memory(self, source_text: str, target_text: str, source_lang: str, target_lang: str, domain: str = "General") -> None:
        now = utc_now()
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO translation_memory(source_text, target_text, source_lang, target_lang, domain, created_at, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(source_text, source_lang, target_lang, domain)
                DO UPDATE SET target_text = excluded.target_text, updated_at = excluded.updated_at, status = 'approved'
                """,
                (source_text, target_text, source_lang, target_lang, domain, now, now),
            )
            self.add_audit(conn, "memory_upsert", source_text, target_text, f"{source_lang}->{target_lang};{domain}")

    def import_memories(self, rows: Iterable[dict]) -> int:
        count = 0
        for row in rows:
            source = (row.get("source") or row.get("source_text") or "").strip()
            target = (row.get("target") or row.get("target_text") or "").strip()
            if not source or not target:
                continue
            self.upsert_memory(
                source,
                target,
                row.get("source_lang", "English"),
                row.get("target_lang", "Spanish"),
                row.get("domain", "General") or "General",
            )
            count += 1
        return count

    def list_glossary(self, source_lang: str, target_lang: str, domain: str = "General") -> list[GlossaryEntry]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT * FROM glossary
                WHERE source_lang = ? AND target_lang = ? AND (domain = ? OR domain = 'General')
                ORDER BY lower(source_term)
                """,
                (source_lang, target_lang, domain),
            ).fetchall()
        return [
            GlossaryEntry(
                row["source_term"], row["target_term"], row["source_lang"], row["target_lang"],
                row["domain"], bool(row["case_sensitive"]), row["notes"]
            )
            for row in rows
        ]

    def upsert_glossary(self, entry: GlossaryEntry) -> None:
        now = utc_now()
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO glossary(source_term, target_term, source_lang, target_lang, domain, case_sensitive, notes, created_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(source_term, source_lang, target_lang, domain)
                DO UPDATE SET target_term = excluded.target_term, case_sensitive = excluded.case_sensitive, notes = excluded.notes
                """,
                (
                    entry.source_term, entry.target_term, entry.source_lang, entry.target_lang, entry.domain,
                    int(entry.case_sensitive), entry.notes, now,
                ),
            )
            self.add_audit(conn, "glossary_upsert", entry.source_term, entry.target_term, entry.domain)

    def list_profiles(self) -> list[StyleProfile]:
        with self.connect() as conn:
            rows = conn.execute("SELECT name, tone, audience, rules FROM style_profiles ORDER BY name").fetchall()
        return [StyleProfile(row["name"], row["tone"], row["audience"], row["rules"]) for row in rows]

    def add_profile(self, profile: StyleProfile) -> None:
        with self.connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO style_profiles(name, tone, audience, rules) VALUES (?, ?, ?, ?)",
                (profile.name, profile.tone, profile.audience, profile.rules),
            )

    def audit_events(self, limit: int = 50) -> list[dict]:
        with self.connect() as conn:
            rows = conn.execute("SELECT * FROM audit_log ORDER BY id DESC LIMIT ?", (limit,)).fetchall()
        return [dict(row) for row in rows]

    def add_audit(self, conn: sqlite3.Connection, event_type: str, source_text: str, target_text: str, metadata: str = "") -> None:
        conn.execute(
            "INSERT INTO audit_log(event_type, source_text, target_text, metadata, created_at) VALUES (?, ?, ?, ?, ?)",
            (event_type, source_text, target_text, metadata, utc_now()),
        )

