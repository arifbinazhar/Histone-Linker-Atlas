from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any
from uuid import uuid4


class Severity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class SegmentStatus(str, Enum):
    NEW = "new"
    TM_EXACT = "tm_exact"
    TM_FUZZY = "tm_fuzzy"
    MACHINE_TRANSLATED = "machine_translated"
    APPROVED = "approved"
    REJECTED = "rejected"


@dataclass(slots=True)
class Segment:
    text: str
    block_type: str
    index: int
    section: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    id: str = field(default_factory=lambda: str(uuid4()))


@dataclass(slots=True)
class QualityIssue:
    segment_id: str
    severity: Severity
    category: str
    message: str
    suggestion: str = ""
    span: tuple[int, int] | None = None


@dataclass(slots=True)
class GlossaryEntry:
    source_term: str
    target_term: str
    source_lang: str
    target_lang: str
    domain: str = "General"
    case_sensitive: bool = False
    notes: str = ""


@dataclass(slots=True)
class StyleProfile:
    name: str
    tone: str = "Formal"
    audience: str = "Enterprise"
    rules: str = "Prefer concise, consistent, terminology-controlled translations."


@dataclass(slots=True)
class TranslationCandidate:
    translation: str
    source: str
    score: float
    origin: str
    memory_id: int | None = None


@dataclass(slots=True)
class TranslationUnit:
    segment: Segment
    target_lang: str
    source_lang: str = "English"
    translation: str = ""
    status: SegmentStatus = SegmentStatus.NEW
    candidates: list[TranslationCandidate] = field(default_factory=list)
    quality_issues: list[QualityIssue] = field(default_factory=list)
    glossary_hits: list[GlossaryEntry] = field(default_factory=list)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")

