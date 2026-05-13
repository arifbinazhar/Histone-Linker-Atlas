from __future__ import annotations

import re
from collections import defaultdict
from difflib import SequenceMatcher

from .models import GlossaryEntry, QualityIssue, Segment, Severity


DOUBLE_SPACE_RE = re.compile(r" {2,}")
MISSING_SPACE_RE = re.compile(r"([a-z])([A-Z])")
REPEATED_PUNCT_RE = re.compile(r"([,.;:!?])\1+")
SPACE_BEFORE_PUNCT_RE = re.compile(r"\s+([,.;:!?])")
UNBALANCED_RE = re.compile(r"[\(\)\[\]{}]")
TOKEN_RE = re.compile(r"\b[A-Za-z][A-Za-z0-9.\-]{2,}\b")


COMMON_WORDS = {
    "the", "and", "for", "with", "from", "this", "that", "are", "you", "your", "will", "shall",
    "document", "translation", "project", "approved", "source", "target", "quality", "review",
    "password", "account", "user", "system", "data", "file", "language", "table", "policy",
}


def analyze_source_quality(segments: list[Segment], glossary: list[GlossaryEntry] | None = None) -> list[QualityIssue]:
    issues: list[QualityIssue] = []
    glossary = glossary or []
    issues.extend(check_segment_level_quality(segments))
    issues.extend(check_terminology_variants(segments, glossary))
    return issues


def check_segment_level_quality(segments: list[Segment]) -> list[QualityIssue]:
    issues: list[QualityIssue] = []
    for segment in segments:
        text = segment.text
        if DOUBLE_SPACE_RE.search(text):
            issues.append(issue(segment, Severity.WARNING, "Spacing", "Multiple spaces found.", "Collapse repeated spaces."))
        if MISSING_SPACE_RE.search(text):
            issues.append(issue(segment, Severity.WARNING, "Spacing", "Possible missing space between words.", "Review camel-cased source text."))
        if REPEATED_PUNCT_RE.search(text):
            issues.append(issue(segment, Severity.ERROR, "Punctuation", "Repeated punctuation found.", "Keep one punctuation mark unless stylistically required."))
        if SPACE_BEFORE_PUNCT_RE.search(text):
            issues.append(issue(segment, Severity.WARNING, "Punctuation", "Space before punctuation.", "Remove the extra space before punctuation."))
        if has_unbalanced_pairs(text):
            issues.append(issue(segment, Severity.ERROR, "Formatting", "Unbalanced brackets or parentheses.", "Fix the missing opening or closing marker."))
        if text and segment.block_type != "heading" and text[-1].isalnum() and len(text.split()) > 5:
            issues.append(issue(segment, Severity.INFO, "Punctuation", "Sentence may be missing terminal punctuation.", "Add a final period if this is a complete sentence."))
        for typo in suspicious_tokens(text):
            issues.append(issue(segment, Severity.WARNING, "Spelling", f"Suspicious token: '{typo}'.", "Check spelling or add the term to the glossary."))
    return issues


def check_terminology_variants(segments: list[Segment], glossary: list[GlossaryEntry]) -> list[QualityIssue]:
    issues: list[QualityIssue] = []
    surface_forms: dict[str, set[str]] = defaultdict(set)
    segment_by_form: dict[str, str] = {}

    glossary_terms = {entry.source_term.lower() for entry in glossary}
    for segment in segments:
        for token in TOKEN_RE.findall(segment.text):
            if token.lower() in COMMON_WORDS or token.lower() in glossary_terms:
                continue
            normalized = normalize_term(token)
            if len(normalized) >= 2:
                surface_forms[normalized].add(token)
                segment_by_form[token] = segment.id

    for _, variants in surface_forms.items():
        if len(variants) < 2:
            continue
        sorted_variants = sorted(variants)
        if likely_term_variants(sorted_variants):
            segment_id = segment_by_form[sorted_variants[0]]
            issues.append(
                QualityIssue(
                    segment_id=segment_id,
                    severity=Severity.WARNING,
                    category="Terminology",
                    message=f"Possible terminology variants: {', '.join(sorted_variants)}.",
                    suggestion="Choose one approved spelling/casing and add it to the glossary.",
                )
            )
    return issues


def suspicious_tokens(text: str) -> list[str]:
    tokens = []
    for token in TOKEN_RE.findall(text):
        lower = token.lower().strip(".")
        if lower in COMMON_WORDS or len(lower) < 8:
            continue
        if re.search(r"(.)\1\1", lower) or SequenceMatcher(None, lower, lower[::-1]).ratio() > 0.86:
            tokens.append(token)
    return tokens[:3]


def has_unbalanced_pairs(text: str) -> bool:
    pairs = {"(": ")", "[": "]", "{": "}"}
    stack: list[str] = []
    for char in text:
        if char in pairs:
            stack.append(pairs[char])
        elif char in pairs.values():
            if not stack or stack.pop() != char:
                return True
    return bool(stack)


def normalize_term(token: str) -> str:
    return re.sub(r"[^a-z0-9]", "", token.lower())


def likely_term_variants(variants: list[str]) -> bool:
    if len({variant.lower() for variant in variants}) > 1:
        return True
    return len(set(variants)) > 1


def issue(segment: Segment, severity: Severity, category: str, message: str, suggestion: str) -> QualityIssue:
    return QualityIssue(segment.id, severity, category, message, suggestion)

