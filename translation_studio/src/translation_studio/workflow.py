from __future__ import annotations

from collections import defaultdict

from .models import SegmentStatus, StyleProfile, TranslationUnit
from .qa import analyze_source_quality
from .storage import StudioStore
from .tm import find_translation_candidates
from .translator import TranslationEngine
from .translator import TranslatorConfig


def prepare_translation_units(
    segments,
    store: StudioStore,
    source_lang: str,
    target_lang: str,
    domain: str,
    profile: StyleProfile,
    auto_translate: bool = True,
    provider: str = "auto",
) -> list[TranslationUnit]:
    glossary = store.list_glossary(source_lang, target_lang, domain)
    memories = store.get_memories(source_lang, target_lang, domain)
    issues = analyze_source_quality(segments, glossary)
    issues_by_segment = defaultdict(list)
    for issue in issues:
        issues_by_segment[issue.segment_id].append(issue)

    engine = TranslationEngine(TranslatorConfig(provider=provider))
    units: list[TranslationUnit] = []
    for segment in segments:
        candidates = find_translation_candidates(segment.text, memories)
        glossary_hits = [
            entry for entry in glossary
            if (entry.source_term in segment.text if entry.case_sensitive else entry.source_term.lower() in segment.text.lower())
        ]
        unit = TranslationUnit(
            segment=segment,
            source_lang=source_lang,
            target_lang=target_lang,
            candidates=candidates,
            quality_issues=issues_by_segment.get(segment.id, []),
            glossary_hits=glossary_hits,
        )
        if candidates and candidates[0].score >= 0.99:
            unit.translation = candidates[0].translation
            unit.status = SegmentStatus.TM_EXACT
        elif candidates and candidates[0].score >= 0.72:
            unit.translation = candidates[0].translation
            unit.status = SegmentStatus.TM_FUZZY
        elif auto_translate:
            unit.translation = engine.translate(segment.text, source_lang, target_lang, glossary_hits, profile)
            unit.status = SegmentStatus.MACHINE_TRANSLATED
        units.append(unit)
    return units
