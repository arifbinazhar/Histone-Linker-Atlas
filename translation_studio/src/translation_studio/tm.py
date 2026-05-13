from __future__ import annotations

from difflib import SequenceMatcher

from .models import TranslationCandidate


def find_translation_candidates(
    source_text: str,
    memories: list[dict],
    min_score: float = 0.55,
    limit: int = 5,
) -> list[TranslationCandidate]:
    if not memories:
        return []

    vector_scores = vector_similarity(source_text, [row["source_text"] for row in memories])
    candidates: list[TranslationCandidate] = []
    for index, memory in enumerate(memories):
        lexical = SequenceMatcher(None, normalize(source_text), normalize(memory["source_text"])).ratio()
        score = max(lexical, vector_scores[index])
        if normalize(source_text) == normalize(memory["source_text"]):
            score = 1.0
        if score >= min_score:
            origin = "tm_exact" if score >= 0.99 else "tm_fuzzy"
            candidates.append(
                TranslationCandidate(
                    translation=memory["target_text"],
                    source=memory["source_text"],
                    score=round(score, 3),
                    origin=origin,
                    memory_id=memory["id"],
                )
            )
    return sorted(candidates, key=lambda item: item.score, reverse=True)[:limit]


def vector_similarity(query: str, corpus: list[str]) -> list[float]:
    try:
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.metrics.pairwise import cosine_similarity
    except ImportError:
        return [SequenceMatcher(None, normalize(query), normalize(text)).ratio() for text in corpus]

    vectorizer = TfidfVectorizer(analyzer="char_wb", ngram_range=(3, 5), lowercase=True)
    matrix = vectorizer.fit_transform([query, *corpus])
    scores = cosine_similarity(matrix[0:1], matrix[1:]).flatten()
    return [float(score) for score in scores]


def normalize(text: str) -> str:
    return " ".join(text.lower().split())

