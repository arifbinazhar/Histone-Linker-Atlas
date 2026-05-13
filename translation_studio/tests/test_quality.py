from translation_studio.models import Segment
from translation_studio.qa import analyze_source_quality
from translation_studio.tm import find_translation_candidates


def test_quality_flags_spacing_and_punctuation():
    segment = Segment("Reset  password now!!", "paragraph", 1)
    issues = analyze_source_quality([segment])
    categories = {issue.category for issue in issues}
    assert "Spacing" in categories
    assert "Punctuation" in categories


def test_translation_memory_exact_match_scores_highest():
    memories = [
        {"id": 1, "source_text": "Reset password", "target_text": "Restablecer contraseña"},
        {"id": 2, "source_text": "Change account email", "target_text": "Cambiar correo de la cuenta"},
    ]
    candidates = find_translation_candidates("Reset password", memories)
    assert candidates[0].score == 1.0
    assert candidates[0].translation == "Restablecer contraseña"

