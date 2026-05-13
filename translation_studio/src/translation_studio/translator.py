from __future__ import annotations

import os
import re
from dataclasses import dataclass

from .models import GlossaryEntry, StyleProfile


LANGUAGE_HINTS = {
    "Spanish": "es",
    "French": "fr",
    "German": "de",
    "Japanese": "ja",
    "Hindi": "hi",
    "Arabic": "ar",
    "Portuguese": "pt",
    "Italian": "it",
}


@dataclass(slots=True)
class TranslatorConfig:
    provider: str = "auto"
    model: str = os.getenv("OPENAI_MODEL", "gpt-4.1-mini")


class TranslationEngine:
    def __init__(self, config: TranslatorConfig | None = None):
        self.config = config or TranslatorConfig()

    def translate(
        self,
        text: str,
        source_lang: str,
        target_lang: str,
        glossary: list[GlossaryEntry],
        profile: StyleProfile,
    ) -> str:
        if os.getenv("OPENAI_API_KEY"):
            try:
                return self._translate_with_openai(text, source_lang, target_lang, glossary, profile)
            except Exception as exc:  # Keep local demos resilient.
                return enforce_glossary(local_fallback_translate(text, target_lang), glossary) + f" [LLM fallback: {exc.__class__.__name__}]"
        return enforce_glossary(local_fallback_translate(text, target_lang), glossary)

    def _translate_with_openai(
        self,
        text: str,
        source_lang: str,
        target_lang: str,
        glossary: list[GlossaryEntry],
        profile: StyleProfile,
    ) -> str:
        from openai import OpenAI

        glossary_text = "\n".join(f"- {entry.source_term} => {entry.target_term}" for entry in glossary) or "None"
        client = OpenAI()
        response = client.chat.completions.create(
            model=self.config.model,
            temperature=0.1,
            messages=[
                {
                    "role": "system",
                    "content": (
                        "You are an enterprise translation engine. Translate faithfully, preserve placeholders, "
                        "numbers, markdown, tags, product names, and formatting. Use the glossary exactly."
                    ),
                },
                {
                    "role": "user",
                    "content": (
                        f"Source language: {source_lang}\n"
                        f"Target language: {target_lang}\n"
                        f"Tone: {profile.tone}\n"
                        f"Audience: {profile.audience}\n"
                        f"Style rules: {profile.rules}\n"
                        f"Glossary:\n{glossary_text}\n\n"
                        f"Translate this segment only:\n{text}"
                    ),
                },
            ],
        )
        translated = response.choices[0].message.content or ""
        return enforce_glossary(translated.strip(), glossary)


def enforce_glossary(text: str, glossary: list[GlossaryEntry]) -> str:
    output = text
    for entry in glossary:
        flags = 0 if entry.case_sensitive else re.IGNORECASE
        if re.search(rf"\b{re.escape(entry.source_term)}\b", output, flags=flags):
            output = re.sub(rf"\b{re.escape(entry.source_term)}\b", entry.target_term, output, flags=flags)
    return output


def local_fallback_translate(text: str, target_lang: str) -> str:
    lexicons = {
        "Spanish": {
            "document": "documento", "translation": "traducción", "quality": "calidad", "password": "contraseña",
            "account": "cuenta", "user": "usuario", "approved": "aprobado", "source": "origen", "target": "destino",
        },
        "French": {
            "document": "document", "translation": "traduction", "quality": "qualité", "password": "mot de passe",
            "account": "compte", "user": "utilisateur", "approved": "approuvé", "source": "source", "target": "cible",
        },
        "German": {
            "document": "Dokument", "translation": "Übersetzung", "quality": "Qualität", "password": "Passwort",
            "account": "Konto", "user": "Benutzer", "approved": "genehmigt", "source": "Quelle", "target": "Ziel",
        },
        "Japanese": {
            "document": "文書", "translation": "翻訳", "quality": "品質", "password": "パスワード",
            "account": "アカウント", "user": "ユーザー", "approved": "承認済み", "source": "原文", "target": "訳文",
        },
    }
    words = lexicons.get(target_lang, {})

    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        replacement = words.get(token.lower())
        if not replacement:
            return token
        return replacement.capitalize() if token[:1].isupper() else replacement

    translated = re.sub(r"\b[A-Za-z]+\b", replace, text)
    return f"{translated} [{target_lang} draft]"

