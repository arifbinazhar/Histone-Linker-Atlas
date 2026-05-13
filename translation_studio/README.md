# Translation Studio

An intelligent, local-first enterprise translation workflow app. It accepts PDF and DOCX uploads, checks source quality, reuses approved translations through a translation memory, enforces glossaries, sends only new segments to an LLM, and captures approval history for continuous learning.

## What This Implements

- Upload and parse PDF, DOCX, and TXT files into translatable segments.
- Preserve document structure labels such as headings, paragraphs, table cells, and list items.
- Validate source quality with spelling, punctuation, spacing, terminology, and formatting checks.
- Store approved translations in SQLite as a translation memory.
- Retrieve exact and fuzzy TM matches with local lexical similarity; if `scikit-learn` is installed, the matcher automatically adds TF-IDF vector similarity.
- Translate new segments with an LLM when configured.
- Enforce glossary terms per source language, target language, and domain.
- Configure style profiles for formal, official, conversational, technical, or social tone.
- Review segments side by side, edit translations, approve/reject, and keep an audit trail.
- Import bilingual CSV/TMX-style data to seed the translation memory.

## Repository Layout

```text
translation_studio/
  app.py                         # Streamlit UI
  requirements.txt               # App dependencies
  .env.example                   # Optional LLM configuration
  README.md
  data/
    .gitkeep
  src/
    translation_studio/
      __init__.py
      models.py                  # Dataclasses and workflow entities
      parsing.py                 # PDF/DOCX/TXT extraction
      qa.py                      # Source quality checks
      storage.py                 # SQLite persistence and audit trail
      tm.py                      # Translation memory search and scoring
      translator.py              # LLM and local fallback translation
      workflow.py                # Orchestration helpers
  tests/
    test_quality.py
```

## Quick Start

```bash
cd translation_studio
python app_server.py
```

Then open:

```text
http://localhost:8501
```

This default server uses only the Python standard library plus the core code in this repository. It can parse TXT files immediately, DOCX files through a built-in XML parser fallback, and PDF files best when `pypdf` is installed.

For the richer Streamlit UI, install optional dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Optional LLM Setup

The app works without an API key, but production-quality translation requires an LLM. It supports OpenAI and Gemini.

1. Copy `.env.example` to `.env`.
2. Set either `GEMINI_API_KEY`, `OPENAI_API_KEY`, or both.
3. Choose `TRANSLATION_PROVIDER=auto`, `gemini`, `openai`, or `local`.

```env
TRANSLATION_PROVIDER=auto

GEMINI_API_KEY=your-gemini-key
GEMINI_MODEL=gemini-2.5-flash

OPENAI_API_KEY=your-openai-key
OPENAI_MODEL=gpt-4.1-mini
```

In `auto` mode, the app tries Gemini first when `GEMINI_API_KEY` exists, then OpenAI when `OPENAI_API_KEY` exists, then the deterministic local fallback translator. The fallback is useful for demos and testing the workflow, but it is not intended to replace a real translation model.

## Bilingual Import Format

Use the "Memory & Glossary" tab to import approved translations from CSV. Required columns:

```csv
source,target,source_lang,target_lang,domain
"Reset password","Restablecer contraseña","English","Spanish","General"
```

## Deployment Notes

This app is local-first and stores data in `translation_studio/data/translation_studio.sqlite3`.

For team deployment:

- Put the SQLite database on persistent storage, or replace `storage.py` with PostgreSQL.
- Put the app behind SSO/reverse proxy authentication.
- Use a managed LLM endpoint and restrict logs for confidential documents.
- Add scheduled exports of approved TM and glossary entries.
