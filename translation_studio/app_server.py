from __future__ import annotations

import csv
import html
import io
import os
import sys
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def load_env_file(path: Path) -> None:
    if not path.exists():
        return
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if key and key not in os.environ:
            os.environ[key] = value


load_env_file(ROOT / ".env")

from translation_studio.models import GlossaryEntry, SegmentStatus
from translation_studio.parsing import parse_uploaded_file
from translation_studio.qa import analyze_source_quality
from translation_studio.storage import StudioStore
from translation_studio.workflow import prepare_translation_units


DB_PATH = ROOT / "data" / "translation_studio.sqlite3"
LANGUAGES = ["Spanish", "Japanese", "French", "German", "Hindi", "Arabic", "Portuguese", "Italian"]
STATE = {"segments": [], "issues": [], "units": []}


class TranslationStudioHandler(BaseHTTPRequestHandler):
    store = StudioStore(DB_PATH)

    def do_GET(self) -> None:
        route = urllib.parse.urlparse(self.path).path
        if route == "/memory.csv":
            self.export_memory()
            return
        self.respond(render_page(self.store))

    def do_POST(self) -> None:
        route = urllib.parse.urlparse(self.path).path
        if route == "/upload":
            self.handle_upload()
        elif route == "/approve":
            self.handle_approve()
        elif route == "/glossary":
            self.handle_glossary()
        elif route == "/import-memory":
            self.handle_import_memory()
        else:
            self.send_error(404)

    def handle_upload(self) -> None:
        fields, files = parse_multipart(self)
        source_lang = fields.get("source_lang", "English")
        target_lang = fields.get("target_lang", "Spanish")
        domain = fields.get("domain", "General") or "General"
        provider = fields.get("provider", os.getenv("TRANSLATION_PROVIDER", "auto")).lower()
        profile_name = fields.get("profile", "Enterprise Formal")
        profile = next((item for item in self.store.list_profiles() if item.name == profile_name), self.store.list_profiles()[0])
        uploaded = files.get("document")
        if not uploaded:
            self.respond(render_page(self.store, "No file was uploaded."), status=400)
            return

        try:
            segments = parse_uploaded_file(io.BytesIO(uploaded["content"]), uploaded["filename"])
            glossary = self.store.list_glossary(source_lang, target_lang, domain)
            issues = analyze_source_quality(segments, glossary)
            units = prepare_translation_units(
                segments, self.store, source_lang, target_lang, domain, profile, auto_translate=True, provider=provider
            )
            STATE.update({
                "segments": segments,
                "issues": issues,
                "units": units,
                "source_lang": source_lang,
                "target_lang": target_lang,
                "domain": domain,
                "provider": provider,
            })
            self.respond(render_page(self.store, f"Parsed {len(segments)} segments and prepared {len(units)} translations."))
        except Exception as exc:
            self.respond(render_page(self.store, f"Upload failed: {exc}"), status=400)

    def handle_approve(self) -> None:
        length = int(self.headers.get("Content-Length", "0"))
        data = urllib.parse.parse_qs(self.rfile.read(length).decode("utf-8", errors="ignore"))
        index = int(data.get("index", ["0"])[0])
        translation = data.get("translation", [""])[0]
        units = STATE.get("units", [])
        if 0 <= index < len(units):
            unit = units[index]
            unit.translation = translation
            unit.status = SegmentStatus.APPROVED
            self.store.upsert_memory(
                unit.segment.text,
                translation,
                STATE.get("source_lang", unit.source_lang),
                STATE.get("target_lang", unit.target_lang),
                STATE.get("domain", "General"),
            )
            self.respond(render_page(self.store, "Approved segment and updated translation memory."))
        else:
            self.respond(render_page(self.store, "Invalid segment index."), status=400)

    def handle_glossary(self) -> None:
        length = int(self.headers.get("Content-Length", "0"))
        data = urllib.parse.parse_qs(self.rfile.read(length).decode("utf-8", errors="ignore"))
        entry = GlossaryEntry(
            data.get("source_term", [""])[0].strip(),
            data.get("target_term", [""])[0].strip(),
            data.get("source_lang", ["English"])[0],
            data.get("target_lang", ["Spanish"])[0],
            data.get("domain", ["General"])[0] or "General",
            bool(data.get("case_sensitive")),
            data.get("notes", [""])[0],
        )
        if entry.source_term and entry.target_term:
            self.store.upsert_glossary(entry)
            self.respond(render_page(self.store, "Glossary entry saved."))
        else:
            self.respond(render_page(self.store, "Source and target terms are required."), status=400)

    def handle_import_memory(self) -> None:
        _, files = parse_multipart(self)
        uploaded = files.get("memory")
        if not uploaded:
            self.respond(render_page(self.store, "No CSV file uploaded."), status=400)
            return
        text = uploaded["content"].decode("utf-8-sig", errors="ignore")
        count = self.store.import_memories(csv.DictReader(io.StringIO(text)))
        self.respond(render_page(self.store, f"Imported {count} translation memory rows."))

    def export_memory(self) -> None:
        memories = self.store.get_memories("English", "Spanish", "General")
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=["source_text", "target_text", "source_lang", "target_lang", "domain"])
        writer.writeheader()
        for row in memories:
            writer.writerow({key: row[key] for key in writer.fieldnames})
        payload = output.getvalue().encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/csv; charset=utf-8")
        self.send_header("Content-Disposition", "attachment; filename=translation_memory.csv")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def respond(self, body: str, status: int = 200) -> None:
        payload = body.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)


def render_page(store: StudioStore, notice: str = "") -> str:
    profiles = store.list_profiles()
    units = STATE.get("units", [])
    issues = STATE.get("issues", [])
    audit = store.audit_events(20)
    target_options = "".join(f"<option>{lang}</option>" for lang in LANGUAGES)
    provider_options = "".join(
        f"<option value='{value}'>{label}</option>"
        for value, label in [
            ("auto", "Auto: Gemini, OpenAI, Local"),
            ("gemini", "Gemini only"),
            ("openai", "OpenAI only"),
            ("local", "Local draft only"),
        ]
    )
    profile_options = "".join(f"<option>{html.escape(profile.name)}</option>" for profile in profiles)
    issue_html = render_issues(issues)
    unit_html = render_units(units)
    audit_html = "".join(
        f"<tr><td>{event['created_at']}</td><td>{html.escape(event['event_type'])}</td><td>{html.escape(event['source_text'][:80])}</td></tr>"
        for event in audit
    )
    notice_html = f"<div class='notice'>{html.escape(notice)}</div>" if notice else ""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Translation Studio</title>
  <style>
    :root {{ color-scheme: light; --ink:#17202a; --muted:#667085; --line:#d0d7de; --brand:#0f766e; --soft:#f6f8fa; }}
    body {{ margin:0; font-family: Inter, Segoe UI, Arial, sans-serif; color:var(--ink); background:#fff; }}
    header {{ padding:24px 32px 16px; border-bottom:1px solid var(--line); }}
    h1 {{ margin:0; font-size:30px; letter-spacing:0; }}
    h2 {{ margin:26px 0 12px; font-size:20px; }}
    main {{ max-width:1220px; margin:0 auto; padding:22px 28px 48px; }}
    .grid {{ display:grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap:12px; }}
    .panel {{ border:1px solid var(--line); border-radius:8px; padding:16px; background:#fff; }}
    .notice {{ border-left:4px solid var(--brand); background:#ecfdf5; padding:12px 14px; margin-bottom:18px; }}
    label {{ display:block; font-weight:650; font-size:13px; margin-bottom:6px; }}
    input, select, textarea {{ width:100%; box-sizing:border-box; border:1px solid var(--line); border-radius:6px; padding:9px 10px; font:inherit; }}
    button {{ background:var(--brand); color:white; border:0; border-radius:6px; padding:10px 14px; font-weight:700; cursor:pointer; }}
    table {{ width:100%; border-collapse:collapse; font-size:14px; }}
    th, td {{ border-bottom:1px solid var(--line); padding:8px; text-align:left; vertical-align:top; }}
    .muted {{ color:var(--muted); }}
    .segment {{ display:grid; grid-template-columns:1fr 1fr; gap:14px; margin-bottom:14px; }}
    .badge {{ display:inline-block; padding:3px 7px; border-radius:999px; background:var(--soft); border:1px solid var(--line); font-size:12px; }}
    .error {{ border-left:4px solid #b42318; }} .warning {{ border-left:4px solid #b54708; }} .info {{ border-left:4px solid #175cd3; }}
    @media (max-width: 850px) {{ .grid, .segment {{ grid-template-columns:1fr; }} main {{ padding:18px; }} }}
  </style>
</head>
<body>
<header>
  <h1>Translation Studio</h1>
  <p class="muted">Upload, validate, reuse memory, enforce glossaries, translate, approve, and learn continuously.</p>
</header>
<main>
  {notice_html}
  <section class="panel">
    <h2>Upload & Configure</h2>
    <form action="/upload" method="post" enctype="multipart/form-data">
      <div class="grid">
        <div><label>Document</label><input type="file" name="document" accept=".pdf,.docx,.txt,.md" required></div>
        <div><label>Source language</label><select name="source_lang"><option>English</option><option>Spanish</option><option>French</option><option>German</option><option>Japanese</option></select></div>
        <div><label>Target language</label><select name="target_lang">{target_options}</select></div>
        <div><label>Domain</label><input name="domain" value="General"></div>
      </div>
      <div class="grid">
        <div><label>Provider</label><select name="provider">{provider_options}</select></div>
        <div><label>Style profile</label><select name="profile">{profile_options}</select></div>
      </div>
      <button type="submit">Parse, QA, Match TM, and Translate</button>
    </form>
  </section>
  <section>
    <h2>Source Quality</h2>
    {issue_html}
  </section>
  <section>
    <h2>Side-by-Side Review</h2>
    {unit_html}
  </section>
  <section class="panel">
    <h2>Glossary</h2>
    <form action="/glossary" method="post">
      <div class="grid">
        <div><label>Source term</label><input name="source_term"></div>
        <div><label>Target term</label><input name="target_term"></div>
        <div><label>Source language</label><input name="source_lang" value="English"></div>
        <div><label>Target language</label><select name="target_lang">{target_options}</select></div>
      </div>
      <p><label>Notes</label><input name="notes"></p>
      <button type="submit">Save Glossary Term</button>
    </form>
  </section>
  <section class="panel">
    <h2>Import Existing Bilingual Memory</h2>
    <form action="/import-memory" method="post" enctype="multipart/form-data">
      <input type="file" name="memory" accept=".csv" required>
      <p class="muted">CSV columns: source,target,source_lang,target_lang,domain</p>
      <button type="submit">Import Memory</button>
    </form>
  </section>
  <section>
    <h2>Audit Trail</h2>
    <table><thead><tr><th>Time</th><th>Event</th><th>Source</th></tr></thead><tbody>{audit_html}</tbody></table>
  </section>
</main>
</body>
</html>"""


def render_issues(issues) -> str:
    if not issues:
        return "<p class='muted'>No quality issues loaded yet.</p>"
    rows = []
    for issue in issues[:80]:
        rows.append(
            f"<div class='panel {issue.severity.value}'><b>{html.escape(issue.severity.value.upper())} - {html.escape(issue.category)}</b>"
            f"<p>{html.escape(issue.message)}</p><p class='muted'>{html.escape(issue.suggestion)}</p></div>"
        )
    return "".join(rows)


def render_units(units) -> str:
    if not units:
        return "<p class='muted'>Upload a document to prepare review segments.</p>"
    blocks = []
    for index, unit in enumerate(units):
        suggestions = "".join(
            f"<li>{candidate.score:.0%}: {html.escape(candidate.translation)}</li>"
            for candidate in unit.candidates
        ) or "<li>No memory suggestion.</li>"
        glossary = ", ".join(entry.source_term for entry in unit.glossary_hits) or "None"
        blocks.append(f"""
        <form class="panel segment" action="/approve" method="post">
          <input type="hidden" name="index" value="{index}">
          <div>
            <span class="badge">{html.escape(unit.status.value)}</span>
            <p><b>Source {unit.segment.index}</b> · {html.escape(unit.segment.block_type)}</p>
            <p>{html.escape(unit.segment.text)}</p>
            <p class="muted">Glossary hits: {html.escape(glossary)}</p>
            <ul>{suggestions}</ul>
          </div>
          <div>
            <label>Reviewed translation</label>
            <textarea name="translation" rows="8">{html.escape(unit.translation)}</textarea>
            <p><button type="submit">Approve & Save to TM</button></p>
          </div>
        </form>""")
    return "".join(blocks)


def parse_multipart(handler: BaseHTTPRequestHandler) -> tuple[dict[str, str], dict[str, dict]]:
    content_type = handler.headers.get("Content-Type", "")
    boundary_marker = "boundary="
    if boundary_marker not in content_type:
        return {}, {}
    boundary = ("--" + content_type.split(boundary_marker, 1)[1]).encode()
    length = int(handler.headers.get("Content-Length", "0"))
    body = handler.rfile.read(length)
    fields: dict[str, str] = {}
    files: dict[str, dict] = {}

    for part in body.split(boundary):
        part = part.strip(b"\r\n")
        if not part or part == b"--" or b"\r\n\r\n" not in part:
            continue
        raw_headers, content = part.split(b"\r\n\r\n", 1)
        content = content.rstrip(b"\r\n")
        headers = raw_headers.decode("utf-8", errors="ignore")
        disposition = next((line for line in headers.split("\r\n") if line.lower().startswith("content-disposition")), "")
        attrs = parse_disposition(disposition)
        name = attrs.get("name")
        if not name:
            continue
        if "filename" in attrs:
            files[name] = {"filename": attrs.get("filename") or "upload.bin", "content": content}
        else:
            fields[name] = content.decode("utf-8", errors="ignore")
    return fields, files


def parse_disposition(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, raw = item.strip().split("=", 1)
            attrs[key] = raw.strip('"')
    return attrs


def run(host: str = "127.0.0.1", port: int = 8501) -> None:
    server = ThreadingHTTPServer((host, port), TranslationStudioHandler)
    print(f"Translation Studio running at http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    run()
