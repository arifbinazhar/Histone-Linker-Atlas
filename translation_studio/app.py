from __future__ import annotations

import csv
import io
import sys
from pathlib import Path
from dataclasses import asdict

import streamlit as st

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

try:
    from dotenv import load_dotenv
    load_dotenv(ROOT / ".env")
except Exception:
    pass

from translation_studio.models import GlossaryEntry, SegmentStatus, StyleProfile
from translation_studio.parsing import parse_uploaded_file
from translation_studio.qa import analyze_source_quality
from translation_studio.storage import StudioStore
from translation_studio.workflow import prepare_translation_units


DB_PATH = ROOT / "data" / "translation_studio.sqlite3"
LANGUAGES = ["Spanish", "Japanese", "French", "German", "Hindi", "Arabic", "Portuguese", "Italian"]


st.set_page_config(page_title="Translation Studio", layout="wide", page_icon="TS")
st.markdown(
    """
    <style>
    .block-container {padding-top: 1.25rem; max-width: 1400px;}
    [data-testid="stMetricValue"] {font-size: 1.65rem;}
    .severity-error {border-left: 4px solid #c62828; padding-left: .75rem;}
    .severity-warning {border-left: 4px solid #ef6c00; padding-left: .75rem;}
    .severity-info {border-left: 4px solid #1565c0; padding-left: .75rem;}
    </style>
    """,
    unsafe_allow_html=True,
)

store = StudioStore(DB_PATH)


def get_state(name, default):
    if name not in st.session_state:
        st.session_state[name] = default
    return st.session_state[name]


def main() -> None:
    st.title("Translation Studio")
    st.caption("Document parsing, quality validation, translation memory, glossary enforcement, LLM translation, and approval in one local workflow.")

    with st.sidebar:
        st.header("Project")
        source_lang = st.selectbox("Source language", ["English", "Spanish", "French", "German", "Japanese"], index=0)
        target_lang = st.selectbox("Target language", LANGUAGES, index=0)
        domain = st.text_input("Domain", "General")
        profiles = store.list_profiles()
        profile_names = [profile.name for profile in profiles]
        selected_profile_name = st.selectbox("Style profile", profile_names)
        profile = next(profile for profile in profiles if profile.name == selected_profile_name)
        provider = st.selectbox("LLM provider", ["auto", "gemini", "openai", "local"], index=0)
        auto_translate = st.toggle("Translate new segments automatically", value=True)

    tab_upload, tab_review, tab_memory, tab_audit = st.tabs(["Upload & QA", "Translate & Approve", "Memory & Glossary", "Audit Trail"])

    with tab_upload:
        uploaded = st.file_uploader("Upload a PDF, DOCX, or TXT document", type=["pdf", "docx", "txt", "md"])
        if uploaded and st.button("Parse and analyze", type="primary"):
            with st.spinner("Extracting structure and validating source quality..."):
                segments = parse_uploaded_file(uploaded, uploaded.name)
                glossary = store.list_glossary(source_lang, target_lang, domain)
                issues = analyze_source_quality(segments, glossary)
                st.session_state["segments"] = segments
                st.session_state["issues"] = issues
                st.session_state["units"] = []
            st.success(f"Parsed {len(segments)} translatable segments and found {len(issues)} source quality observations.")

        segments = get_state("segments", [])
        issues = get_state("issues", [])
        if segments:
            render_quality_dashboard(segments, issues)

    with tab_review:
        segments = get_state("segments", [])
        if not segments:
            st.info("Upload and parse a document first.")
        else:
            if st.button("Run TM matching and translation", type="primary"):
                with st.spinner("Searching translation memory and translating new segments..."):
                    st.session_state["units"] = prepare_translation_units(
                        segments,
                        store,
                        source_lang,
                        target_lang,
                        domain,
                        profile,
                        auto_translate,
                        provider,
                    )
            units = get_state("units", [])
            if units:
                render_review(units, store, source_lang, target_lang, domain)

    with tab_memory:
        render_memory_glossary(store, source_lang, target_lang, domain)
        render_profile_editor(store)

    with tab_audit:
        st.subheader("Audit Trail")
        events = store.audit_events()
        if events:
            st.dataframe(events, use_container_width=True, hide_index=True)
        else:
            st.info("No audit events yet. Approvals and glossary updates will appear here.")


def render_quality_dashboard(segments, issues) -> None:
    st.subheader("Source Quality")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Segments", len(segments))
    col2.metric("Errors", sum(1 for item in issues if item.severity.value == "error"))
    col3.metric("Warnings", sum(1 for item in issues if item.severity.value == "warning"))
    col4.metric("Info", sum(1 for item in issues if item.severity.value == "info"))

    issue_by_segment = {}
    for item in issues:
        issue_by_segment.setdefault(item.segment_id, []).append(item)

    for segment in segments[:80]:
        with st.expander(f"{segment.index}. {segment.block_type} - {segment.text[:90]}"):
            st.write(segment.text)
            segment_issues = issue_by_segment.get(segment.id, [])
            if not segment_issues:
                st.success("No quality issues detected.")
            for item in segment_issues:
                st.markdown(
                    f"<div class='severity-{item.severity.value}'><b>{item.severity.value.upper()} - {item.category}</b><br>{item.message}<br><small>{item.suggestion}</small></div>",
                    unsafe_allow_html=True,
                )


def render_review(units, store: StudioStore, source_lang: str, target_lang: str, domain: str) -> None:
    st.subheader("Side-by-Side Review")
    approved = 0
    rejected = 0

    for index, unit in enumerate(units):
        with st.container(border=True):
            header_cols = st.columns([5, 2, 2])
            header_cols[0].markdown(f"**Segment {unit.segment.index}** · `{unit.segment.block_type}`")
            header_cols[1].markdown(f"`{unit.status.value.replace('_', ' ').title()}`")
            header_cols[2].markdown(f"**Glossary hits:** {len(unit.glossary_hits)}")

            left, right = st.columns(2)
            with left:
                st.caption("Source")
                st.write(unit.segment.text)
                if unit.candidates:
                    st.caption("TM suggestions")
                    for candidate in unit.candidates:
                        st.write(f"{candidate.score:.0%} · {candidate.translation}")
                if unit.quality_issues:
                    st.caption("Quality issues")
                    for issue in unit.quality_issues:
                        st.warning(f"{issue.category}: {issue.message}")
            with right:
                st.caption("Target")
                key = f"translation_{unit.segment.id}"
                unit.translation = st.text_area("Translation", value=unit.translation, key=key, label_visibility="collapsed", height=120)
                action_cols = st.columns(3)
                if action_cols[0].button("Approve", key=f"approve_{unit.segment.id}", type="primary"):
                    store.upsert_memory(unit.segment.text, unit.translation, source_lang, target_lang, domain)
                    unit.status = SegmentStatus.APPROVED
                    approved += 1
                    st.toast("Approved and saved to translation memory.")
                if action_cols[1].button("Reject", key=f"reject_{unit.segment.id}"):
                    unit.status = SegmentStatus.REJECTED
                    rejected += 1
                if action_cols[2].button("Use TM", key=f"tm_{unit.segment.id}", disabled=not unit.candidates):
                    unit.translation = unit.candidates[0].translation
                    st.rerun()

    export_rows = [{"source": unit.segment.text, "target": unit.translation, "status": unit.status.value} for unit in units]
    st.download_button("Download bilingual CSV", to_csv(export_rows), "translation_review.csv", "text/csv")


def render_memory_glossary(store: StudioStore, source_lang: str, target_lang: str, domain: str) -> None:
    st.subheader("Translation Memory")
    memories = store.get_memories(source_lang, target_lang, domain)
    st.caption(f"{len(memories)} approved memory records for {source_lang} to {target_lang}.")
    if memories:
        st.dataframe(memories, use_container_width=True, hide_index=True)

    import_file = st.file_uploader("Import bilingual CSV", type=["csv"], key="memory_import")
    if import_file and st.button("Import approved translations"):
        rows = list(csv.DictReader(io.StringIO(import_file.getvalue().decode("utf-8-sig"))))
        count = store.import_memories(rows)
        st.success(f"Imported {count} approved translation units.")
        st.rerun()

    st.subheader("Glossary")
    with st.form("glossary_form", clear_on_submit=True):
        cols = st.columns([2, 2, 1])
        source_term = cols[0].text_input("Source term")
        target_term = cols[1].text_input("Approved target term")
        case_sensitive = cols[2].checkbox("Case-sensitive")
        notes = st.text_input("Notes")
        if st.form_submit_button("Save glossary term", type="primary") and source_term and target_term:
            store.upsert_glossary(GlossaryEntry(source_term, target_term, source_lang, target_lang, domain, case_sensitive, notes))
            st.success("Glossary entry saved.")
            st.rerun()

    glossary = store.list_glossary(source_lang, target_lang, domain)
    if glossary:
        st.dataframe([asdict(entry) for entry in glossary], use_container_width=True, hide_index=True)


def render_profile_editor(store: StudioStore) -> None:
    st.subheader("Style & Tone Profiles")
    with st.form("profile_form", clear_on_submit=True):
        cols = st.columns(3)
        name = cols[0].text_input("Profile name")
        tone = cols[1].selectbox("Tone", ["Formal", "Official", "Conversational", "Technical", "Social"])
        audience = cols[2].text_input("Audience", "Enterprise")
        rules = st.text_area("Style rules", "Preserve terminology, placeholders, numbers, and formatting.")
        if st.form_submit_button("Save style profile") and name:
            store.add_profile(StyleProfile(name, tone, audience, rules))
            st.success("Style profile saved.")
            st.rerun()


def to_csv(rows: list[dict]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=["source", "target", "status"])
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


if __name__ == "__main__":
    main()
