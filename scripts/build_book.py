"""Compile chapters/*.md into a GitBook-style static HTML site under docs/book/."""

from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import TextLexer, get_lexer_by_name

REPO_ROOT = Path(__file__).resolve().parent.parent
CHAPTERS_DIR = REPO_ROOT / "chapters"
BOOK_DIR = REPO_ROOT / "docs" / "book"
ASSETS_DIR = BOOK_DIR / "assets"

REPO_URL = "https://github.com/zhanxw/ComplexHeatmap-reference-py"
PKG_REPO_URL = "https://github.com/zhanxw/ComplexHeatmap-py"
BOOK_TITLE = "ComplexHeatmap-py Complete Reference"
BOOK_SUBTITLE = "Python Port via ComplexHeatmap-py"

CHAPTERS: list[tuple[str, str, str, str]] = [
    ("", "index.md", "index.html", "About"),
    ("1", "01-introduction.md", "introduction.html", "Introduction"),
    ("2", "02-single-heatmap.md", "a-single-heatmap.html", "A Single Heatmap"),
    ("3", "03-heatmap-annotations.md", "heatmap-annotations.html", "Heatmap Annotations"),
    ("4", "04-a-list-of-heatmaps.md", "a-list-of-heatmaps.html", "A List of Heatmaps"),
    ("5", "05-legends.md", "legends.html", "Legends"),
    ("6", "06-heatmap-decoration.md", "heatmap-decoration.html", "Heatmap Decoration"),
    ("7", "07-oncoprint.md", "oncoprint.html", "OncoPrint"),
    ("8", "08-upset-plot.md", "upset-plot.html", "UpSet plot"),
    ("9", "09-interactive-complexheatmap.md", "interactive-complexheatmap.html", "Interactive ComplexHeatmap"),
    ("10", "10-integrate-with-other-packages.md", "integrate-with-other-packages.html", "Integrate with other packages"),
    ("11", "11-other-high-level-plots.md", "other-high-level-plots.html", "Other High-level Plots"),
    ("12", "12-three-dimensional-complexheatmap.md", "three-dimensional-complexheatmap.html", "Three-dimensional ComplexHeatmap"),
    ("13", "13-genome-level-heatmap.md", "genome-level-heatmap.html", "Genome-level heatmap"),
    ("14", "14-more-examples.md", "more-examples.html", "More Examples"),
    ("15", "15-other-tricks.md", "other-tricks.html", "Other Tricks"),
    ("16", "16-publication-examples.md", "publication-examples.html", "Publication Examples"),
]


@dataclass
class SectionItem:
    level: int
    number: str
    title: str
    anchor: str


@dataclass
class CompiledChapter:
    num: str
    md_name: str
    html_name: str
    title: str
    body_html: str
    sections: list[SectionItem] = field(default_factory=list)
    search_entries: list[dict[str, str]] = field(default_factory=list)


@dataclass
class GalleryItem:
    src: str
    filename: str
    caption: str
    chapter_key: str
    chapter_label: str
    chapter_html: str
    section_label: str
    target_url: str


def _find_pandoc() -> str:
    for candidate in ("/opt/homebrew/bin/pandoc", shutil.which("pandoc")):
        if candidate and Path(candidate).exists():
            return str(candidate)
    raise RuntimeError("pandoc executable not found on system.")


def _pandoc_to_html(md_text: str, pandoc_bin: str) -> str:
    fmt = (
        "markdown+pipe_tables+backtick_code_blocks+fenced_code_blocks"
        "+header_attributes+auto_identifiers-tex_math_dollars"
    )
    cmd = [pandoc_bin, "-f", fmt, "-t", "html5", "--syntax-highlighting=none"]
    proc = subprocess.run(cmd, input=md_text, text=True, capture_output=True)
    if proc.returncode != 0:
        # Fallback for older pandoc without --syntax-highlighting=none
        cmd_fb = [pandoc_bin, "-f", fmt, "-t", "html5", "--no-highlight"]
        proc = subprocess.run(cmd_fb, input=md_text, text=True, capture_output=True, check=True)
    return proc.stdout


_CODE_BLOCK_RE = re.compile(
    r'<pre(?:\s+class="([^"]*)")?><code(?:\s+class="([^"]*)")?>(.*?)</code></pre>',
    re.DOTALL,
)


def _highlight_code_blocks(raw_html: str) -> str:
    formatter = HtmlFormatter(nowrap=True)

    def _repl(match: re.Match[str]) -> str:
        pre_cls = match.group(1) or ""
        code_cls = match.group(2) or ""
        code_raw = html.unescape(match.group(3))

        lang = ""
        for cls in (pre_cls, code_cls):
            for tok in cls.split():
                if tok.startswith("language-"):
                    lang = tok[len("language-"):]
                    break
                elif tok not in ("sourceCode",):
                    lang = tok
                    break
            if lang:
                break

        try:
            lexer = get_lexer_by_name(lang) if lang else TextLexer()
        except Exception:
            lexer = TextLexer()

        highlighted = highlight(code_raw, lexer, formatter).rstrip("\n")
        return (
            '<div class="code-wrapper">'
            '<pre class="codehilite"><code>'
            f"{highlighted}"
            "</code></pre></div>"
        )

    return _CODE_BLOCK_RE.sub(_repl, raw_html)


_HEADING_RE = re.compile(r"<h([1-3])([^>]*)>(.*?)</h\1>", re.DOTALL)
_ATTR_ID_RE = re.compile(r'id="([^"]+)"')
_LEADING_NUM_RE = re.compile(r"^\s*\d+(?:\.\d+)*\.?\s+")
_TAG_STRIP_RE = re.compile(r"<[^>]+>")


def _slugify(text: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9\s-]", "", text).strip().lower()
    s = re.sub(r"[\s_]+", "-", s)
    return s or "section"


def _strip_tags(fragment: str) -> str:
    return html.unescape(_TAG_STRIP_RE.sub("", fragment)).strip()


def _process_headings_and_sections(
    body_html: str, chap_num: str, html_name: str, default_title: str
) -> tuple[str, str, list[SectionItem], list[dict[str, str]]]:
    sections: list[SectionItem] = []
    seen_ids: set[str] = set()
    sec_counters = [0, 0]  # h2, h3
    chapter_title = default_title

    matches = list(_HEADING_RE.finditer(body_html))
    chunks_out: list[str] = []
    search_entries: list[dict[str, str]] = []
    last_end = 0

    heading_records: list[tuple[int, int, int, str, str, str]] = []

    for m in matches:
        level = int(m.group(1))
        attrs = m.group(2) or ""
        inner = m.group(3)
        plain = _strip_tags(inner)
        plain_clean = _LEADING_NUM_RE.sub("", plain).strip() or plain

        id_m = _ATTR_ID_RE.search(attrs)
        anchor = id_m.group(1) if id_m else _slugify(plain_clean)
        base_anchor = anchor
        suffix = 2
        while anchor in seen_ids:
            anchor = f"{base_anchor}-{suffix}"
            suffix += 1
        seen_ids.add(anchor)

        num_label = ""
        if level == 1:
            chapter_title = plain_clean
            if chap_num:
                num_label = chap_num
        elif level == 2:
            sec_counters[0] += 1
            sec_counters[1] = 0
            num_label = f"{chap_num}.{sec_counters[0]}" if chap_num else ""
            sections.append(SectionItem(level=2, number=num_label, title=plain_clean, anchor=anchor))
        elif level == 3:
            if sec_counters[0] == 0:
                sec_counters[0] = 1
            sec_counters[1] += 1
            num_label = (
                f"{chap_num}.{sec_counters[0]}.{sec_counters[1]}"
                if chap_num
                else ""
            )
            sections.append(SectionItem(level=3, number=num_label, title=plain_clean, anchor=anchor))

        # Rebuild heading HTML cleanly without duplicating leading numbers inside inner HTML
        inner_clean = _LEADING_NUM_RE.sub("", inner, count=1) if _LEADING_NUM_RE.match(plain) else inner
        num_span = f'<span class="header-section-number">{html.escape(num_label)}</span> ' if num_label else ""
        anchor_link = f'<a class="anchor-section" aria-label="Anchor link to section" href="#{html.escape(anchor)}">#</a>'
        new_heading = f'<h{level} id="{html.escape(anchor)}">{num_span}{inner_clean}{anchor_link}</h{level}>'

        heading_records.append((m.start(), m.end(), level, num_label, plain_clean, anchor))
        chunks_out.append(body_html[last_end:m.start()])
        chunks_out.append(new_heading)
        last_end = m.end()

    chunks_out.append(body_html[last_end:])
    processed_html = "".join(chunks_out)

    # Build search index entries sliced by headings
    if not heading_records:
        text_plain = re.sub(r"\s+", " ", _strip_tags(body_html)).strip()
        search_entries.append({
            "url": html_name,
            "title": f"{chap_num} {chapter_title}".strip(),
            "text": text_plain[:4000],
        })
    else:
        for idx, rec in enumerate(heading_records):
            start_pos = rec[1]
            end_pos = heading_records[idx + 1][0] if idx + 1 < len(heading_records) else len(body_html)
            seg_text = re.sub(r"\s+", " ", _strip_tags(body_html[start_pos:end_pos])).strip()
            num_label = rec[3]
            plain_clean = rec[4]
            anchor = rec[5]
            full_title = f"{num_label} {plain_clean}".strip()
            url = html_name if rec[2] == 1 else f"{html_name}#{anchor}"
            search_entries.append({
                "url": url,
                "title": full_title,
                "text": seg_text[:4000],
            })

    return processed_html, chapter_title, sections, search_entries


class _GalleryExtractor(HTMLParser):
    """Extract ordered local figures with chapter and section context from compiled chapter HTML."""

    def __init__(self, chap_num: str, chap_title: str, html_name: str) -> None:
        super().__init__()
        self.chap_key = chap_num if chap_num else "0"
        self.chap_label = f"{chap_num}. {chap_title}" if chap_num else chap_title
        self.html_name = html_name
        self.items: list[GalleryItem] = []

        self._in_heading: str | None = None
        self._cur_heading_id: str = ""
        self._heading_buf: list[str] = []
        self._in_anchor_hash = False
        self._current_level = 1
        self._current_anchor: str = ""
        self._current_heading_text: str = self.chap_label

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = dict(attrs)
        if tag in ("h1", "h2", "h3"):
            self._in_heading = tag
            self._cur_heading_id = attr_map.get("id") or ""
            self._heading_buf = []
            self._in_anchor_hash = False
        elif self._in_heading and tag == "a":
            cls = attr_map.get("class") or ""
            if "anchor-section" in cls.split():
                self._in_anchor_hash = True
        elif tag == "img":
            src = (attr_map.get("src") or "").strip()
            if not src or src.startswith(("http://", "https://", "data:")):
                return
            filename = Path(src).name
            raw_alt = attr_map.get("alt") or ""
            caption = re.sub(r"\s+", " ", raw_alt).strip()
            if not caption:
                caption = Path(src).stem.replace("_", " ").replace("-", " ")
            if self._current_level in (2, 3) and self._current_anchor:
                target_url = f"{self.html_name}#{self._current_anchor}"
            else:
                target_url = self.html_name
            self.items.append(
                GalleryItem(
                    src=src,
                    filename=filename,
                    caption=caption,
                    chapter_key=self.chap_key,
                    chapter_label=self.chap_label,
                    chapter_html=self.html_name,
                    section_label=self._current_heading_text or self.chap_label,
                    target_url=target_url,
                )
            )

    def handle_endtag(self, tag: str) -> None:
        if self._in_heading and tag == "a":
            self._in_anchor_hash = False
        elif self._in_heading == tag:
            txt = re.sub(r"\s+", " ", "".join(self._heading_buf)).strip()
            self._current_level = int(tag[1])
            self._current_anchor = self._cur_heading_id
            if txt:
                self._current_heading_text = txt
            self._in_heading = None

    def handle_data(self, data: str) -> None:
        if self._in_heading and not self._in_anchor_hash:
            self._heading_buf.append(data)


def _extract_gallery_items(
    processed_html: str, chap_num: str, chap_title: str, html_name: str
) -> list[GalleryItem]:
    parser = _GalleryExtractor(chap_num, chap_title, html_name)
    parser.feed(processed_html)
    parser.close()
    return parser.items


def _build_gallery_chapter(all_items: list[GalleryItem]) -> CompiledChapter:
    """Build the Figure Gallery CompiledChapter aggregating all book figures."""
    groups: list[tuple[str, str, str, list[GalleryItem]]] = []
    by_key: dict[str, tuple[str, str, list[GalleryItem]]] = {}
    for it in all_items:
        if it.chapter_key not in by_key:
            grp_list: list[GalleryItem] = []
            by_key[it.chapter_key] = (it.chapter_label, it.chapter_html, grp_list)
            groups.append((it.chapter_key, it.chapter_label, it.chapter_html, grp_list))
        by_key[it.chapter_key][2].append(it)

    total = len(all_items)
    pills_html = [
        f'<button type="button" class="gallery-pill active" data-ch-filter="all">All ({total})</button>'
    ]
    for ch_key, ch_label, _ch_html, items in groups:
        short_lbl = f"Ch. {ch_key}" if ch_key != "0" else "About"
        pills_html.append(
            f'<button type="button" class="gallery-pill" data-ch-filter="{html.escape(ch_key)}" '
            f'title="{html.escape(ch_label)}">{html.escape(short_lbl)} ({len(items)})</button>'
        )

    parts: list[str] = [
        '<div id="gallery-root">',
        '<h1 id="figure-gallery">Figure Gallery'
        '<a class="anchor-section" aria-label="Anchor link to section" href="#figure-gallery">#</a></h1>',
        f"<p>Browse all <strong>{total}</strong> executable example figures generated with "
        f'<a href="{html.escape(PKG_REPO_URL)}" target="_blank" rel="noopener"><strong>ComplexHeatmap-py</strong></a>. '
        "Filter by keyword, function name, or chapter; click any thumbnail to zoom cleanly and jump directly "
        "to its full Python recipe and documentation.</p>",
        '<div class="gallery-controls">',
        '  <div class="gallery-search-row">',
        '    <input id="gallery-filter-input" type="search" placeholder="Filter figures by keyword, function, topic, or filename (e.g. oncoprint, upset, anno_textbox, 3d)..." aria-label="Filter figures">',
        f'    <span id="gallery-count-badge">{total} / {total} figures</span>',
        '    <div class="gallery-size-btns" role="group" aria-label="Thumbnail grid density">',
        '      <button type="button" class="gallery-size-btn" data-grid-size="compact">Compact</button>',
        '      <button type="button" class="gallery-size-btn active" data-grid-size="medium">Medium</button>',
        '      <button type="button" class="gallery-size-btn" data-grid-size="large">Large</button>',
        "    </div>",
        "  </div>",
        '  <div class="gallery-chapter-pills">',
        "\n".join(pills_html),
        "  </div>",
        "</div>",
    ]

    sections: list[SectionItem] = []
    search_entries: list[dict[str, str]] = [
        {
            "url": "gallery.html",
            "title": "Figure Gallery",
            "text": f"Interactive gallery of all {total} example figures from ComplexHeatmap-py Complete Reference.",
        }
    ]

    for ch_key, ch_label, ch_html, items in groups:
        anchor_id = f"gallery-ch-{ch_key}"
        sections.append(SectionItem(level=2, number="", title=f"{ch_label} ({len(items)})", anchor=anchor_id))
        search_entries.append(
            {
                "url": f"gallery.html#{anchor_id}",
                "title": f"Gallery: {ch_label}",
                "text": " ".join(f"{it.caption} {it.section_label} {it.filename}" for it in items)[:4000],
            }
        )
        parts.append(f'<section class="gallery-group" data-ch-group="{html.escape(ch_key)}">')
        parts.append(
            f'<h2 id="{html.escape(anchor_id)}">'
            f'<span class="gallery-group-title">{html.escape(ch_label)} '
            f'<span class="gallery-group-count">({len(items)})</span></span>'
            f'<a class="gallery-chapter-open" href="{html.escape(ch_html)}">Read chapter &rarr;</a>'
            f'<a class="anchor-section" aria-label="Anchor link to section" href="#{html.escape(anchor_id)}">#</a>'
            f"</h2>"
        )
        parts.append('<div class="gallery-grid">')
        for it in items:
            haystack = f"{it.caption} {it.section_label} {it.chapter_label} {it.filename}".lower()
            parts.append(
                f'<figure class="gallery-card" data-ch="{html.escape(ch_key)}" data-search="{html.escape(haystack)}">'
                f'<button type="button" class="gallery-thumb-btn" '
                f'data-full-src="{html.escape(it.src)}" '
                f'data-caption="{html.escape(it.caption)}" '
                f'data-chapter="{html.escape(it.chapter_label)}" '
                f'data-section="{html.escape(it.section_label)}" '
                f'data-url="{html.escape(it.target_url)}" '
                f'data-file="{html.escape(it.filename)}" '
                f'aria-label="Inspect {html.escape(it.caption)}">'
                f'<img src="{html.escape(it.src)}" alt="{html.escape(it.caption)}" loading="lazy" decoding="async">'
                f"</button>"
                f'<figcaption class="gallery-card-body">'
                f'<div class="gallery-card-caption">{html.escape(it.caption)}</div>'
                f'<div class="gallery-card-meta">'
                f'<a class="gallery-card-link" href="{html.escape(it.target_url)}">{html.escape(it.section_label)} &rarr;</a>'
                f'<code class="gallery-card-file">{html.escape(it.filename)}</code>'
                f"</div>"
                f"</figcaption>"
                f"</figure>"
            )
        parts.append("</div>")
        parts.append("</section>")

    parts.extend(
        [
            '<p id="gallery-empty-msg" class="gallery-empty" hidden>No figures match your current filter.</p>',
            '<div id="gallery-lightbox" class="gallery-lightbox" role="dialog" aria-modal="true" aria-label="Figure preview" hidden>',
            '  <div class="lightbox-backdrop"></div>',
            '  <div class="lightbox-dialog">',
            '    <button id="lightbox-close" class="lightbox-close" type="button" aria-label="Close preview">&times;</button>',
            '    <button id="lightbox-prev" class="lightbox-nav lightbox-prev" type="button" aria-label="Previous figure">&#8249;</button>',
            '    <div class="lightbox-img-stage"><img id="lightbox-img" src="" alt=""></div>',
            '    <button id="lightbox-next" class="lightbox-nav lightbox-next" type="button" aria-label="Next figure">&#8250;</button>',
            '    <div class="lightbox-footer">',
            '      <div class="lightbox-text">',
            '        <div id="lightbox-counter" class="lightbox-counter"></div>',
            '        <div id="lightbox-caption" class="lightbox-caption"></div>',
            '        <div id="lightbox-meta" class="lightbox-meta"></div>',
            "      </div>",
            '      <div class="lightbox-actions">',
            '        <a id="lightbox-jump" class="lightbox-btn-primary" href="#">View Section &amp; Python Code &rarr;</a>',
            '        <a id="lightbox-raw" class="lightbox-btn-secondary" href="#" target="_blank" rel="noopener">Open PNG &#8599;</a>',
            "      </div>",
            "    </div>",
            "  </div>",
            "</div>",
            "</div>",
        ]
    )

    return CompiledChapter(
        num="",
        md_name="",
        html_name="gallery.html",
        title="Figure Gallery",
        body_html="\n".join(parts),
        sections=sections,
        search_entries=search_entries,
    )


def _build_sidebar_html(chapters: list[CompiledChapter], current_html: str) -> str:
    items: list[str] = []
    for ch in chapters:
        is_active = ch.html_name == current_html
        cls = "chapter active" if is_active else "chapter"
        num_html = f'<span class="header-section-number">{html.escape(ch.num)}</span>' if ch.num else ""
        sub_html = ""
        if ch.sections:
            lis = []
            for sec in ch.sections:
                sec_cls = f"level-{sec.level}"
                sec_num = f'<span class="header-section-number">{html.escape(sec.number)}</span> ' if sec.number else ""
                lis.append(
                    f'<li class="{sec_cls}"><a href="{html.escape(ch.html_name)}#{html.escape(sec.anchor)}">'
                    f"{sec_num}{html.escape(sec.title)}</a></li>"
                )
            sub_html = "<ul>" + "\n".join(lis) + "</ul>"
        items.append(
            f'<li class="{cls}"><a href="{html.escape(ch.html_name)}">'
            f"{num_html}<span>{html.escape(ch.title)}</span></a>{sub_html}</li>"
        )
    return "\n".join(items)


def _render_page(
    ch: CompiledChapter,
    chapters: list[CompiledChapter],
    prev_ch: CompiledChapter | None,
    next_ch: CompiledChapter | None,
) -> str:
    sidebar_items = _build_sidebar_html(chapters, ch.html_name)
    full_ch_title = f"Chapter {ch.num} {ch.title}" if ch.num else ch.title
    page_title = f"{full_ch_title} | {BOOK_TITLE}"

    prev_arrow = (
        f'<a href="{html.escape(prev_ch.html_name)}" class="navigation navigation-prev" '
        f'aria-label="Previous chapter: {html.escape(prev_ch.title)}" title="Previous: {html.escape(prev_ch.title)}">&#8249;</a>'
        if prev_ch
        else ""
    )
    next_arrow = (
        f'<a href="{html.escape(next_ch.html_name)}" class="navigation navigation-next" '
        f'aria-label="Next chapter: {html.escape(next_ch.title)}" title="Next: {html.escape(next_ch.title)}">&#8250;</a>'
        if next_ch
        else ""
    )

    prev_btn = (
        f'<a href="{html.escape(prev_ch.html_name)}">&larr; {html.escape((prev_ch.num + " " if prev_ch.num else "") + prev_ch.title)}</a>'
        if prev_ch
        else "<span></span>"
    )
    next_btn = (
        f'<a href="{html.escape(next_ch.html_name)}">{html.escape((next_ch.num + " " if next_ch.num else "") + next_ch.title)} &rarr;</a>'
        if next_ch
        else "<span></span>"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{html.escape(page_title)}</title>
  <meta name="description" content="Python port of ComplexHeatmap Complete Reference using ComplexHeatmap-py.">
  <link rel="stylesheet" href="assets/gitbook.css">
</head>
<body>
  <div class="book">
    <nav class="book-summary" aria-label="Table of Contents">
      <div class="book-summary-header">
        <a class="book-summary-title" href="index.html">{html.escape(BOOK_TITLE)}</a>
        <span class="book-summary-subtitle">{html.escape(BOOK_SUBTITLE)}</span>
      </div>
      <div class="book-search-box">
        <input id="book-search-input" type="search" placeholder="Search book..." aria-label="Search book">
      </div>
      <ul class="summary">
        {sidebar_items}
        <li class="divider"></li>
      </ul>
      <div class="summary-footer">
        <div>Powered by <a href="{html.escape(PKG_REPO_URL)}" target="_blank" rel="noopener"><strong>ComplexHeatmap-py</strong></a></div>
        <div><a href="{html.escape(PKG_REPO_URL)}" target="_blank" rel="noopener">Library GitHub</a> &middot; <a href="{html.escape(REPO_URL)}" target="_blank" rel="noopener">Book Repo</a> &middot; <a href="https://jokergoo.github.io/ComplexHeatmap-reference/book/" target="_blank" rel="noopener">Original R Book</a></div>
      </div>
    </nav>

    <div class="book-body">
      <header class="book-header">
        <div class="header-left">
          <button id="sidebar-toggle" class="toolbar-btn" type="button" title="Toggle Table of Contents" aria-label="Toggle Table of Contents">&#9776;</button>
          <button id="search-toggle-btn" class="toolbar-btn" type="button" title="Search" aria-label="Search">&#128269;</button>
          <h1 class="header-title"><a href="index.html">{html.escape(BOOK_TITLE)}</a></h1>
        </div>
        <div class="header-right">
          <div class="dropdown">
            <button id="font-settings-btn" class="toolbar-btn" type="button" title="Font and Theme Settings" aria-label="Font Settings">A</button>
            <div id="font-settings-menu" class="dropdown-menu">
              <div class="dropdown-row">
                <button id="font-size-dec" type="button">A-</button>
                <button id="font-size-inc" type="button">A+</button>
              </div>
              <div class="dropdown-row">
                <button type="button" data-font="sans">Sans</button>
                <button type="button" data-font="serif">Serif</button>
              </div>
              <div class="dropdown-row">
                <button type="button" data-theme="white">White</button>
                <button type="button" data-theme="sepia">Sepia</button>
                <button type="button" data-theme="night">Night</button>
              </div>
            </div>
          </div>
          <a class="toolbar-btn" href="gallery.html" title="Browse all example figures">Gallery</a>
          <a class="toolbar-btn" href="{html.escape(PKG_REPO_URL)}" target="_blank" rel="noopener" title="ComplexHeatmap-py Library on GitHub">ComplexHeatmap-py</a>
          <a class="toolbar-btn" href="{html.escape(REPO_URL)}" target="_blank" rel="noopener" title="Book Source on GitHub">Book Repo</a>
        </div>
      </header>

      <div class="page-wrapper">
        <main class="page-inner{' page-inner-wide' if ch.html_name == 'gallery.html' else ''}">
          <div id="search-results" aria-live="polite"></div>
          <section class="normal markdown-section">
            {ch.body_html}
          </section>
          <nav class="chapter-footer-nav" aria-label="Chapter Navigation">
            {prev_btn}
            {next_btn}
          </nav>
        </main>
        {prev_arrow}
        {next_arrow}
      </div>
    </div>
  </div>
  <script src="assets/gitbook.js"></script>
</body>
</html>
"""


def main() -> None:
    pandoc_bin = _find_pandoc()
    BOOK_DIR.mkdir(parents=True, exist_ok=True)
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)

    compiled: list[CompiledChapter] = []
    all_search: list[dict[str, str]] = []
    all_gallery_items: list[GalleryItem] = []

    for num, md_name, html_name, default_title in CHAPTERS:
        md_path = CHAPTERS_DIR / md_name
        if not md_path.exists():
            raise FileNotFoundError(f"Missing chapter markdown: {md_path}")
        md_text = md_path.read_text(encoding="utf-8")
        raw_html = _pandoc_to_html(md_text, pandoc_bin)
        hl_html = _highlight_code_blocks(raw_html)
        body_html, title, sections, search_entries = _process_headings_and_sections(
            hl_html, num, html_name, default_title
        )
        all_gallery_items.extend(_extract_gallery_items(body_html, num, title, html_name))
        compiled.append(
            CompiledChapter(
                num=num,
                md_name=md_name,
                html_name=html_name,
                title=title,
                body_html=body_html,
                sections=sections,
                search_entries=search_entries,
            )
        )
        all_search.extend(search_entries)

    gallery_ch = _build_gallery_chapter(all_gallery_items)
    compiled.append(gallery_ch)
    all_search.extend(gallery_ch.search_entries)

    for idx, ch in enumerate(compiled):
        prev_ch = compiled[idx - 1] if idx > 0 else None
        next_ch = compiled[idx + 1] if idx + 1 < len(compiled) else None
        page_html = _render_page(ch, compiled, prev_ch, next_ch)
        out_path = BOOK_DIR / ch.html_name
        out_path.write_text(page_html, encoding="utf-8")

    search_path = ASSETS_DIR / "search_index.json"
    search_path.write_text(json.dumps(all_search, indent=2), encoding="utf-8")
    print(
        f"Built {len(compiled)} book HTML pages ({len(all_gallery_items)} gallery figures) "
        f"and {len(all_search)} search entries in {BOOK_DIR}"
    )


if __name__ == "__main__":
    main()
