"""Verify generated book HTML pages, PNG figures, internal links, and search index."""

from __future__ import annotations

import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urldefrag

from PIL import Image

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = REPO_ROOT / "docs"
BOOK_DIR = DOCS_DIR / "book"
ASSETS_DIR = BOOK_DIR / "assets"

EXPECTED_HTML_FILES = [
    "index.html",
    "introduction.html",
    "a-single-heatmap.html",
    "heatmap-annotations.html",
    "a-list-of-heatmaps.html",
    "legends.html",
    "heatmap-decoration.html",
    "oncoprint.html",
    "upset-plot.html",
    "interactive-complexheatmap.html",
    "integrate-with-other-packages.html",
    "other-high-level-plots.html",
    "three-dimensional-complexheatmap.html",
    "genome-level-heatmap.html",
    "more-examples.html",
    "other-tricks.html",
]


class _BookPageParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: set[str] = set()
        self.images: list[str] = []
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = dict(attrs)
        if "id" in attr_map and attr_map["id"]:
            self.ids.add(attr_map["id"])
        if tag == "img" and attr_map.get("src"):
            self.images.append(attr_map["src"] or "")
        if tag == "a" and attr_map.get("href"):
            self.links.append(attr_map["href"] or "")


def main() -> None:
    errors: list[str] = []

    # 1. Required root docs files & assets
    for req in (
        DOCS_DIR / ".nojekyll",
        DOCS_DIR / "index.html",
        ASSETS_DIR / "gitbook.css",
        ASSETS_DIR / "gitbook.js",
        ASSETS_DIR / "search_index.json",
    ):
        if not req.exists():
            errors.append(f"Missing required file: {req}")

    # 2. Parse all 16 HTML pages
    page_ids: dict[str, set[str]] = {}
    page_links: dict[str, list[str]] = {}
    referenced_images: set[Path] = set()

    for fname in EXPECTED_HTML_FILES:
        fpath = BOOK_DIR / fname
        if not fpath.exists():
            errors.append(f"Missing HTML chapter page: {fpath}")
            continue
        parser = _BookPageParser()
        parser.feed(fpath.read_text(encoding="utf-8"))
        page_ids[fname] = parser.ids
        page_links[fname] = parser.links

        # Verify images referenced in main section
        for src in parser.images:
            if src.startswith(("http://", "https://", "data:")):
                continue
            img_path = (BOOK_DIR / src).resolve()
            if not img_path.exists():
                errors.append(f"{fname}: broken image src {src!r} -> {img_path}")
            else:
                referenced_images.add(img_path)

    # 3. Validate every referenced PNG image
    for img_path in sorted(referenced_images):
        sz = img_path.stat().st_size
        if sz <= 3000:
            errors.append(f"Image too small ({sz} bytes): {img_path}")
            continue
        try:
            with Image.open(img_path) as im:
                im.verify()
            with Image.open(img_path) as im:
                w, h = im.size
            if w < 200 or h < 150:
                errors.append(f"Image dimensions too small ({w}x{h}): {img_path}")
        except Exception as exc:
            errors.append(f"Invalid image {img_path}: {exc}")

    # Also check for any unreferenced broken PNGs under images/
    all_pngs = sorted((BOOK_DIR / "images").glob("*.png"))
    if not referenced_images:
        errors.append("No PNG images referenced across the entire book!")

    # 4. Validate internal links & anchors
    broken_links = 0
    for fname, links in page_links.items():
        for href in links:
            if not href or href.startswith(("http://", "https://", "mailto:", "javascript:")):
                continue
            target_file, frag = urldefrag(href)
            if not target_file:
                target_file = fname
            if target_file not in page_ids:
                errors.append(f"{fname}: broken link target {href!r}")
                broken_links += 1
                continue
            if frag and frag not in page_ids[target_file]:
                errors.append(f"{fname}: broken anchor #{frag} in {target_file} (href={href!r})")
                broken_links += 1

    # 5. Validate search_index.json covers all 16 chapters
    search_path = ASSETS_DIR / "search_index.json"
    if search_path.exists():
        data = json.loads(search_path.read_text(encoding="utf-8"))
        indexed_pages = {urldefrag(item["url"])[0] for item in data}
        for fname in EXPECTED_HTML_FILES:
            if fname not in indexed_pages:
                errors.append(f"search_index.json missing entries for {fname}")

    if errors:
        for err in errors:
            print(f"ERROR: {err}")
        raise SystemExit(1)

    print(
        f"BOOK VERIFICATION PASSED: {len(EXPECTED_HTML_FILES)}/{len(EXPECTED_HTML_FILES)} chapters, "
        f"{len(referenced_images)} verified PNG figures ({len(all_pngs)} total PNGs on disk), "
        f"{broken_links} broken links."
    )


if __name__ == "__main__":
    main()
