"""Generate Chapter 1 figures: Conceptual layouts from design.R using Viewport positioning."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from complexheatmap import (
    Viewport,
    convertHeight,
    gpar,
    grid_rect,
    grid_text,
    grobHeight,
    popViewport,
    pushViewport,
    textGrob,
    unit,
)

from _common import reset_state, save_Standalone


def single_heatmap_layout() -> None:
    """Port of design.R::single_heatmap_layout() using explicit Viewport coordinates."""
    th_mm = convertHeight(grobHeight(textGrob("a", gp=gpar(fontsize=12))) * 3.0, "mm", valueOnly=True)
    band = unit(th_mm * 3.2, "mm")

    pushViewport(
        Viewport(
            width=unit(1.0, "npc") - unit(4.0, "mm"),
            height=unit(1.0, "npc") - unit(4.0, "mm"),
            gp=gpar(fontsize=11),
        )
    )
    try:
        # Top (row 1, col 2): bottom-to-top closest to body -> outermost
        pushViewport(
            Viewport(
                x=band,
                y=unit(1.0, "npc"),
                width=unit(1.0, "npc") - band * 2.0,
                height=band,
                just=("left", "top"),
            )
        )
        top_names = ["column annotations", "column names/labels", "dendrogram", "title"]
        for i, lbl in enumerate(top_names, start=1):
            pushViewport(Viewport(y=i / 4.0, height=1.0 / 4.0, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=10.5))
            popViewport()
        popViewport()

        # Bottom (row 3, col 2): top-to-bottom closest to body -> outermost
        pushViewport(
            Viewport(
                x=band,
                y=unit(0.0, "npc"),
                width=unit(1.0, "npc") - band * 2.0,
                height=band,
                just=("left", "bottom"),
            )
        )
        bot_names = list(reversed(top_names))
        for i, lbl in enumerate(bot_names, start=1):
            pushViewport(Viewport(y=i / 4.0, height=1.0 / 4.0, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=10.5))
            popViewport()
        popViewport()

        # Left (row 2, col 1): outermost -> closest to body
        pushViewport(
            Viewport(
                x=unit(0.0, "npc"),
                y=band,
                width=band,
                height=unit(1.0, "npc") - band * 2.0,
                just=("left", "bottom"),
            )
        )
        left_names = ["title", "dendrogram", "row names/labels", "row annotations"]
        for i, lbl in enumerate(left_names, start=1):
            pushViewport(Viewport(x=i / 4.0, width=1.0 / 4.0, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=10.5))
            popViewport()
        popViewport()

        # Right (row 2, col 3): closest to body -> outermost
        pushViewport(
            Viewport(
                x=unit(1.0, "npc"),
                y=band,
                width=band,
                height=unit(1.0, "npc") - band * 2.0,
                just=("right", "bottom"),
            )
        )
        right_names = list(reversed(left_names))
        for i, lbl in enumerate(right_names, start=1):
            pushViewport(Viewport(x=i / 4.0, width=1.0 / 4.0, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=10.5))
            popViewport()
        popViewport()

        # Center (row 2, col 2): 2x2 heatmap body slices
        pushViewport(
            Viewport(
                x=band,
                y=band,
                width=unit(1.0, "npc") - band * 2.0,
                height=unit(1.0, "npc") - band * 2.0,
                just=("left", "bottom"),
            )
        )
        slice_labels = [
            ["matrix, row slice 2\ncolumn slice 1", "matrix, row slice 1\ncolumn slice 1"],
            ["matrix, row slice 2\ncolumn slice 2", "matrix, row slice 1\ncolumn slice 2"],
        ]
        for i in range(1, 3):
            for j in range(1, 3):
                pushViewport(Viewport(x=i / 2.0, y=j / 2.0, width=0.5, height=0.5, just=("right", "top")))
                pushViewport(Viewport(width=unit(1.0, "npc") - unit(4.0, "mm"), height=unit(1.0, "npc") - unit(4.0, "mm")))
                grid_rect(gp=gpar(col="red", fill="#FF000020"))
                grid_text(slice_labels[i - 1][j - 1], gp=gpar(fontsize=10.5))
                popViewport()
                popViewport()
        popViewport()
    finally:
        popViewport()


def heatmap_list_layout(direction: str = "horizontal") -> None:
    """Port of design.R::heatmap_list_layout() using explicit Viewport coordinates."""
    th_mm = convertHeight(grobHeight(textGrob("a", gp=gpar(fontsize=12))) * 3.0, "mm", valueOnly=True)
    band = unit(th_mm * 2.0, "mm")

    pushViewport(
        Viewport(
            width=unit(1.0, "npc") - unit(4.0, "mm"),
            height=unit(1.0, "npc") - unit(4.0, "mm"),
            gp=gpar(fontsize=11),
        )
    )
    try:
        # Top
        pushViewport(
            Viewport(
                x=band,
                y=unit(1.0, "npc"),
                width=unit(1.0, "npc") - band * 2.0,
                height=band,
                just=("left", "top"),
            )
        )
        for i, lbl in enumerate(["title", "legends"], start=1):
            pushViewport(Viewport(y=i / 2.0, height=0.5, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=11))
            popViewport()
        popViewport()

        # Bottom
        pushViewport(
            Viewport(
                x=band,
                y=unit(0.0, "npc"),
                width=unit(1.0, "npc") - band * 2.0,
                height=band,
                just=("left", "bottom"),
            )
        )
        for i, lbl in enumerate(["legends", "title"], start=1):
            pushViewport(Viewport(y=i / 2.0, height=0.5, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=11))
            popViewport()
        popViewport()

        # Left
        pushViewport(
            Viewport(
                x=unit(0.0, "npc"),
                y=band,
                width=band,
                height=unit(1.0, "npc") - band * 2.0,
                just=("left", "bottom"),
            )
        )
        for i, lbl in enumerate(["legends", "title"], start=1):
            pushViewport(Viewport(x=i / 2.0, width=0.5, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=11))
            popViewport()
        popViewport()

        # Right
        pushViewport(
            Viewport(
                x=unit(1.0, "npc"),
                y=band,
                width=band,
                height=unit(1.0, "npc") - band * 2.0,
                just=("right", "bottom"),
            )
        )
        for i, lbl in enumerate(["title", "legends"], start=1):
            pushViewport(Viewport(x=i / 2.0, width=0.5, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=11))
            popViewport()
        popViewport()

        # Center heatmap_list region
        pushViewport(
            Viewport(
                x=band,
                y=band,
                width=unit(1.0, "npc") - band * 2.0,
                height=unit(1.0, "npc") - band * 2.0,
                just=("left", "bottom"),
            )
        )
        if direction == "horizontal":
            items = ["heatmap1", "heatmap2", "row annotation"]
            for i, lbl in enumerate(items, start=1):
                pushViewport(Viewport(x=i / 3.0, width=1.0 / 3.0, just="right"))
                pushViewport(Viewport(width=unit(1.0, "npc") - unit(4.0, "mm"), height=unit(1.0, "npc") - unit(4.0, "mm")))
                grid_rect(gp=gpar(col="green", fill="#00FF0020"))
                grid_text(lbl, rot=90, gp=gpar(fontsize=11))
                popViewport()
                popViewport()
        else:
            items = list(reversed(["heatmap1", "heatmap2", "column annotation"]))
            for i, lbl in enumerate(items, start=1):
                pushViewport(Viewport(y=i / 3.0, height=1.0 / 3.0, just="top"))
                pushViewport(Viewport(width=unit(1.0, "npc") - unit(4.0, "mm"), height=unit(1.0, "npc") - unit(4.0, "mm")))
                grid_rect(gp=gpar(col="red", fill="#FF000020"))
                grid_text(lbl, gp=gpar(fontsize=11))
                popViewport()
                popViewport()
        popViewport()
    finally:
        popViewport()


def main() -> None:
    reset_state()
    save_Standalone(single_heatmap_layout, "ch01_single_layout.png", width=7.5, height=5.5, dpi=140)
    reset_state()

    save_Standalone(lambda: heatmap_list_layout("horizontal"), "ch01_list_horizontal.png", width=4.8, height=4.8, dpi=140)
    reset_state()

    save_Standalone(lambda: heatmap_list_layout("vertical"), "ch01_list_vertical.png", width=4.8, height=4.8, dpi=140)
    reset_state()


if __name__ == "__main__":
    main()
