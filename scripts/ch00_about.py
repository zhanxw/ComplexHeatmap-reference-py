"""Generate Chapter 0 figure: ComplexHeatmap-py hex-badge logo."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import complexheatmap as ch
from complexheatmap import (
    Heatmap,
    Viewport,
    add_transparency,
    draw,
    gpar,
    grid_draw,
    grid_grabExpr,
    grid_polygon,
    grid_text,
    popViewport,
    pushViewport,
    rand_color,
    unit,
)

from _common import reset_state, save_Standalone


def render_logo() -> None:
    """Port of ComplexHeatmap/logo.R adapted from showcase_gallery.py."""
    lines = [
        "oooooo                  oooooo  oo          oo  oo",
        "oooooo                  oooooo  oo          oo  oo",
        "oo      oooooo  oo  oo  oo  oo  oo  oooooo  oooooo",
        "oo      oooooo  oooooo  oo  oo  oo  oooooo  oooooo",
        "oo      oo  oo  oooooo  oooo    oo  oo  oo    oo  ",
        "oo      oo  oo  oooooo  oooo    oo  oo  oo    oo  ",
        "oooooo  oooooo  oo  oo  oo      oo  oooo    oo  oo",
        "oooooo  oooooo  oo  oo  oo      oo  oooo    oo  oo",
    ]
    lines2 = [
        "oo  oo                    oo                    oooooo",
        "oo  oo                    oo                    oooooo",
        "oo  oo  oooooo    oooo  oooooo  oo  oo  oooooo  oo  oo",
        "oo  oo  oooooo    oooo  oooooo  oooooo  oooooo  oo  oo",
        "oooooo  oo  oo  oo  oo    oo    oooooo  oo  oo  oooo  ",
        "oooooo  oo  oo  oo  oo    oo    oooooo  oo  oo  oooo  ",
        "oo  oo  oooo    oooooo    oo    oo  oo  oooo    oo    ",
        "oo  oo  oooo    oooooo    oo    oo  oo  oooo    oo    ",
    ]

    mat1 = np.zeros((8, max(len(ln) for ln in lines)), dtype=int)
    for i, ln in enumerate(lines):
        for j, ch_c in enumerate(ln):
            if not ch_c.isspace():
                mat1[i, j] = 1
    mat1 = np.hstack([np.zeros((8, 2), dtype=int), mat1])

    mat2 = np.zeros((8, max(len(ln) for ln in lines2)), dtype=int)
    for i, ln in enumerate(lines2):
        for j, ch_c in enumerate(ln):
            if not ch_c.isspace():
                mat2[i, j] = 1

    if mat1.shape[1] < mat2.shape[1]:
        mat1 = np.hstack([mat1, np.zeros((8, mat2.shape[1] - mat1.shape[1]), dtype=int)])
    elif mat2.shape[1] < mat1.shape[1]:
        mat2 = np.hstack([mat2, np.zeros((8, mat1.shape[1] - mat2.shape[1]), dtype=int)])

    mat = np.vstack([mat1, np.zeros((2, mat1.shape[1]), dtype=int), mat2])
    mat = np.vstack([np.zeros((16, mat.shape[1]), dtype=int), mat, np.zeros((16, mat.shape[1]), dtype=int)])
    mat = np.hstack([np.zeros((mat.shape[0], 6), dtype=int), mat, np.zeros((mat.shape[0], 6), dtype=int)])
    mat[12:14, 5:-5] = 1
    mat[-14:-12, 5:-5] = 1

    rng = np.random.default_rng(123)
    bg_palette = rand_color(32, luminosity="light", transparency=0.82)

    def col_fun(x: np.ndarray) -> np.ndarray:
        arr = np.asarray(x, dtype=float)
        cols = np.empty(arr.shape, dtype=object)
        for idx, val in np.ndenumerate(arr):
            if val == 1:
                cols[idx] = add_transparency("#4DAF4A", float(rng.uniform(0.0, 0.3)))
            else:
                cols[idx] = bg_palette[int(rng.integers(0, len(bg_palette)))]
        return cols

    col_fun.breaks = [0.0, 1.0]  # type: ignore[attr-defined]

    ht = Heatmap(
        mat,
        name="logo",
        col=col_fun,
        rect_gp=gpar(col="white", lwd=0.4),
        cluster_rows=False,
        cluster_columns=False,
        show_heatmap_legend=False,
    )
    g = grid_grabExpr(
        lambda: draw(ht, padding=unit([0, 0, 0, 0], "mm"), newpage=False),
        width_in=6.0,
        height_in=4.5,
    )

    def _draw_hex() -> None:
        pushViewport(
            Viewport(
                xscale=(0.0, 2.0),
                yscale=(0.0, 2.0),
                width=unit(0.92, "snpc"),
                height=unit(0.92, "snpc"),
            )
        )
        ang = np.arange(6) * (math.pi / 3.0) + math.pi / 6.0
        hx = np.cos(ang) + 1.0
        hy = np.sin(ang) + 1.0
        grid_polygon(hx, hy, default_units="native", gp=gpar(fill="#F9FFF9", col="#4DAF4A", lwd=5))

        pushViewport(Viewport(x=0.5, y=0.525, width=0.64, height=0.64 * mat.shape[0] / mat.shape[1]))
        grid_draw(g)
        popViewport()

        grid_text(
            "Python Port  |  v0.2.0",
            x=unit(1.0, "native"),
            y=unit(0.30, "native"),
            gp=gpar(fontsize=10.5, fontface="bold", col="#2B6E2A"),
        )
        grid_polygon(hx, hy, default_units="native", gp=gpar(fill="transparent", col="#4DAF4A", lwd=7))
        popViewport()

    save_Standalone(_draw_hex, "ch00_logo.png", width=6.0, height=6.0, dpi=140)
    reset_state()


def main() -> None:
    reset_state()
    render_logo()


if __name__ == "__main__":
    main()
