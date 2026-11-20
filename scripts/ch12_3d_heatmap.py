"""Generate figures for Chapter 12: Three-dimensional ComplexHeatmap."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

import complexheatmap as ch
from _common import reset_state, save_Standalone, save_ht


def render_motivation() -> None:
    reset_state()
    rng = np.random.default_rng(123)
    mat = pd.DataFrame(rng.normal(size=(50, 10)), columns=list("abcdefghij"))

    # 1. densityHeatmap
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    dh = ch.densityHeatmap(mat, title="Density heatmap")
    save_ht(dh, "ch12_motivation_density.png", width=7.2, height=5.0, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 2. 2D frequencyHeatmap
    reset_state()
    fh = ch.frequencyHeatmap(mat, title="Frequency heatmap (2D)")
    save_ht(fh, "ch12_motivation_freq.png", width=7.0, height=5.0, padding=ch.unit([4, 10, 4, 4], "mm"))

    # 3. 3D frequencyHeatmap
    reset_state()
    ch.ht_opt.TITLE_PADDING = ch.unit([10, 2], "mm")
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(8, "mm")
    fh3d = ch.frequencyHeatmap(
        mat,
        use_3d=True,
        title="Frequency heatmap (use_3d=True)",
        top_annotation=ch.HeatmapAnnotation(
            _spacer=ch.anno_empty(border=False, height=ch.unit(6, "mm")),
            show_annotation_name=False,
        ),
        bar_max_length=ch.unit(8, "mm"),
    )
    save_ht(fh3d, "ch12_motivation_freq3d.png", width=7.5, height=5.5, padding=ch.unit([4, 10, 6, 6], "mm"))


def render_bar3d_demos() -> None:
    # 1. Single 3D bar
    reset_state()

    def _draw_single() -> None:
        ch.bar3D(x=0.5, y=0.5, w=0.2, h=0.2, l=ch.unit(1, "cm"), theta=60)

    save_Standalone(_draw_single, "ch12_bar3d_single.png", width=3.5, height=3.0)

    # 2. Projection diagram (w, h, l, theta)
    reset_state()

    def _draw_diagram() -> None:
        ch.bar3D(x=0.45, y=0.42, w=0.22, h=0.22, l=ch.unit(1.4, "cm"), theta=60)
        ch.grid_text("w", 0.45, 0.29, just="top", gp=ch.gpar(fontfamily="serif", fontface="italic", fontsize=12))
        ch.grid_text("h", 0.32, 0.42, just="right", gp=ch.gpar(fontfamily="serif", fontface="italic", fontsize=12))
        ch.grid_text(
            "l",
            ch.unit(0.56, "npc") + ch.unit(0.7, "cm") * 0.5 + ch.unit(3.5, "mm"),
            ch.unit(0.31, "npc") + ch.unit(0.7, "cm") * math.sqrt(3) * 0.5,
            gp=ch.gpar(fontfamily="serif", fontface="italic", fontsize=12),
        )
        ch.grid_segments(
            ch.unit(0.34, "npc") + ch.unit(4.5, "mm"),
            ch.unit(0.31, "npc"),
            ch.unit(0.34, "npc") + ch.unit(4.5, "mm") * 0.5,
            ch.unit(0.31, "npc") + ch.unit(4.5, "mm") * math.sqrt(3) * 0.5,
        )
        ch.grid_text(
            "\u03b8",
            ch.unit(0.34, "npc") + ch.unit(7.5, "mm") * math.sqrt(3) * 0.5,
            ch.unit(0.31, "npc") + ch.unit(7.5, "mm") * 0.5,
            gp=ch.gpar(fontfamily="serif", fontface="italic", fontsize=12),
        )

    save_Standalone(_draw_diagram, "ch12_bar3d_diagram.png", width=4.0, height=3.5)

    # 3. Multiple fill colors (add_luminance shading across 3 faces)
    reset_state()

    def _draw_colors() -> None:
        ch.bar3D(
            x=np.linspace(0.2, 0.8, 4),
            y=0.35,
            w=ch.unit(6, "mm"),
            h=ch.unit(6, "mm"),
            l=ch.unit(1.1, "cm"),
            fill=["red", "green", "blue", "purple"],
        )

    save_Standalone(_draw_colors, "ch12_bar3d_colors.png", width=6.0, height=2.2)

    # 4. Different projection angles theta
    reset_state()

    def _draw_thetas() -> None:
        xs = np.linspace(0.18, 0.78, 4)
        thetas = [20.0, 40.0, 60.0, 80.0]
        for x_pos, th in zip(xs, thetas):
            ch.bar3D(
                x=float(x_pos),
                y=0.32,
                w=ch.unit(6, "mm"),
                h=ch.unit(6, "mm"),
                l=ch.unit(1.1, "cm"),
                theta=th,
                fill="#4DAF4A",
            )
            ch.grid_text(f"theta={int(th)}\u00b0", x=float(x_pos), y=0.14, just="top", gp=ch.gpar(fontsize=9.5))

    save_Standalone(_draw_thetas, "ch12_bar3d_theta.png", width=6.4, height=2.4)


def render_heatmap3d_examples() -> None:
    rng = np.random.default_rng(7)
    mat = pd.DataFrame(
        rng.uniform(0.0, 1.0, size=(10, 10)),
        index=list("ABCDEFGHIJ"),
        columns=list("abcdefghij"),
    )

    # 1. Default Heatmap3D (showing 3D bars with default padding)
    reset_state()
    h1 = ch.Heatmap3D(mat, name="mat", column_title="This is a 3D heatmap")
    save_ht(h1, "ch12_heatmap3d_basic.png", width=7.5, height=6.0)

    # 2. Adjusted padding + _spacer for 3D bar overhang
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(8, "mm")
    ch.ht_opt.TITLE_PADDING = ch.unit([10, 2], "mm")
    h2 = ch.Heatmap3D(
        mat,
        name="mat",
        column_title="3D heatmap with adjusted padding",
        top_annotation=ch.HeatmapAnnotation(
            _spacer=ch.anno_empty(border=False, height=ch.unit(7, "mm")),
            show_annotation_name=False,
        ),
        bar_max_length=ch.unit(8, "mm"),
    )
    save_ht(h2, "ch12_heatmap3d_padded.png", width=8.0, height=6.4, padding=ch.unit([4, 4, 6, 6], "mm"))

    # 3. Split Heatmap3D with bar geometry controls
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(8, "mm")
    ch.ht_opt.TITLE_PADDING = ch.unit([8, 2], "mm")
    h3 = ch.Heatmap3D(
        mat,
        name="score",
        col=ch.colorRamp2([0, 0.5, 1], ["#2B8CBE", "#FEE08B", "#D7301F"]),
        row_split=["Group 1"] * 5 + ["Group 2"] * 5,
        column_split=["Cond A"] * 5 + ["Cond B"] * 5,
        bar_rel_width=0.75,
        bar_rel_height=0.75,
        bar_max_length=ch.unit(7.5, "mm"),
        top_annotation=ch.HeatmapAnnotation(
            _spacer=ch.anno_empty(border=False, height=ch.unit(6, "mm")),
            show_annotation_name=False,
        ),
        row_gap=ch.unit(3, "mm"),
        column_gap=ch.unit(3, "mm"),
    )
    save_ht(
        h3,
        "ch12_heatmap3d_split.png",
        width=8.5,
        height=6.8,
        column_title="Split 3D Heatmap",
        padding=ch.unit([4, 4, 6, 6], "mm"),
    )


def render_heatmap3d_measles() -> None:
    reset_state()
    mat = ch.load_dataset("measles")
    year_text = [str(c) if int(c) % 10 == 0 else "" for c in mat.columns]
    ha_column = ch.HeatmapAnnotation(
        year=ch.anno_text(year_text, rot=0, location=ch.unit(1, "npc"), just="top", gp=ch.gpar(fontsize=9)),
        show_annotation_name=False,
    )
    col_fun = ch.colorRamp2([0, 800, 1000, 127000], ["white", "cornflowerblue", "yellow", "red"])
    with ch.ht_opt.local(TITLE_PADDING=ch.unit([14, 2], "mm"), HEATMAP_LEGEND_PADDING=ch.unit(7, "mm")):
        h_measles = ch.Heatmap3D(
            mat,
            name="cases",
            col=col_fun,
            cluster_columns=False,
            show_row_dend=False,
            show_column_names=False,
            row_names_side="left",
            row_names_gp=ch.gpar(fontsize=7.5),
            column_title="Measles cases in US states 1930-2001\nVaccine introduced 1961",
            top_annotation=ch.HeatmapAnnotation(
                _spacer=ch.anno_empty(border=False, height=ch.unit(8, "mm")),
                show_annotation_name=False,
            ),
            bottom_annotation=ha_column,
            heatmap_legend_param={"at": [0, 5e4, 1e5, 1.5e5], "labels": ["0", "50k", "100k", "150k"]},
            bar_rel_width=1.0,
            bar_rel_height=1.0,
            bar_max_length=ch.unit(1.5, "cm"),
        )
        save_ht(h_measles, "ch12_heatmap3d_measles.png", width=11.0, height=10.0, padding=ch.unit([4, 4, 6, 6], "mm"))


def main() -> None:
    render_motivation()
    render_bar3d_demos()
    render_heatmap3d_examples()
    render_heatmap3d_measles()
    reset_state()
    print("Chapter 12 figures generated successfully.")


if __name__ == "__main__":
    main()
