"""Generate all Chapter 4 figures: A List of Heatmaps."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import complexheatmap as ch
from complexheatmap import (
    Heatmap,
    HeatmapAnnotation,
    Viewport,
    anno_barplot,
    anno_points,
    anno_text,
    colorRamp2,
    gpar,
    grid_draw,
    grid_grabExpr,
    ht_opt,
    popViewport,
    pushViewport,
    rowAnnotation,
    unit,
    v,
)

from _common import reset_state, save_Standalone, save_ht


def _make_base_data():
    rng = np.random.default_rng(123)
    mat1 = pd.DataFrame(
        np.vstack([rng.normal(2, 1, (8, 10)), rng.normal(-2, 1, (4, 10))]),
        index=[f"R{i}" for i in range(1, 13)],
        columns=[f"C{j}" for j in range(1, 11)],
    )
    mat2 = pd.DataFrame(
        np.vstack([rng.uniform(1, 3, (6, 10)), rng.uniform(0, 2, (6, 10))]),
        index=[f"R{i}" for i in range(1, 13)],
        columns=[f"C{j}" for j in range(1, 11)],
    )
    le = rng.choice(["a", "b", "c"], 12)
    ind = rng.permutation(12)
    mat1 = mat1.iloc[ind]
    mat2 = mat2.iloc[ind]
    le = le[ind]

    col_rnorm = colorRamp2([-3, 0, 3], ["green", "white", "red"])
    col_runif = colorRamp2([0, 3], ["white", "orange"])
    col_letters = {"a": "pink", "b": "purple", "c": "blue"}
    return mat1, mat2, le, col_rnorm, col_runif, col_letters


def render_intro_and_4_1_to_4_5() -> None:
    mat1, mat2, le, col_rnorm, col_runif, col_letters = _make_base_data()

    ht1 = Heatmap(mat1, name="rnorm")
    ht2 = Heatmap(mat2, name="runif")
    ht3 = Heatmap(le, name="letters")
    save_ht(ht1 + ht2 + ht3, "ch04_intro_concat.png", width=8.8, height=5.0)
    reset_state()

    # 4.1 Titles
    ht1_t = Heatmap(mat1, name="rnorm", col=col_rnorm, row_title="Heatmap 1", column_title="Heatmap 1")
    ht2_t = Heatmap(mat2, name="runif", col=col_runif, row_title="Heatmap 2", column_title="Heatmap 2")
    ht3_t = Heatmap(le, name="letters", col=col_letters)
    save_ht(
        ht1_t + ht2_t + ht3_t,
        "ch04_titles.png",
        width=9.0,
        height=5.4,
        row_title="Three heatmaps, row title",
        row_title_gp=gpar(col="red"),
        column_title="Three heatmaps, column title",
        column_title_gp=gpar(fontsize=16),
    )
    reset_state()

    # 4.2 Size of heatmaps
    ht1_c = Heatmap(mat1, name="rnorm", col=col_rnorm)
    ht2_s = Heatmap(mat2, name="runif", col=col_runif, width=unit(4, "cm"))
    ht3_s = Heatmap(le, name="letters", col=col_letters, width=unit(5, "mm"))
    save_ht(ht1_c + ht2_s + ht3_s, "ch04_size_partial.png", width=8.8, height=5.0)
    reset_state()

    ht1_a = Heatmap(mat1, name="rnorm", col=col_rnorm, width=unit(4, "cm"))
    ht2_a = Heatmap(mat2, name="runif", col=col_runif, width=unit(6, "cm"))
    ht3_a = Heatmap(le, name="letters", col=col_letters, width=unit(1, "cm"))
    save_ht(ht1_a + ht2_a + ht3_a, "ch04_size_all_abs.png", width=8.8, height=5.0)
    reset_state()

    ht1_n = Heatmap(mat1, name="rnorm", col=col_rnorm, width=unit(6, "null"))
    ht2_n = Heatmap(mat2, name="runif", col=col_runif, width=unit(4, "null"))
    ht3_n = Heatmap(le, name="letters", col=col_letters, width=unit(1, "null"))
    save_ht(ht1_n + ht2_n + ht3_n, "ch04_size_null.png", width=8.8, height=5.0)
    reset_state()

    # 4.3 Gap between heatmaps
    ht2_c = Heatmap(mat2, name="runif", col=col_runif)
    ht3_c = Heatmap(le, name="letters", col=col_letters)
    hl_gap = ht1_c + ht2_c + ht3_c
    save_ht(hl_gap, "ch04_gap_single.png", width=9.0, height=5.0, ht_gap=unit(1, "cm"))
    reset_state()
    save_ht(hl_gap, "ch04_gap_vector.png", width=9.0, height=5.0, ht_gap=unit([3, 10], "mm"))
    reset_state()

    # 4.4 Automatic adjustment to main heatmap
    ht1_km = Heatmap(mat1, name="rnorm", col=col_rnorm, row_km=2)
    save_ht(ht2_c + ht1_km + ht3_c, "ch04_auto_adjust_first.png", width=8.8, height=5.0)
    reset_state()

    save_ht(ht2_c + ht1_km + ht3_c, "ch04_main_heatmap_name.png", width=8.8, height=5.0, main_heatmap="rnorm")
    reset_state()

    save_ht(
        ht2_c + ht1_km + ht3_c,
        "ch04_main_heatmap_sides.png",
        width=8.8,
        height=5.0,
        main_heatmap="rnorm",
        row_dend_side="right",
        row_sub_title_side="left",
    )
    reset_state()

    ht1_nc = Heatmap(mat1, name="rnorm", col=col_rnorm, cluster_rows=False)
    save_ht(ht1_nc + ht2_c + ht3_c, "ch04_main_no_cluster.png", width=8.8, height=5.0)
    reset_state()

    save_ht(ht1_c + ht2_c + ht3_c, "ch04_auto_adjust_false.png", width=9.2, height=5.0, auto_adjust=False)
    reset_state()

    # 4.5 Control main heatmap in draw()
    ht1_ov = Heatmap(mat1, name="rnorm", col=col_rnorm, row_km=2, cluster_rows=False)
    save_ht(
        ht1_ov + ht2_c + ht3_c,
        "ch04_draw_override.png",
        width=8.8,
        height=5.2,
        row_km=1,
        row_split=le,
        cluster_rows=True,
    )
    reset_state()


def render_4_6_annotations_adjusted() -> None:
    mat1, mat2, le, col_rnorm, col_runif, col_letters = _make_base_data()
    ht2_c = Heatmap(mat2, name="runif", col=col_runif)
    ht3_c = Heatmap(le, name="letters", col=col_letters)

    ha1 = HeatmapAnnotation(foo1=np.arange(1, 11), annotation_name_side="left")
    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, top_annotation=ha1) + ht2_c + ht3_c,
        "ch04_anno_adj_dend.png",
        width=8.8,
        height=5.4,
    )
    reset_state()

    ha1_c = HeatmapAnnotation(foo1=np.arange(1, 11), bar1=anno_points(np.arange(1, 11)), annotation_name_side="left")
    ha2_c = HeatmapAnnotation(bar2=anno_barplot(np.arange(1, 11)))
    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, top_annotation=ha1_c)
        + Heatmap(mat2, name="runif", col=col_runif, top_annotation=ha2_c)
        + ht3_c,
        "ch04_anno_adj_complex.png",
        width=9.0,
        height=5.5,
        ht_gap=unit([6, 2], "mm"),
    )
    reset_state()

    ha2_tall = HeatmapAnnotation(bar2=anno_barplot(np.arange(1, 11), height=unit(2, "cm")))
    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, top_annotation=ha1)
        + Heatmap(mat2, name="runif", col=col_runif, top_annotation=ha2_tall)
        + ht3_c,
        "ch04_anno_adj_simple_vs_complex.png",
        width=9.0,
        height=5.6,
        ht_gap=unit([6, 2], "mm"),
    )
    reset_state()

    ha2_mat = HeatmapAnnotation(bar2=pd.DataFrame({"b1": np.arange(1, 11), "b2": np.arange(11, 21), "b3": np.arange(21, 31)}))
    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, top_annotation=ha1)
        + Heatmap(mat2, name="runif", col=col_runif, top_annotation=ha2_mat)
        + ht3_c,
        "ch04_anno_adj_unequal_simple.png",
        width=9.0,
        height=5.5,
    )
    reset_state()

    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, bottom_annotation=ha1_c) + ht2_c + ht3_c,
        "ch04_anno_adj_bottom.png",
        width=9.0,
        height=5.5,
    )
    reset_state()


def render_4_7_to_4_10() -> None:
    mat1, mat2, le, col_rnorm, col_runif, col_letters = _make_base_data()

    # 4.7 Concatenate with annotations
    ha_r = rowAnnotation(foo=np.arange(1, 13), bar=anno_barplot(np.arange(1, 13), width=unit(4, "cm")))
    ht1 = Heatmap(mat1, name="rnorm", col=col_rnorm, row_km=2)
    save_ht(ht1 + ha_r, "ch04_concat_row_anno.png", width=8.0, height=5.2)
    reset_state()

    save_ht(
        ht1 + rowAnnotation(foo=np.arange(1, 13)) + rowAnnotation(bar=anno_barplot(np.arange(1, 13), width=unit(4, "cm"))),
        "ch04_concat_sep_anno.png",
        width=8.0,
        height=5.2,
    )
    reset_state()

    save_ht(
        Heatmap(mat1, name="rnorm", col=col_rnorm, row_km=2, right_annotation=ha_r),
        "ch04_concat_right_anno.png",
        width=8.0,
        height=5.2,
    )
    reset_state()

    save_ht(
        ht1 + ha_r + rowAnnotation(rn=anno_text(mat1.index, location=unit(0, "npc"), just="left")),
        "ch04_concat_text_anno.png",
        width=8.4,
        height=5.2,
    )
    reset_state()

    save_ht(
        rowAnnotation(foo=np.arange(1, 13))
        + ht1
        + rowAnnotation(bar=anno_barplot(np.arange(1, 13), width=unit(4, "cm")))
        + Heatmap(mat2, name="runif", col=col_runif),
        "ch04_concat_arbitrary.png",
        width=9.6,
        height=5.2,
    )
    reset_state()

    # 4.8 Concatenate only annotations
    anno_list = rowAnnotation(foo=np.arange(1, 13)) + rowAnnotation(bar=anno_barplot(np.arange(1, 13), width=unit(4, "cm")))
    save_ht(anno_list, "ch04_anno_only.png", width=6.0, height=5.2, padding=unit([12, 2, 2, 2], "mm"))
    reset_state()

    save_ht(
        rowAnnotation(bar=anno_barplot(np.arange(1, 13), width=unit(4, "cm"))) + None,
        "ch04_anno_only_single.png",
        width=5.6,
        height=5.2,
        padding=unit([12, 2, 2, 2], "mm"),
    )
    reset_state()

    save_ht(
        Heatmap(np.empty((12, 0)), cluster_rows=False, row_dend_reorder=False, show_heatmap_legend=False, width=unit(0, "mm")) + anno_list,
        "ch04_anno_only_split.png",
        width=6.0,
        height=5.4,
        row_split=np.repeat(["A", "B"], 6),
        padding=unit([12, 2, 2, 2], "mm"),
    )
    reset_state()

    # 4.9 Vertical concatenation
    mat1t = mat1.T
    mat2t = mat2.T
    ht1_v = Heatmap(mat1t, name="rnorm", col=col_rnorm, row_title="rnorm")
    ht2_v = Heatmap(mat2t, name="runif", col=col_runif, row_title="runif")
    ht3_v = Heatmap(pd.DataFrame([le], index=["letters"]), name="letters", col=col_letters)
    save_ht(ht1_v %v% ht2_v %v% ht3_v, "ch04_vertical_basic.png", width=6.2, height=6.8)
    reset_state()

    save_ht(ht1_v %v% ht2_v %v% ht3_v, "ch04_vertical_col_km.png", width=6.2, height=6.8, column_km=2)
    reset_state()

    ha_v = HeatmapAnnotation(foo=anno_barplot(np.arange(1, 13), height=unit(2, "cm")))
    save_ht(ht1_v %v% ha_v %v% ht2_v %v% ht3_v, "ch04_vertical_with_anno.png", width=6.2, height=7.6, column_km=2)
    reset_state()

    ht1_vk = Heatmap(mat1t, name="rnorm", col=col_rnorm, row_km=2)
    ht2_vk = Heatmap(mat2t, name="runif", col=col_runif, row_km=2)
    save_ht(ht1_vk %v% ha_v %v% ht2_vk %v% ht3_v, "ch04_vertical_row_km.png", width=6.2, height=7.8, column_km=2)
    reset_state()

    ht1_va = Heatmap(
        mat1t,
        name="rnorm",
        col=col_rnorm,
        row_km=2,
        left_annotation=rowAnnotation(foo1=np.arange(1, 11), bar1=anno_barplot(np.arange(1, 11))),
    )
    ha_vr = HeatmapAnnotation(foo=anno_barplot(np.arange(1, 13), height=unit(2, "cm"), axis_param={"side": "right"}))
    ht2_va = Heatmap(mat2t, name="runif", col=col_runif, row_km=2, left_annotation=rowAnnotation(foo2=np.arange(1, 11)))
    save_ht(ht1_va %v% ha_vr %v% ht2_va %v% ht3_v, "ch04_vertical_row_and_col_anno.png", width=6.6, height=7.8, column_km=2)
    reset_state()

    # 4.10 Subset heatmap list
    ht1_s = Heatmap(
        mat1,
        name="rnorm",
        col=col_rnorm,
        left_annotation=rowAnnotation(foo1=np.arange(1, 13), bar1=anno_points(np.arange(1, 13))),
    )
    ht2_s = Heatmap(mat2, name="runif", col=col_runif)
    ha_s = rowAnnotation(foo2=anno_barplot(np.arange(1, 13)), bar2=np.arange(12, 0, -1))
    ht_list_s = ht1_s + ht2_s + ha_s
    save_ht(ht_list_s[:6, ["rnorm", "bar2"]], "ch04_subset.png", width=6.8, height=4.4)
    reset_state()


def render_4_13_to_4_15() -> None:
    mat1, mat2, le, col_rnorm, col_runif, col_letters = _make_base_data()

    # 4.13 Global options
    with ht_opt.local(
        heatmap_column_names_gp=gpar(fontface="italic"),
        heatmap_column_title_gp=gpar(fontsize=10),
        legend_border="black",
        heatmap_border=True,
        annotation_border=True,
    ):
        ht1_o = Heatmap(mat1, name="ht1", column_title="Heatmap 1", top_annotation=HeatmapAnnotation(foo=np.arange(1, 11)))
        ht2_o = Heatmap(mat2, name="ht2", column_title="Heatmap 2", top_annotation=HeatmapAnnotation(bar=np.arange(1, 11)))
        save_ht(ht1_o + ht2_o, "ch04_ht_opt.png", width=8.6, height=5.4)
    reset_state()

    # 4.14 Annotation blank space & extension
    rng = np.random.default_rng(123)
    m = rng.standard_normal((10, 10))
    ht_ext_right = Heatmap(
        m,
        name="mat",
        top_annotation=HeatmapAnnotation(foo=anno_points(np.arange(1, 11))),
        show_row_dend=False,
    )
    ht_ext_left = Heatmap(
        m,
        name="mat",
        top_annotation=HeatmapAnnotation(foo=anno_points(np.arange(1, 11)), annotation_name_side="left"),
        row_title="fooooooooooo",
        show_row_dend=False,
    )
    g_l = grid_grabExpr(
        lambda: ch.draw(ht_ext_right, row_title="fooooooooooo", column_title="Default (annotation_name_side='right')", newpage=False),
        width_in=5.0,
        height_in=4.4,
    )
    g_r = grid_grabExpr(
        lambda: ch.draw(ht_ext_left, column_title="annotation_name_side='left'", newpage=False),
        width_in=5.0,
        height_in=4.4,
    )

    def _draw_ext_pair() -> None:
        pushViewport(Viewport(x=0.25, width=0.5))
        grid_draw(g_l)
        popViewport()
        pushViewport(Viewport(x=0.75, width=0.5))
        grid_draw(g_r)
        popViewport()

    save_Standalone(_draw_ext_pair, "ch04_adj_ext_comparison.png", width=10.0, height=4.4)
    reset_state()

    save_ht(
        ht_ext_right,
        "ch04_adj_ext_short.png",
        width=6.4,
        height=2.6,
        row_title="fooooooooooo",
        column_title="Short heatmap with automatic annotation extension",
    )
    reset_state()

    # 4.15 Manually increase space around plot via padding
    rnames = [f"R{i}" for i in range(1, 11)]
    rnames[0] = "a long long long long long row name"
    m2 = pd.DataFrame(m, index=rnames, columns=[f"C{j}" for j in range(1, 11)])
    ht_pad = Heatmap(m2, name="mat", row_names_side="left", show_row_dend=False)
    gp_l = grid_grabExpr(lambda: ch.draw(ht_pad, column_title="Default padding", newpage=False), width_in=5.0, height_in=4.5)
    gp_r = grid_grabExpr(
        lambda: ch.draw(ht_pad, padding=unit([2, 20, 2, 2], "mm"), column_title="padding = [2, 20, 2, 2] mm", newpage=False),
        width_in=5.0,
        height_in=4.5,
    )

    def _draw_pad_pair() -> None:
        pushViewport(Viewport(x=0.25, width=0.5))
        grid_draw(gp_l)
        popViewport()
        pushViewport(Viewport(x=0.75, width=0.5))
        grid_draw(gp_r)
        popViewport()

    save_Standalone(_draw_pad_pair, "ch04_padding.png", width=10.0, height=4.5)
    reset_state()


def main() -> None:
    reset_state()
    render_intro_and_4_1_to_4_5()
    render_4_6_annotations_adjusted()
    render_4_7_to_4_10()
    render_4_13_to_4_15()
    print("Chapter 4 figures generated successfully.")


if __name__ == "__main__":
    main()
