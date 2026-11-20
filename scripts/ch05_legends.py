"""Generate all figures for Chapter 5: Legends."""

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
    Legend,
    anno_lines,
    anno_points,
    anno_summary,
    colorRamp2,
    convertHeight,
    convertWidth,
    draw,
    gpar,
    grid,
    gt_render,
    packLegend,
    popViewport,
    pushViewport,
    rowAnnotation,
    unit,
    v,
    viewport,
)

from _common import reset_state, save_Standalone, save_ht

MONTH_NAMES = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]


def _save_lgd(lgd: ch.Legends, filename: str, pad_in: float = 0.45) -> None:
    w_in = max(2.2, float(convertWidth(lgd.width, "inch", valueOnly=True)) + pad_in)
    h_in = max(1.8, float(convertHeight(lgd.height, "inch", valueOnly=True)) + pad_in)
    save_Standalone(lambda: lgd.draw(), filename, width=w_in, height=h_in)


def main() -> None:
    col_fun = colorRamp2([0, 0.5, 1], ["blue", "white", "red"])

    # ---- 5.1 Continuous legends ----
    reset_state()
    lgd = Legend(col_fun=col_fun, title="foo")
    _save_lgd(lgd, "ch05_01_cont_default.png")

    # Drawing a standalone legend at explicit viewport coordinates
    reset_state()
    def _draw_positions() -> None:
        pushViewport(viewport(width=0.9, height=0.9))
        grid.rect()
        draw(lgd, x=unit(1, "cm"), y=unit(1, "cm"), just=["left", "bottom"])
        draw(lgd, x=unit(0.5, "npc"), y=unit(0.5, "npc"))
        draw(lgd, x=unit(1, "npc"), y=unit(1, "npc"), just=["right", "top"])
        popViewport()
    save_Standalone(_draw_positions, "ch05_02_draw_positions.png", width=6.5, height=4.8)

    # Preview of heatmap_legend_param & annotation_legend_param
    reset_state()
    rng = np.random.default_rng(123)
    m_preview = rng.standard_normal((10, 10))
    ht_prev = Heatmap(
        m_preview,
        name="mat",
        heatmap_legend_param={
            "title": "rnorm",
            "at": [-2, 0, 2],
            "labels": ["neg_two", "zero", "pos_two"],
        },
        top_annotation=HeatmapAnnotation(
            foo=np.arange(1, 11),
            annotation_legend_param={"foo": {"title": "foo_top_anno"}},
        ),
    )
    save_ht(ht_prev, "ch05_03_heatmap_anno_param.png", width=5.6, height=4.8)

    # Continuous vertical variants
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.25, 0.5, 0.75, 1]), "ch05_04_cont_at.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.5, 1], labels=["low", "median", "high"]), "ch05_05_cont_labels.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", legend_height=unit(6, "cm")), "ch05_06_cont_height.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", grid_width=unit(1, "cm")), "ch05_07_cont_width.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", labels_gp=gpar(col="red", font=3)), "ch05_08_cont_labels_gp.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", border="red"), "ch05_09_cont_border.png")
    reset_state()
    _save_lgd(
        Legend(col_fun=col_fun, title="foooooooo", title_position="lefttop-rot", legend_height=unit(4, "cm")),
        "ch05_10_cont_title_lefttop_rot.png",
    )
    reset_state()
    _save_lgd(
        Legend(col_fun=col_fun, title="foooooooo", title_position="leftcenter-rot", legend_height=unit(4, "cm")),
        "ch05_11_cont_title_leftcenter_rot.png",
    )

    # Mathematical formulas (Matplotlib mathtext)
    reset_state()
    _save_lgd(
        Legend(
            col_fun=col_fun,
            title=r"$\hat{\beta} = (X^T X)^{-1} X^T y$",
            at=[0, 0.25, 0.5, 0.75, 1],
            labels=[r"$\alpha$", r"$\beta$", r"$\gamma$", r"$\delta$", r"$\epsilon$"],
        ),
        "ch05_12_cont_math.png",
        pad_in=0.7,
    )

    # Rich-text title via gt_render
    reset_state()
    col_fun_gt = colorRamp2([-3, 0, 3], ["blue", "white", "red"])
    _save_lgd(
        Legend(
            col_fun=col_fun_gt,
            title=gt_render("<span style='color:orange'>**Legend title**</span>"),
            title_gp=gpar(box_fill="grey"),
            at=[-3, 0, 3],
            labels=["negative three", "zero", "positive three"],
        ),
        "ch05_13_cont_gridtext.png",
        pad_in=0.6,
    )

    # Horizontal continuous legends
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", direction="horizontal"), "ch05_14_horiz_default.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.25, 0.5, 0.75, 1], direction="horizontal"), "ch05_15_horiz_at.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.5, 1], labels=["low", "median", "high"], direction="horizontal"), "ch05_16_horiz_labels.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", legend_width=unit(6, "cm"), direction="horizontal"), "ch05_17_horiz_width.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", labels_gp=gpar(col="red", font=3), direction="horizontal"), "ch05_18_horiz_labels_gp.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", labels_rot=45, direction="horizontal"), "ch05_19_horiz_rot45.png", pad_in=0.65)
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foooooooo", direction="horizontal", title_position="topcenter"), "ch05_20_horiz_topcenter.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foooooooo", direction="horizontal", title_position="lefttop"), "ch05_21_horiz_lefttop.png")

    # Unequal interval breaks & reversed breaks
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1]), "ch05_22_unequal_vert.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.3, 1], legend_height=unit(4, "cm")), "ch05_23_unequal_no_shift.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1], direction="horizontal"), "ch05_24_unequal_horiz.png", pad_in=0.7)
    reset_state()
    _save_lgd(
        Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1], direction="horizontal", title_position="lefttop", labels_rot=90),
        "ch05_25_unequal_horiz_rot90.png",
        pad_in=0.7,
    )
    reset_state()
    _save_lgd(
        Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.5, 0.75, 1], labels=["minimal", "q10", "median", "q75", "maximal"], direction="horizontal", title_position="lefttop"),
        "ch05_26_unequal_horiz_lefttop.png",
        pad_in=0.8,
    )
    reset_state()
    _save_lgd(Legend(col_fun=col_fun, title="foo", at=[1, 0.8, 0.6, 0.4, 0.2, 0]), "ch05_27_reversed_at.png")

    # break_dist
    col_fun_prop = colorRamp2([0, 0.05, 0.1, 0.5, 1], ["green", "white", "red", "black", "blue"])
    reset_state()
    _save_lgd(Legend(col_fun=col_fun_prop, title="Prop"), "ch05_28_break_dist_default.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun_prop, title="Prop", at=[0, 0.05, 0.1, 0.5, 1]), "ch05_29_break_dist_manual_at.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun_prop, title="Prop", break_dist=1), "ch05_30_break_dist_equal.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun_prop, title="Prop", break_dist=[1, 1, 3, 3]), "ch05_31_break_dist_1133.png")
    reset_state()
    _save_lgd(Legend(col_fun=col_fun_prop, title="Prop", break_dist=[1, 1, 3, 3], legend_height=unit(4, "cm")), "ch05_32_break_dist_1133_h4.png")
    reset_state()
    col_fun2 = colorRamp2([0, 0.1, 0.1 + 1e-6, 1], ["white", "red", "yellow", "blue"])
    _save_lgd(
        Legend(col_fun=col_fun2, title="Prop", at=[0, 0.05, 0.1, 0.5, 1], break_dist=[1, 1, 3, 3], legend_height=unit(4, "cm")),
        "ch05_33_break_dist_twocolor.png",
    )

    # ---- 5.2 Discrete legends ----
    reset_state()
    _save_lgd(Legend(at=list(range(1, 7)), title="foo", legend_gp=gpar(fill=list(range(1, 7)))), "ch05_34_disc_at.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7)))), "ch05_35_disc_months.png")
    reset_state()
    at_seq = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    _save_lgd(Legend(at=at_seq, title="foo", legend_gp=gpar(fill=col_fun(at_seq).tolist())), "ch05_36_disc_from_cont.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7))), title_position="lefttop"), "ch05_37_disc_lefttop.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7))), title_position="leftcenter-rot"), "ch05_38_disc_leftcenter_rot.png")
    reset_state()
    _save_lgd(Legend(at=list(range(1, 7)), legend_gp=gpar(fill=list(range(1, 7))), title="foo", grid_height=unit(1, "cm"), grid_width=unit(5, "mm")), "ch05_39_disc_grid_size.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foo", labels_gp=gpar(col="red", fontsize=14)), "ch05_40_disc_labels_gp.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foo", title_gp=gpar(col="red", fontsize=14)), "ch05_41_disc_title_gp.png")
    reset_state()
    _save_lgd(
        Legend(
            title=gt_render("<span style='color:orange'>**Legend title**</span>"),
            title_gp=gpar(box_fill="grey"),
            at=[-3, 0, 3],
            labels=["negative three", "zero", "positive three"],
            legend_gp=gpar(fill=[1, 2, 3]),
        ),
        "ch05_42_disc_gridtext.png",
        pad_in=0.6,
    )
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foo", border="red"), "ch05_43_disc_border.png")

    # Multi-column / multi-row discrete legends
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3), "ch05_44_disc_ncol3.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, title_position="topcenter"), "ch05_45_disc_ncol3_topcenter.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, by_row=True), "ch05_46_disc_ncol3_byrow.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, gap=unit(1, "cm")), "ch05_47_disc_ncol3_gap.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, row_gap=unit(5, "mm")), "ch05_48_disc_ncol3_rowgap.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", nrow=3), "ch05_49_disc_nrow3.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foooooo", nrow=1, title_position="lefttop-rot"), "ch05_50_disc_nrow1_lefttoprot.png")
    reset_state()
    _save_lgd(Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foooooo", nrow=1, title_position="leftcenter"), "ch05_51_disc_nrow1_leftcenter.png")

    # Points, lines, boxplot, pch 26..28, multiline, custom graphics
    reset_state()
    _save_lgd(
        Legend(labels=MONTH_NAMES[:6], title="foo", type="points", pch=list(range(1, 7)), legend_gp=gpar(col=list(range(1, 7))), background="#FF8080"),
        "ch05_52_disc_points_num.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=MONTH_NAMES[:6], title="foo", type="points", pch=list("abcdef"), legend_gp=gpar(col=list(range(1, 7))), background="white"),
        "ch05_53_disc_points_char.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=MONTH_NAMES[:6], title="foo", type="lines", legend_gp=gpar(col=list(range(1, 7)), lty=list(range(1, 7))), grid_width=unit(1, "cm")),
        "ch05_54_disc_lines.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=MONTH_NAMES[:6], title="foo", type="boxplot", legend_gp=gpar(fill=list(range(1, 7)))),
        "ch05_55_disc_boxplot.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=[f"pch = {p}" for p in (26, 27, 28)], type="points", pch=[26, 27, 28]),
        "ch05_56_disc_pch26_28.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=["aaaaa\naaaaa", "bbbbb\nbbbbb", "c", "d"], legend_gp=gpar(fill=[1, 2, 3, 4])),
        "ch05_57_disc_multiline.png",
    )
    reset_state()
    _save_lgd(
        Legend(labels=["aaaaa\naaaaa", "c", "d", "bbbbb\nbbbbb"], legend_gp=gpar(fill=[1, 2, 3, 4]), nrow=2),
        "ch05_58_disc_multiline_nrow2.png",
    )
    reset_state()
    _save_lgd(
        Legend(
            labels=list("abcd"),
            graphics=[
                lambda x, y, w, h: grid.rect(x, y, w * 0.33, h, gp=gpar(fill="red")),
                lambda x, y, w, h: grid.rect(x, y, w, h * 0.33, gp=gpar(fill="blue")),
                lambda x, y, w, h: grid.text("A", x, y, gp=gpar(col="darkgreen")),
                lambda x, y, w, h: grid.points(x, y, gp=gpar(col="orange"), pch=16),
            ],
        ),
        "ch05_59_disc_custom_graphics.png",
    )

    # ---- 5.3 A list of legends (packLegend) ----
    lgd1 = Legend(at=list(range(1, 7)), legend_gp=gpar(fill=list(range(1, 7))), title="legend1")
    lgd2 = Legend(col_fun=col_fun, title="legend2", at=[0, 0.25, 0.5, 0.75, 1])
    lgd3 = Legend(labels=MONTH_NAMES[:3], legend_gp=gpar(fill=[7, 8, 9]), title="legend3")

    reset_state()
    pd_vert = packLegend(lgd1, lgd2, lgd3)
    _save_lgd(pd_vert, "ch05_60_pack_vertical.png")

    reset_state()
    pd_horiz = packLegend(lgd1, lgd2, lgd3, direction="horizontal")
    _save_lgd(pd_horiz, "ch05_61_pack_horizontal.png")

    reset_state()
    pd_mh = packLegend(lgd1, lgd3, lgd2, lgd3, lgd2, lgd1, max_height=unit(10, "cm"), column_gap=unit(1, "cm"))
    _save_lgd(pd_mh, "ch05_62_pack_max_height.png")

    reset_state()
    lgd1_row = Legend(at=list(range(1, 7)), legend_gp=gpar(fill=list(range(1, 7))), title="legend1", nrow=1)
    lgd2_row = Legend(col_fun=col_fun, title="legend2", at=[0, 0.25, 0.5, 0.75, 1], direction="horizontal")
    pd_mw = packLegend(
        lgd1_row, lgd2_row, lgd3, lgd1_row, lgd2_row, lgd3,
        max_width=unit(10, "cm"), direction="horizontal",
        column_gap=unit(5, "mm"), row_gap=unit(1, "cm"),
    )
    _save_lgd(pd_mw, "ch05_63_pack_max_width.png")

    reset_state()
    pd_pos = packLegend(lgd1_row, lgd2_row, lgd3, direction="horizontal")
    def _draw_packed_pos() -> None:
        pushViewport(viewport(width=0.8, height=0.8))
        grid.rect()
        draw(pd_pos, x=unit(1, "cm"), y=unit(1, "cm"), just=["left", "bottom"])
        draw(pd_pos, x=unit(1, "npc"), y=unit(1, "npc"), just=["right", "top"])
        popViewport()
    save_Standalone(_draw_packed_pos, "ch05_64_pack_draw_positions.png", width=7.5, height=5.0)

    # ---- 5.4 Heatmap and annotation legends ----
    rng = np.random.default_rng(123)
    m = rng.standard_normal((10, 10))

    reset_state()
    ht_hlp = Heatmap(
        m, name="mat",
        heatmap_legend_param={
            "at": [-2, 0, 2],
            "labels": ["low", "zero", "high"],
            "title": "Some values",
            "legend_height": unit(4, "cm"),
            "title_position": "lefttop-rot",
        },
    )
    save_ht(ht_hlp, "ch05_65_hm_legend_param.png", width=5.6, height=4.8)

    reset_state()
    ha_alp = HeatmapAnnotation(
        foo=rng.uniform(size=10),
        bar=rng.choice(["f", "m"], 10, replace=True),
        annotation_legend_param={
            "foo": {"title": "Fooooooh", "at": [0, 0.5, 1], "labels": ["zero", "median", "one"]},
            "bar": {"title": "Baaaaaaar", "at": ["f", "m"], "labels": ["Female", "Male"]},
        },
    )
    save_ht(Heatmap(m, name="mat", top_annotation=ha_alp), "ch05_66_anno_legend_param.png", width=5.6, height=4.8)

    reset_state()
    ha1 = HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True))
    ha2 = HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
    ht_hgrp = Heatmap(m, name="mat1", top_annotation=ha1) + rowAnnotation(sth=rng.uniform(size=10)) + Heatmap(m, name="mat2", top_annotation=ha2)
    save_ht(ht_hgrp, "ch05_67_horiz_grouping.png", width=8.8, height=5.4, padding=unit([10, 2, 2, 2], "mm"))

    reset_state()
    ha1_v = HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True), annotation_name_side="left")
    ha2_v = HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
    ht_vgrp = Heatmap(m, name="mat1", top_annotation=ha1_v) % v % Heatmap(m, name="mat2", top_annotation=ha2_v, right_annotation=rowAnnotation(sth=np.arange(1, 11)))
    save_ht(ht_vgrp, "ch05_68_vert_grouping.png", width=6.0, height=8.0, padding=unit([10, 2, 2, 2], "mm"))

    reset_state()
    ha_hide = HeatmapAnnotation(
        foo=rng.uniform(size=10),
        bar=rng.choice(["f", "m"], 10, replace=True),
        show_legend=[True, False],
        annotation_name_side="left",
    )
    ht_hide = Heatmap(m, name="mat1", top_annotation=ha_hide) + Heatmap(m, name="mat2", show_heatmap_legend=False)
    save_ht(ht_hide, "ch05_69_hide_legends.png", width=8.0, height=5.0)

    reset_state()
    ha_m1 = HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True))
    ha_m2 = rowAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(list("abc"), 10, replace=True))
    ha_m3 = rowAnnotation(foo3=rng.uniform(size=10), bar3=rng.choice(MONTH_NAMES[:3], 10, replace=True))
    ht_list_m = Heatmap(m, name="mat1", top_annotation=ha_m1) + Heatmap(m, name="mat2", left_annotation=ha_m2) + ha_m3
    save_ht(ht_list_m, "ch05_70_merge_legends.png", width=8.5, height=5.4, merge_legends=True)

    reset_state()
    save_ht(ht_list_m, "ch05_71_legend_grouping_original.png", width=8.5, height=5.4, legend_grouping="original")

    reset_state()
    ht_cbd = Heatmap(
        m, name="mat",
        heatmap_legend_param={"color_bar": "discrete"},
        top_annotation=HeatmapAnnotation(
            foo=np.arange(1, 11),
            annotation_legend_param={"foo": {"color_bar": "discrete"}},
        ),
    )
    save_ht(ht_cbd, "ch05_72_color_bar_discrete.png", width=5.6, height=4.8)

    reset_state()
    chr_vec = rng.choice(list("abc"), 10, replace=True)
    ht_ord = Heatmap(
        m,
        name="mat",
        top_annotation=HeatmapAnnotation(
            chr=chr_vec,
            fa_ordered=chr_vec,
            col={"fa_ordered": {"c": "#1B9E77", "a": "#D95F02", "b": "#7570B3"}},
            annotation_legend_param={"fa_ordered": {"at": ["b", "c", "a"]}},
        ),
    )
    save_ht(ht_ord, "ch05_73_level_order.png", width=5.8, height=4.8)

    # ---- 5.5 Add customized legends ----
    reset_state()
    ha_pt = HeatmapAnnotation(pt=anno_points(np.arange(1, 11), gp=gpar(col=np.repeat([2, 3], 5)), height=unit(2, "cm")), show_annotation_name=False)
    ha_ln = HeatmapAnnotation(ln=anno_lines(np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)]), gp=gpar(col=[4, 5], lty=[1, 2]), height=unit(2, "cm")), show_annotation_name=False)
    ht_cust = (
        Heatmap(m, name="mat1", top_annotation=ha_pt)
        + Heatmap(m, name="mat2", top_annotation=ha_ln)
        + Heatmap(m[:, [0]], name="mat3", top_annotation=HeatmapAnnotation(summary=anno_summary(gp=gpar(fill=[2, 3]))), width=unit(1, "cm"))
    )
    save_ht(ht_cust, "ch05_74_custom_annos_no_legend.png", width=7.5, height=5.4, ht_gap=unit(7, "mm"), row_km=2)

    reset_state()
    lgd_list = [
        Legend(labels=["red", "green"], title="pt", type="points", pch=16, legend_gp=gpar(col=[2, 3])),
        Legend(labels=["darkblue", "lightblue"], title="ln", type="lines", legend_gp=gpar(col=[4, 5], lty=[1, 2])),
        Legend(labels=["group1", "group2"], title="km", type="boxplot", legend_gp=gpar(fill=[2, 3])),
    ]
    save_ht(
        ht_cust, "ch05_75_custom_annos_with_legend.png",
        width=8.2, height=5.4, ht_gap=unit(7, "mm"), row_km=2, annotation_legend_list=lgd_list,
    )

    # ---- 5.6 The side of legends ----
    reset_state()
    ha_s1 = HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True), annotation_name_side="left")
    ha_s2 = HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
    ht_sides = Heatmap(m, name="mat1", top_annotation=ha_s1) + rowAnnotation(sth=rng.uniform(size=10), annotation_name_side="top") + Heatmap(m, name="mat2", top_annotation=ha_s2)
    save_ht(
        ht_sides,
        "ch05_76_legend_sides.png",
        width=9.2,
        height=6.0,
        heatmap_legend_side="left",
        annotation_legend_side="bottom",
        ht_gap=unit(4, "mm"),
        padding=unit([2, 2, 2, 12], "mm"),
    )

    reset_state()
    ha_h1 = HeatmapAnnotation(
        foo1=rng.uniform(size=10),
        bar1=rng.choice(["f", "m"], 10, replace=True),
        annotation_name_side="left",
        annotation_legend_param={"foo1": {"legend_direction": "horizontal"}, "bar1": {"nrow": 1}},
    )
    ha_h2 = HeatmapAnnotation(
        foo2=rng.uniform(size=10),
        bar2=rng.choice(["f", "m"], 10, replace=True),
        annotation_legend_param={"foo2": {"legend_direction": "horizontal"}, "bar2": {"nrow": 1}},
    )
    ht_all_horiz = (
        Heatmap(m, name="mat1", top_annotation=ha_h1, heatmap_legend_param={"legend_direction": "horizontal"})
        + rowAnnotation(sth=rng.uniform(size=10), annotation_name_side="top", annotation_legend_param={"sth": {"legend_direction": "horizontal"}})
        + Heatmap(m, name="mat2", top_annotation=ha_h2, heatmap_legend_param={"legend_direction": "horizontal"})
    )
    save_ht(
        ht_all_horiz,
        "ch05_77_all_horizontal_bottom.png",
        width=9.2,
        height=6.4,
        merge_legends=True,
        heatmap_legend_side="bottom",
        annotation_legend_side="bottom",
        ht_gap=unit(4, "mm"),
        padding=unit([2, 10, 2, 12], "mm"),
    )

    print("Generated 77 figures for Chapter 5.")


if __name__ == "__main__":
    main()
