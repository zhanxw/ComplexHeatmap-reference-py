"""Generate all Chapter 3 figures: Heatmap Annotations."""

from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy.stats import gaussian_kde

import complexheatmap as ch
from complexheatmap import (
    AnnotationFunction,
    Heatmap,
    HeatmapAnnotation,
    Legend,
    Viewport,
    anno_barplot,
    anno_block,
    anno_boxplot,
    anno_customize,
    anno_density,
    anno_empty,
    anno_histogram,
    anno_horizon,
    anno_image,
    anno_joyplot,
    anno_lines,
    anno_link,
    anno_mark,
    anno_numeric,
    anno_points,
    anno_simple,
    anno_summary,
    anno_text,
    anno_textbox,
    c,
    colorRamp2,
    column_order,
    decorate_annotation,
    gpar,
    grid_boxplot,
    grid_circle,
    grid_draw,
    grid_grabExpr,
    grid_lines,
    grid_points,
    grid_raster,
    grid_rect,
    grid_roundrect,
    grid_segments,
    grid_text,
    gt_render,
    max_text_width,
    popViewport,
    pushViewport,
    rand_color,
    rowAnnotation,
    seekViewport,
    textbox_grob,
    unit,
)

from _common import reset_state, save_Standalone, save_decorated, save_ht

MONTH_NAMES = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]


def save_anno(
    ha: HeatmapAnnotation,
    filename: str,
    index: list[int] | None = None,
    width: float | None = None,
    height: float | None = None,
    dpi: int = 260,
    margin_cm: float = 3.6,
) -> Path:
    """Render a standalone HeatmapAnnotation at proportional dimensions."""
    ext = ha.extended.to_mm("x")
    if ha.which == "column":
        h_mm = ch.convertHeight(ha.height, "mm", valueOnly=True)
        fig_w = width if width is not None else 7.2
        fig_h = height if height is not None else max(1.35, (h_mm + ext[0] + ext[2] + 24.0) / 25.4)
    else:
        w_mm = ch.convertWidth(ha.width, "mm", valueOnly=True)
        fig_w = width if width is not None else max(3.2, (w_mm + ext[1] + ext[3] + 36.0) / 25.4)
        fig_h = height if height is not None else 4.8

    def _draw() -> None:
        if ha.which == "column":
            pushViewport(
                Viewport(
                    width=unit(1, "npc") - unit(margin_cm, "cm") - unit(float(ext[1] + ext[3]), "mm"),
                    height=ha.height,
                )
            )
        else:
            pushViewport(
                Viewport(
                    height=unit(1, "npc") - unit(2.2, "cm") - unit(float(ext[0] + ext[2]), "mm"),
                    width=ha.width,
                )
            )
        try:
            ha.draw(index=index, test=False)
        finally:
            popViewport()

    out = save_Standalone(_draw, filename, width=fig_w, height=fig_h, dpi=dpi)
    reset_state()
    return out


def _generate_demo_icons(out_dir: Path) -> list[str]:
    """Create 10 crisp offline RGBA PNG icons for anno_image and anno_customize."""
    palette = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3",
    ]
    paths: list[str] = []
    size = 128
    pad = 16
    c_mid = size // 2
    for idx, color in enumerate(palette):
        im = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        draw = ImageDraw.Draw(im)
        kind = idx % 10
        if kind == 0:  # solid circle
            draw.ellipse([pad, pad, size - pad, size - pad], fill=color)
        elif kind == 1:  # ring
            draw.ellipse([pad, pad, size - pad, size - pad], outline=color, width=14)
        elif kind == 2:  # rounded square
            draw.rounded_rectangle([pad, pad, size - pad, size - pad], radius=22, fill=color)
        elif kind == 3:  # diamond
            draw.polygon([(c_mid, pad), (size - pad, c_mid), (c_mid, size - pad), (pad, c_mid)], fill=color)
        elif kind == 4:  # triangle
            draw.polygon([(c_mid, pad), (size - pad, size - pad), (pad, size - pad)], fill=color)
        elif kind == 5:  # 5-point star
            pts = []
            for k in range(10):
                ang = -math.pi / 2 + k * math.pi / 5
                r = (size // 2 - pad) if k % 2 == 0 else (size // 2 - pad) * 0.42
                pts.append((c_mid + r * math.cos(ang), c_mid + r * math.sin(ang)))
            draw.polygon(pts, fill=color)
        elif kind == 6:  # plus cross
            t = 18
            draw.rectangle([c_mid - t, pad, c_mid + t, size - pad], fill=color)
            draw.rectangle([pad, c_mid - t, size - pad, c_mid + t], fill=color)
        elif kind == 7:  # hexagon
            pts = [
                (c_mid + 46 * math.cos(k * math.pi / 3), c_mid + 46 * math.sin(k * math.pi / 3))
                for k in range(6)
            ]
            draw.polygon(pts, fill=color)
        elif kind == 8:  # mini bar chart icon
            draw.rectangle([20, 68, 42, 110], fill=color)
            draw.rectangle([52, 34, 74, 110], fill=color)
            draw.rectangle([84, 50, 106, 110], fill=color)
        else:  # checkmark inside circle
            draw.ellipse([pad, pad, size - pad, size - pad], outline=color, width=10)
            draw.line([(36, 68), (56, 88), (94, 42)], fill=color, width=12)
        p = out_dir / f"icon_{idx}.png"
        im.save(p)
        paths.append(str(p))
    return paths


def render_intro() -> None:
    rng = np.random.default_rng(123)
    mat = pd.DataFrame(
        rng.standard_normal((10, 10)),
        index=[f"R{i}" for i in range(1, 11)],
        columns=[f"C{j}" for j in range(1, 11)],
    )
    column_ha = HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=anno_barplot(rng.uniform(size=10)))
    row_ha = rowAnnotation(foo2=rng.uniform(size=10), bar2=anno_barplot(rng.uniform(size=10)))

    ht1 = Heatmap(mat, name="mat", top_annotation=column_ha, right_annotation=row_ha)
    save_ht(ht1, "ch03_intro_top_right.png", width=6.4, height=5.2)
    reset_state()

    ht2 = Heatmap(mat, name="mat", bottom_annotation=column_ha, left_annotation=row_ha)
    save_ht(ht2, "ch03_intro_bottom_left.png", width=6.4, height=5.2)
    reset_state()


def render_3_1_simple() -> None:
    rng = np.random.default_rng(123)

    save_anno(HeatmapAnnotation(foo=np.arange(1, 11)), "ch03_simple_numeric.png")
    bars = rng.choice(["a", "b", "c"], 10)
    save_anno(HeatmapAnnotation(bar=bars), "ch03_simple_discrete.png")

    col_fun = colorRamp2([0, 5, 10], ["blue", "white", "red"])
    save_anno(HeatmapAnnotation(foo=np.arange(1, 11), col={"foo": col_fun}), "ch03_simple_colfun.png")

    disc_col = {"a": "red", "b": "green", "c": "blue"}
    save_anno(HeatmapAnnotation(bar=bars, col={"bar": disc_col}), "ch03_simple_discrete_col.png")

    save_anno(
        HeatmapAnnotation(foo=np.arange(1, 11), bar=bars, col={"foo": col_fun, "bar": disc_col}),
        "ch03_simple_multi.png",
    )

    foo_na = np.array([1, 2, 3, 4, np.nan, 6, 7, 8, 9, 10], dtype=float)
    bar_na = [None] + list(rng.choice(["a", "b", "c"], 9))
    save_anno(
        HeatmapAnnotation(foo=foo_na, bar=bar_na, col={"foo": col_fun, "bar": disc_col}, na_col="black"),
        "ch03_simple_na.png",
    )

    save_anno(
        HeatmapAnnotation(foo=np.arange(1, 11), bar=bars, col={"foo": col_fun, "bar": disc_col}, gp=gpar(col="black")),
        "ch03_simple_gp.png",
    )

    df_mat = pd.DataFrame({"a": rng.uniform(size=10), "b": rng.uniform(size=10)})
    save_anno(HeatmapAnnotation(foo=df_mat), "ch03_simple_matrix_named.png")
    save_anno(
        HeatmapAnnotation(foo=np.column_stack([rng.uniform(size=10), rng.uniform(size=10)])),
        "ch03_simple_matrix_unnamed.png",
    )

    anno_df = pd.DataFrame({"foo": np.arange(1, 11), "bar": bars})
    save_anno(
        HeatmapAnnotation(df=anno_df, col={"foo": col_fun, "bar": disc_col}),
        "ch03_simple_df.png",
    )
    save_anno(
        HeatmapAnnotation(df=anno_df, foo2=rng.standard_normal(10), col={"foo": col_fun, "bar": disc_col}),
        "ch03_simple_df_mixed.png",
    )

    save_anno(
        HeatmapAnnotation(
            foo=np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)]),
            bar=bars,
            col={"foo": col_fun, "bar": disc_col},
            border=True,
        ),
        "ch03_simple_border.png",
    )

    save_anno(
        HeatmapAnnotation(
            foo=pd.DataFrame({"a": np.arange(1, 11), "b": np.arange(10, 0, -1)}),
            bar=bars,
            col={"foo": col_fun, "bar": disc_col},
            simple_anno_size=unit(1, "cm"),
        ),
        "ch03_simple_size.png",
    )


def render_3_2_anno_simple() -> None:
    rng = np.random.default_rng(123)

    save_anno(
        HeatmapAnnotation(
            foo=anno_simple(np.arange(1, 11), pch=[1] * 10, pt_gp=gpar(col="red"), pt_size=unit(np.arange(1, 11), "mm"))
        ),
        "ch03_anno_simple_pt_size.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_simple(np.arange(1, 11), pch=list(range(1, 11)))),
        "ch03_anno_simple_pch_vec.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_simple(np.arange(1, 11), pch=rng.choice(["a", "b", "c"], 10).tolist())),
        "ch03_anno_simple_pch_letters.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_simple(np.arange(1, 11), pch=[1, 2, 3, 4, None, 6, 7, 8, None, 10])),
        "ch03_anno_simple_pch_na.png",
    )

    mat2 = np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)])
    save_anno(
        HeatmapAnnotation(foo=anno_simple(mat2, pch=np.tile([1, 2], (10, 1)))),
        "ch03_anno_simple_mat_pch_col.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_simple(mat2, pch=np.arange(1, 11))),
        "ch03_anno_simple_mat_pch_row.png",
    )
    pch_mat = np.arange(1, 21, dtype=object).reshape(10, 2)
    pch_mat.ravel()[rng.choice(20, 10, replace=False)] = None
    save_anno(
        HeatmapAnnotation(foo=anno_simple(mat2, pch=pch_mat)),
        "ch03_anno_simple_mat_pch_mat.png",
    )

    pvalue = 10.0 ** (-rng.uniform(0, 3, 10))
    is_sig = pvalue < 0.01
    pch_sig = ["*" if s else None for s in is_sig]
    pvalue_col_fun = colorRamp2([0, 2, 3], ["green", "white", "red"])
    ha_pv = HeatmapAnnotation(
        pvalue=anno_simple(-np.log10(pvalue), col=pvalue_col_fun, pch=pch_sig),
        annotation_name_side="left",
    )
    ht_pv = Heatmap(rng.standard_normal((10, 10)), name="mat", top_annotation=ha_pv)
    lgd_pvalue = Legend(title="p-value", col_fun=pvalue_col_fun, at=[0, 1, 2, 3], labels=["1", "0.1", "0.01", "0.001"])
    lgd_sig = Legend(pch="*", type="points", labels=["< 0.01"])
    save_ht(
        ht_pv,
        "ch03_anno_simple_pvalue.png",
        width=6.6,
        height=5.0,
        annotation_legend_list=[lgd_pvalue, lgd_sig],
        padding=unit([2, 10, 2, 2], "mm"),
    )
    reset_state()

    save_anno(
        HeatmapAnnotation(foo=anno_simple(np.arange(1, 11), height=unit(2, "cm"))),
        "ch03_anno_simple_height.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_simple(mat2, simple_anno_size=unit(2, "cm"))),
        "ch03_anno_simple_sas.png",
    )


def render_3_3_empty() -> None:
    rng = np.random.default_rng(123)

    save_anno(
        HeatmapAnnotation(foo=anno_empty(border=True)),
        "ch03_anno_empty_basic.png",
        index=list(range(10)),
    )

    def random_text(n: int) -> list[str]:
        letters = list("abcdefghijklmnopqrstuvwxyz")
        return ["".join(rng.choice(letters, int(rng.integers(4, 11)))) for _ in range(n)]

    text_list = [random_text(4) for _ in range(4)]
    flat_txt = [w for grp in text_list for w in grp]
    ha = rowAnnotation(foo=anno_empty(border=False, width=max_text_width(flat_txt) + unit(4, "mm")))
    ht = Heatmap(rng.standard_normal((100, 10)), name="mat", row_km=4, right_annotation=ha)

    def _dec_slices() -> None:
        for i in range(1, 5):
            with decorate_annotation("foo", slice=i):
                grid_rect(x=0, width=unit(2, "mm"), gp=gpar(fill=i, col=None), just="left")
                grid_text("\n".join(text_list[i - 1]), x=unit(4, "mm"), just="left")

    save_decorated(ht, "ch03_anno_empty_slices.png", _dec_slices, width=6.5, height=5.8)
    reset_state()

    ha_pts = HeatmapAnnotation(foo=anno_empty(border=True, height=unit(3, "cm")))
    ht_pts = Heatmap(rng.standard_normal((10, 10)), name="mat", top_annotation=ha_pts)
    value = rng.uniform(size=10)

    def _dec_pts(res) -> None:
        co = column_order(res)
        with decorate_annotation("foo"):
            x_pos = np.arange(1, 11)
            val_ord = value[co]
            pushViewport(Viewport(xscale=(0.5, 10.5), yscale=(0, 1)))
            grid_lines([0.5, 10.5], [0.5, 0.5], gp=gpar(lty=2), default_units="native")
            grid_points(
                x_pos,
                val_ord,
                pch=16,
                size=unit(2, "mm"),
                gp=gpar(col=["red" if v > 0.5 else "blue" for v in val_ord]),
                default_units="native",
            )
            ch.grid_yaxis(at=[0, 0.5, 1])
            popViewport()

    save_decorated(ht_pts, "ch03_anno_empty_points.png", _dec_pts, width=6.4, height=5.4)
    reset_state()


def render_3_4_block() -> None:
    rng = np.random.default_rng(123)
    m10 = rng.standard_normal((10, 10))

    save_ht(
        Heatmap(m10, name="mat", top_annotation=HeatmapAnnotation(foo=anno_block(gp=gpar(fill=[2, 3, 4]))), column_km=3),
        "ch03_anno_block_basic.png",
        width=6.0,
        height=4.8,
    )
    reset_state()

    save_ht(
        Heatmap(
            m10,
            name="mat",
            top_annotation=HeatmapAnnotation(
                foo=anno_block(gp=gpar(fill=[2, 3, 4]), labels=["group1", "group2", "group3"], labels_gp=gpar(col="white", fontsize=10))
            ),
            column_km=3,
            left_annotation=rowAnnotation(
                foo=anno_block(gp=gpar(fill=[2, 3, 4]), labels=["group1", "group2", "group3"], labels_gp=gpar(col="white", fontsize=10))
            ),
            row_km=3,
        ),
        "ch03_anno_block_labels.png",
        width=6.2,
        height=5.0,
    )
    reset_state()

    mat2 = rng.standard_normal((50, 50))
    split = np.repeat(np.arange(1, 6), 10)
    ha_step1 = HeatmapAnnotation(
        empty=anno_empty(border=False),
        foo=anno_block(gp=gpar(fill=list(range(2, 7))), labels=list("ABCDE")),
    )
    save_ht(
        Heatmap(mat2, name="mat2", column_split=split, top_annotation=ha_step1, column_title=None),
        "ch03_anno_block_multi_step1.png",
        width=6.4,
        height=5.2,
    )
    reset_state()

    def group_block_anno(group: list[int], empty_anno: str, gp_box=None, label: str | None = None, label_gp=None) -> None:
        vp1 = seekViewport(f"annotation_{empty_anno}_{min(group)}")
        vp2 = seekViewport(f"annotation_{empty_anno}_{max(group)}")
        vp_glob = seekViewport("global")
        x0 = unit(vp1.x0_mm - vp_glob.x0_mm, "mm")
        y0 = unit(vp1.y0_mm - vp_glob.y0_mm, "mm")
        w = unit((vp2.x0_mm + vp2.w_mm) - vp1.x0_mm, "mm")
        h = unit((vp2.y0_mm + vp2.h_mm) - vp1.y0_mm, "mm")
        grid_rect(x0, y0, width=w, height=h, just=["left", "bottom"], gp=gp_box or gpar())
        if label is not None:
            grid_text(label, x=x0 + w * 0.5, y=y0 + h * 0.5, gp=label_gp or gpar())

    ha_step2 = HeatmapAnnotation(
        empty=anno_empty(border=False, height=unit(8, "mm")),
        foo=anno_block(gp=gpar(fill=list(range(2, 7))), labels=list("ABCDE")),
    )
    ht_multi = Heatmap(mat2, name="mat2", column_split=split, top_annotation=ha_step2, column_title=None)
    save_decorated(
        ht_multi,
        "ch03_anno_block_multi_one.png",
        lambda: group_block_anno([1, 2, 3], "empty", gp_box=gpar(fill="red"), label="group 1", label_gp=gpar(col="white")),
        width=6.4,
        height=5.2,
    )
    reset_state()

    save_decorated(
        ht_multi,
        "ch03_anno_block_multi_two.png",
        lambda: (
            group_block_anno([1, 2, 3], "empty", gp_box=gpar(fill="red"), label="group 1", label_gp=gpar(col="white")),
            group_block_anno([4, 5], "empty", gp_box=gpar(fill="blue"), label="group 2", label_gp=gpar(col="white")),
        ),
        width=6.4,
        height=5.2,
    )
    reset_state()

    col_map = {"1": "red", "2": "blue", "A": "green", "B": "orange"}

    def _pf(index, levels):
        lv = [str(levels)] if isinstance(levels, str) else [str(v) for v in levels]
        grid_rect(gp=gpar(fill=col_map.get(lv[-1], "grey"), col="black"))
        txt = ",".join(lv) + f"\n{len(index)} rows"
        grid_text(txt, 0.5, 0.5, rot=0, gp=gpar(col=col_map.get(lv[0], "black")))

    save_ht(
        Heatmap(m10, name="mat", row_km=2, row_split=rng.choice(["A", "B"], 10))
        + rowAnnotation(foo=anno_block(panel_fun=_pf, width=unit(3, "cm"))),
        "ch03_anno_block_panel_fun.png",
        width=7.2,
        height=5.2,
    )
    reset_state()

    sp_ab = rng.choice(["A", "B"], 10)
    align_a = {"A": [i for i, v in enumerate(sp_ab) if v == "A"]}

    def _pf_count(index, nm):
        grid_rect()
        grid_text(f"{len(index)} rows", 0.5, 0.5)

    save_ht(
        Heatmap(m10, name="mat", row_split=sp_ab)
        + rowAnnotation(foo=anno_block(align_to=align_a, panel_fun=_pf_count, width=unit(3, "cm"))),
        "ch03_anno_block_align_to_split.png",
        width=6.4,
        height=5.0,
    )
    reset_state()

    save_ht(
        Heatmap(m10, name="mat", cluster_rows=False)
        + rowAnnotation(foo=anno_block(align_to={"A": [1, 2, 3], "B": [6, 7, 8]}, panel_fun=_pf_count, width=unit(3, "cm"))),
        "ch03_anno_block_align_to_nosplit.png",
        width=6.4,
        height=5.0,
    )
    reset_state()


def render_3_5_image(icon_paths: list[str]) -> None:
    save_anno(HeatmapAnnotation(foo=anno_image(icon_paths)), "ch03_anno_image_basic.png")
    save_anno(
        HeatmapAnnotation(foo=anno_image(icon_paths, gp=gpar(fill=list(range(1, 11)), col="black"))),
        "ch03_anno_image_gp.png",
    )
    save_anno(HeatmapAnnotation(foo=anno_image(icon_paths, border="red")), "ch03_anno_image_border.png")
    save_anno(HeatmapAnnotation(foo=anno_image(icon_paths, space=unit(3, "mm"))), "ch03_anno_image_space.png")
    partial = ["", ""] + list(icon_paths[2:])
    save_anno(HeatmapAnnotation(foo=anno_image(partial)), "ch03_anno_image_partial.png")


def render_3_6_points() -> None:
    rng = np.random.default_rng(123)
    save_anno(HeatmapAnnotation(foo=anno_points(rng.uniform(size=10))), "ch03_anno_points_vec.png")
    save_anno(
        HeatmapAnnotation(foo=anno_points(rng.uniform(size=(10, 2)), pch=[1, 2], gp=gpar(col=[2, 3]))),
        "ch03_anno_points_mat.png",
    )
    save_anno(
        HeatmapAnnotation(
            foo=anno_points(
                rng.uniform(size=10),
                ylim=(0, 1),
                axis_param={"side": "right", "at": [0, 0.5, 1], "labels": ["zero", "half", "one"]},
            )
        ),
        "ch03_anno_points_axis.png",
    )
    save_anno(
        rowAnnotation(
            foo=anno_points(
                rng.uniform(size=10),
                which="row",
                ylim=(0, 1),
                width=unit(2, "cm"),
                axis_param={"side": "bottom", "at": [0, 0.5, 1], "labels": ["zero", "half", "one"], "labels_rot": 45},
            )
        ),
        "ch03_anno_points_row_rot.png",
        width=3.6,
        height=4.8,
    )
    save_anno(
        HeatmapAnnotation(foo=anno_points(rng.uniform(size=10), height=unit(2, "cm"))),
        "ch03_anno_points_height.png",
    )


def render_3_7_lines() -> None:
    rng = np.random.default_rng(123)
    save_anno(HeatmapAnnotation(foo=anno_lines(rng.uniform(size=10))), "ch03_anno_lines_vec.png")
    mat_l = np.column_stack([np.tile(np.arange(1, 6), 2), np.tile(np.arange(5, 0, -1), 2)])
    save_anno(
        HeatmapAnnotation(foo=anno_lines(mat_l, gp=gpar(col=[2, 3]), add_points=True, pt_gp=gpar(col=[5, 6]), pch=[1, 16])),
        "ch03_anno_lines_mat.png",
    )
    save_anno(HeatmapAnnotation(foo=anno_lines(rng.uniform(size=10), smooth=True)), "ch03_anno_lines_smooth.png")


def render_3_8_barplot() -> None:
    rng = np.random.default_rng(123)
    save_anno(HeatmapAnnotation(foo=anno_barplot(np.arange(1, 11))), "ch03_anno_barplot_vec.png")
    save_anno(HeatmapAnnotation(foo=anno_barplot(np.arange(1, 11), bar_width=1)), "ch03_anno_barplot_width1.png")
    save_anno(
        HeatmapAnnotation(foo=anno_barplot(np.arange(1, 11), gp=gpar(fill=list(range(1, 11))))),
        "ch03_anno_barplot_gp.png",
    )
    save_anno(HeatmapAnnotation(foo=anno_barplot(np.arange(-5, 6), baseline="min")), "ch03_anno_barplot_baseline_min.png")
    save_anno(HeatmapAnnotation(foo=anno_barplot(np.arange(-5, 6), baseline=0)), "ch03_anno_barplot_baseline_0.png")

    mat_b = np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)])
    save_anno(HeatmapAnnotation(foo=anno_barplot(mat_b)), "ch03_anno_barplot_stacked.png")
    save_anno(
        HeatmapAnnotation(foo=anno_barplot(mat_b, gp=gpar(fill=[2, 3], col=[2, 3]))),
        "ch03_anno_barplot_stacked_gp.png",
    )

    m_prop = rng.uniform(size=(10, 4))
    m_prop = m_prop / m_prop.sum(axis=1, keepdims=True)
    save_anno(
        HeatmapAnnotation(foo=anno_barplot(m_prop, gp=gpar(fill=[2, 3, 4, 5]), bar_width=1, height=unit(6, "cm"))),
        "ch03_anno_barplot_prop.png",
    )

    ha_rev = (
        rowAnnotation(
            axis_reverse=anno_barplot(
                m_prop, gp=gpar(fill=[2, 3, 4, 5]), axis_param={"direction": "reverse"}, bar_width=1, width=unit(4, "cm")
            )
        )
        + rowAnnotation(axis_normal=anno_barplot(m_prop, gp=gpar(fill=[2, 3, 4, 5]), bar_width=1, width=unit(4, "cm")))
    )
    save_ht(ha_rev, "ch03_anno_barplot_reverse.png", width=5.6, height=5.8, ht_gap=unit(4, "mm"))
    reset_state()

    save_anno(
        HeatmapAnnotation(foo=anno_barplot(np.arange(1, 11), add_numbers=True, height=unit(1, "cm"))),
        "ch03_anno_barplot_numbers.png",
    )
    save_anno(HeatmapAnnotation(foo=anno_barplot(mat_b, beside=True)), "ch03_anno_barplot_beside.png")
    save_anno(HeatmapAnnotation(foo=anno_barplot(mat_b, beside=True, attach=True)), "ch03_anno_barplot_attach.png")


def render_3_9_boxplot() -> None:
    rng = np.random.default_rng(12345)
    m = rng.standard_normal((10, 10))
    save_anno(HeatmapAnnotation(foo=anno_boxplot(m, height=unit(4, "cm"))), "ch03_anno_boxplot_basic.png")
    save_anno(
        HeatmapAnnotation(foo=anno_boxplot(m, height=unit(4, "cm"), gp=gpar(fill=list(range(1, 11))))),
        "ch03_anno_boxplot_gp.png",
    )
    save_anno(
        HeatmapAnnotation(foo=anno_boxplot(m, height=unit(4, "cm"), box_width=0.9, outline=False)),
        "ch03_anno_boxplot_no_outline.png",
    )


def render_3_10_histogram() -> None:
    rng = np.random.default_rng(123)
    m = rng.standard_normal((10, 100))
    save_anno(rowAnnotation(foo=anno_histogram(m)), "ch03_anno_histogram_basic.png")
    save_anno(rowAnnotation(foo=anno_histogram(m, n_breaks=20)), "ch03_anno_histogram_breaks.png")
    save_anno(
        rowAnnotation(foo=anno_histogram(m, gp=gpar(fill=list(range(1, 11))))),
        "ch03_anno_histogram_gp.png",
    )


def render_3_11_density() -> None:
    rng = np.random.default_rng(123)
    m = rng.standard_normal((10, 100))
    save_anno(rowAnnotation(foo=anno_density(m)), "ch03_anno_density_lines.png")
    save_anno(
        rowAnnotation(foo=anno_density(m, joyplot_scale=2, gp=gpar(fill="#CCCCCC80"))),
        "ch03_anno_density_joy.png",
    )
    save_anno(
        rowAnnotation(foo=anno_density(m, type="violin", gp=gpar(fill=list(range(1, 11))))),
        "ch03_anno_density_violin.png",
    )
    m2 = rng.standard_normal((50, 10))
    save_anno(
        rowAnnotation(foo=anno_density(m2, type="heatmap", width=unit(6, "cm"))),
        "ch03_anno_density_heatmap.png",
    )
    save_anno(
        rowAnnotation(foo=anno_density(m2, type="heatmap", width=unit(6, "cm"), heatmap_colors=["white", "orange"])),
        "ch03_anno_density_heatmap_orange.png",
    )


def render_3_12_joyplot() -> None:
    rng = np.random.default_rng(123)
    m1 = rng.standard_normal((100, 10))
    gx = np.linspace(-4, 4, 200)
    lt1 = [pd.DataFrame({"x": gx, "y": gaussian_kde(m1[:, j])(gx)}) for j in range(m1.shape[1])]
    save_anno(
        rowAnnotation(foo=anno_joyplot(lt1, width=unit(4, "cm"), gp=gpar(fill=list(range(1, 11))), transparency=0.75)),
        "ch03_anno_joyplot_fill.png",
    )

    m2 = rng.standard_normal((100, 50))
    lt2 = [pd.DataFrame({"x": gx, "y": gaussian_kde(m2[:, j])(gx)}) for j in range(m2.shape[1])]
    save_anno(
        rowAnnotation(foo=anno_joyplot(lt2, width=unit(4, "cm"), gp=gpar(fill=None), scale=4)),
        "ch03_anno_joyplot_lines.png",
    )


def render_3_13_horizon() -> None:
    rng = np.random.default_rng(123)
    lt = [np.cumprod(1 + rng.uniform(-x / 100.0, x / 100.0, 1000)) - 1.0 for x in range(1, 21)]
    save_anno(rowAnnotation(foo=anno_horizon(lt)), "ch03_anno_horizon_basic.png")
    save_anno(
        rowAnnotation(foo=anno_horizon(lt, gp=gpar(pos_fill="orange", neg_fill="darkgreen"))),
        "ch03_anno_horizon_col.png",
    )
    save_anno(
        rowAnnotation(
            foo=anno_horizon(
                lt,
                gp=gpar(pos_fill=np.repeat(["orange", "red"], 10), neg_fill=np.repeat(["darkgreen", "blue"], 10)),
            )
        ),
        "ch03_anno_horizon_vec_col.png",
    )
    save_anno(rowAnnotation(foo=anno_horizon(lt, negative_from_top=True)), "ch03_anno_horizon_neg_top.png")
    save_anno(rowAnnotation(foo=anno_horizon(lt, gap=unit(1, "mm"))), "ch03_anno_horizon_gap.png")


def render_3_14_text() -> None:
    rng = np.random.default_rng(123)
    save_anno(
        rowAnnotation(foo=anno_text(MONTH_NAMES, gp=gpar(fontsize=[i + 5 for i in range(12)]))),
        "ch03_anno_text_fontsize.png",
    )
    save_anno(
        rowAnnotation(foo=anno_text(MONTH_NAMES, location=1, rot=30, just="right", gp=gpar(fontsize=[i + 5 for i in range(12)]))),
        "ch03_anno_text_rot.png",
    )
    save_anno(
        rowAnnotation(foo=anno_text(MONTH_NAMES, location=0.5, just="center")),
        "ch03_anno_text_center.png",
    )
    save_anno(
        rowAnnotation(
            foo=anno_text(
                MONTH_NAMES,
                location=0.5,
                just="center",
                gp=gpar(fill=np.repeat([2, 3, 4], 4), col="white", border="black"),
                width=max_text_width(MONTH_NAMES) * 1.2,
            )
        ),
        "ch03_anno_text_bg.png",
    )

    rich_txt = [f"<span style='color:red'>**{ch_l}**</span> *{ch_l.lower()}* <span style='color:blue'>{ch_l}</span>" for ch_l in "ABCDEFGHIJ"]
    save_anno(
        rowAnnotation(
            foo=anno_text(
                gt_render(rich_txt, align_widths=True, r=unit(2, "pt"), padding=unit([2, 2, 2, 2], "pt")),
                gp=gpar(box_col="blue", box_lwd=2),
                just="right",
                location=unit(1, "npc"),
            )
        ),
        "ch03_anno_text_gt_render.png",
    )

    m = rng.standard_normal((10, 10))
    save_ht(
        Heatmap(m)
        + rowAnnotation(
            month=anno_text(MONTH_NAMES[:10], just="center", location=unit(0.5, "npc"), show_name=True),
            annotation_name_rot=0,
        ),
        "ch03_anno_text_show_name.png",
        width=6.6,
        height=4.8,
    )
    reset_state()


def render_3_15_numeric() -> None:
    rng = np.random.default_rng(123)
    x = np.round(rng.uniform(size=10), 3)
    idx = list(range(10))

    save_anno(rowAnnotation(numeric=anno_numeric(x), annotation_name_rot=0), "ch03_anno_numeric_basic.png", index=idx)
    save_anno(
        rowAnnotation(numeric=anno_numeric(x, bg_gp=gpar(fill="orange", col="black")), annotation_name_rot=0),
        "ch03_anno_numeric_bg.png",
        index=idx,
    )
    save_anno(
        rowAnnotation(numeric=anno_numeric(x, bg_gp=gpar(fill="orange", col="black"), align_to="right"), annotation_name_rot=0),
        "ch03_anno_numeric_right.png",
        index=idx,
    )
    save_anno(
        rowAnnotation(
            numeric=anno_numeric(x, labels_format=lambda v: [f"{val * 100:.2f}%" for val in v]),
            annotation_name_rot=0,
        ),
        "ch03_anno_numeric_pct.png",
        index=idx,
    )

    x_div = np.round(rng.uniform(-1, 1, 10), 3)
    div_cols = ["green" if v < 0 else "red" for v in x_div]
    save_anno(
        rowAnnotation(numeric=anno_numeric(x_div, bg_gp=gpar(fill=div_cols, col=div_cols)), annotation_name_rot=0),
        "ch03_anno_numeric_diverging.png",
        index=idx,
    )

    x_pv = 10.0 ** (-rng.uniform(1, 6, 10))
    save_anno(
        rowAnnotation(
            pvalue=anno_numeric(
                x_pv,
                rg=(float(x_pv.min()), 1.0),
                x_convert=lambda v: -np.log10(v),
                labels_format=lambda v: [f"{val:.2e}" for val in v],
            ),
            annotation_name_rot=0,
        ),
        "ch03_anno_numeric_pvalue.png",
        index=idx,
    )


def render_3_16_customize(icon_path: str) -> None:
    rng = np.random.default_rng(123)
    x = ["a", "a", "a", "b", "b", "b", "c", "c", "d", "d"]
    icon_arr = np.asarray(Image.open(icon_path).convert("RGBA"))
    graphics = {
        "a": lambda px, py, w, h: grid_rect(px, py, w * 0.8, h * 0.33, gp=gpar(fill="red")),
        "b": lambda px, py, w, h: grid_text("A", px, py, gp=gpar(col="darkgreen")),
        "c": lambda px, py, w, h: grid_points(px, py, gp=gpar(col="orange"), pch=16),
        "d": lambda px, py, w, h: grid_raster(icon_arr, px, py, width=unit(0.8, "snpc"), height=unit(0.8, "snpc") * icon_arr.shape[0] / icon_arr.shape[1]),
    }

    m = rng.standard_normal((10, 10))
    save_ht(
        Heatmap(
            m,
            name="mat",
            top_annotation=HeatmapAnnotation(foo=anno_customize(x, graphics=graphics)),
            right_annotation=rowAnnotation(bar=anno_customize(x, graphics=graphics)),
        ),
        "ch03_anno_customize_basic.png",
        width=6.4,
        height=5.2,
        padding=unit([8, 2, 2, 2], "mm"),
    )
    reset_state()

    save_ht(
        Heatmap(
            m,
            name="mat",
            top_annotation=HeatmapAnnotation(foo=anno_customize(x, graphics=graphics)),
            right_annotation=rowAnnotation(bar=anno_customize(x, graphics=graphics)),
            column_split=x,
            row_split=x,
        ),
        "ch03_anno_customize_split.png",
        width=6.6,
        height=5.4,
        padding=unit([8, 2, 2, 2], "mm"),
    )
    reset_state()

    ht = Heatmap(m, name="mat", top_annotation=HeatmapAnnotation(foo=anno_customize(x, graphics=graphics)))
    lgd = Legend(title="foo", at=list(graphics.keys()), graphics=graphics)
    save_ht(ht, "ch03_anno_customize_legend.png", width=6.4, height=5.0, annotation_legend_list=lgd)
    reset_state()


def render_3_17_mark() -> None:
    rng = np.random.default_rng(123)
    m = pd.DataFrame(rng.standard_normal((100, 10)), index=[str(i) for i in range(1, 101)])
    ha = rowAnnotation(foo=anno_mark(at=[0, 1, 2, 3, 19, 59, 96, 97, 98, 99], labels=MONTH_NAMES[:10]))
    ht1 = Heatmap(m, name="mat", cluster_rows=False, row_names_side="left", row_names_gp=gpar(fontsize=4)) + ha
    ht2 = Heatmap(m, name="mat", cluster_rows=False, row_names_side="left", row_names_gp=gpar(fontsize=4), row_km=4) + ha

    g1 = grid_grabExpr(lambda: ch.draw(ht1, newpage=False), width_in=4.8, height_in=5.5)
    g2 = grid_grabExpr(lambda: ch.draw(ht2, newpage=False), width_in=4.8, height_in=5.5)

    def _draw_pair() -> None:
        pushViewport(Viewport(x=0.25, width=0.5))
        grid_draw(g1)
        popViewport()
        pushViewport(Viewport(x=0.75, width=0.5))
        grid_draw(g2)
        popViewport()

    save_Standalone(_draw_pair, "ch03_anno_mark.png", width=9.6, height=5.5)
    reset_state()


def render_3_18_zoom() -> None:
    rng = np.random.default_rng(123)
    m = rng.standard_normal((100, 10))
    subgroup = rng.choice(["a", "b", "c"], 100, p=[1 / 16, 5 / 16, 10 / 16])
    rg = (float(m.min()), float(m.max()))

    def panel_fun(index, nm):
        pushViewport(Viewport(xscale=rg, yscale=(0, 2)))
        grid_rect()
        ch.grid_xaxis(gp=gpar(fontsize=8))
        grid_boxplot(m[index, :].ravel(), pos=1, direction="horizontal")
        popViewport()

    align_cat = {k: [i for i, v in enumerate(subgroup) if v == k] for k in ("a", "b", "c")}
    anno1 = anno_link(align_to=align_cat, which="row", panel_fun=panel_fun, size=unit(1.5, "cm"), gap=unit(1, "cm"), width=unit(4, "cm"))
    save_ht(
        Heatmap(m, name="mat", row_split=subgroup) + rowAnnotation(foo=anno1),
        "ch03_anno_zoom_cat.png",
        width=7.5,
        height=6.4,
        padding=unit([12, 2, 12, 2], "mm"),
    )
    reset_state()

    align_sub = {k: [i for i, v in enumerate(subgroup) if v == k] for k in ("a", "b")}
    anno2 = anno_link(align_to=align_sub, which="row", panel_fun=panel_fun, size=unit(1.5, "cm"), gap=unit(1, "cm"), width=unit(4, "cm"))
    save_ht(
        Heatmap(m, name="mat", row_split=subgroup) + rowAnnotation(foo=anno2),
        "ch03_anno_zoom_subset.png",
        width=7.5,
        height=6.4,
        padding=unit([12, 2, 12, 2], "mm"),
    )
    reset_state()

    anno3 = anno_link(align_to=[list(range(20))], which="row", panel_fun=panel_fun, size=unit(2, "cm"), gap=unit(1, "cm"), width=unit(4, "cm"))
    save_ht(
        Heatmap(m, name="mat", right_annotation=rowAnnotation(foo=anno3), cluster_rows=False),
        "ch03_anno_zoom_nosplit.png",
        width=7.0,
        height=5.6,
    )
    reset_state()


def render_3_19_textbox() -> None:
    rng = np.random.default_rng(123)
    letters = list("abcdefghijklmnopqrstuvwxyz")

    def random_text(n: int, n_words: int = 1) -> list[str]:
        out = []
        for _ in range(n):
            k = int(rng.integers(1, n_words + 1)) if n_words > 1 else 1
            ws = ["".join(rng.choice(letters, int(rng.integers(4, 11)))) for _ in range(k)]
            out.append(" ".join(ws))
        return out

    words = random_text(10)
    fs_sorted = np.sort(rng.uniform(5, 20, 10))
    panels_words = [
        ("Default", textbox_grob(words)),
        ("gp=gpar(col=1..10)", textbox_grob(words, gp=gpar(col=list(range(1, 11))))),
        ("Varying fontsize", textbox_grob(words, gp=gpar(fontsize=rng.uniform(5, 20, 10)))),
        ("Background", textbox_grob(words, background_gp=gpar(fill="#CCCCCC", col="#808080"))),
        ("round_corners=True", textbox_grob(words, background_gp=gpar(fill="#CCCCCC", col=None), round_corners=True)),
        ("padding=4mm", textbox_grob(words, background_gp=gpar(fill="#CCCCCC", col="#808080"), padding=unit(4, "mm"))),
        ("line_space=5mm", textbox_grob(words, line_space=unit(5, "mm"))),
        ("text_space=5mm", textbox_grob(words, text_space=unit(5, "mm"))),
        ("first_text_from='bottom'", textbox_grob(words, gp=gpar(fontsize=fs_sorted), first_text_from="bottom")),
    ]

    def _draw_words_grid() -> None:
        ncols, nrows = 3, 3
        for idx, (title, gb) in enumerate(panels_words):
            r, col_i = divmod(idx, ncols)
            pushViewport(Viewport(x=unit((col_i + 0.5) / ncols, "npc"), y=unit(1.0 - (r + 0.5) / nrows, "npc"), width=unit(1.0 / ncols, "npc"), height=unit(1.0 / nrows, "npc")))
            grid_text(title, x=0.5, y=unit(0.92, "npc"), just="top", gp=gpar(fontsize=9, fontface="bold"))
            pushViewport(Viewport(y=unit(0.44, "npc"), height=unit(0.75, "npc")))
            grid_draw(gb)
            popViewport()
            popViewport()

    save_Standalone(_draw_words_grid, "ch03_textbox_grob_words.png", width=9.0, height=5.8)
    reset_state()

    sentences = random_text(5, 8)
    fs_s = rng.uniform(6, 15, 5)
    panels_sent = [
        ("Default sentences", textbox_grob(sentences, gp=gpar(col=list(range(1, 6)), fontsize=fs_s))),
        ("word_wrap=True", textbox_grob(sentences, gp=gpar(col=list(range(1, 6)), fontsize=fs_s), word_wrap=True, max_width=unit(8, "cm"))),
        ("add_new_line=True", textbox_grob(sentences, gp=gpar(col=list(range(1, 6)), fontsize=fs_s), add_new_line=True)),
        ("add_new_line + word_wrap", textbox_grob(sentences, gp=gpar(col=list(range(1, 6)), fontsize=fs_s), add_new_line=True, word_wrap=True, max_width=unit(8, "cm"))),
    ]

    def _draw_sent_grid() -> None:
        ncols, nrows = 2, 2
        for idx, (title, gb) in enumerate(panels_sent):
            r, col_i = divmod(idx, ncols)
            pushViewport(Viewport(x=unit((col_i + 0.5) / ncols, "npc"), y=unit(1.0 - (r + 0.5) / nrows, "npc"), width=unit(1.0 / ncols, "npc"), height=unit(1.0 / nrows, "npc")))
            grid_text(title, x=0.5, y=unit(0.94, "npc"), just="top", gp=gpar(fontsize=10, fontface="bold"))
            pushViewport(Viewport(y=unit(0.45, "npc"), height=unit(0.80, "npc")))
            grid_draw(gb)
            popViewport()
            popViewport()

    save_Standalone(_draw_sent_grid, "ch03_textbox_grob_sentences.png", width=9.2, height=6.4)
    reset_state()

    mat = rng.standard_normal((100, 10))
    split = np.repeat(list("abcdefghij"), 10)
    text_d = {k: random_text(10) for k in list("abcdefghij")}

    save_ht(
        Heatmap(mat, name="mat", cluster_rows=False, row_split=split)
        + rowAnnotation(textbox=anno_textbox(split, text_d, max_width=unit(9.0, "cm"))),
        "ch03_anno_textbox_right.png",
        width=9.2,
        height=8.2,
        padding=unit([6, 2, 6, 2], "mm"),
    )
    reset_state()

    save_ht(
        rowAnnotation(textbox=anno_textbox(split, text_d, side="left", max_width=unit(9.0, "cm")))
        + Heatmap(mat, name="mat", cluster_rows=False, row_split=split, row_title=None),
        "ch03_anno_textbox_left.png",
        width=9.2,
        height=8.2,
        padding=unit([6, 2, 6, 2], "mm"),
    )
    reset_state()

    split5 = np.repeat(list("abcde"), 20)
    sent_d = {k: random_text(3, 5) for k in list("abcde")}
    save_ht(
        Heatmap(mat, name="mat", row_split=split5)
        + rowAnnotation(textbox=anno_textbox(split5, sent_d, word_wrap=True, add_new_line=True, max_width=unit(10.0, "cm"))),
        "ch03_anno_textbox_sentences.png",
        width=9.4,
        height=7.4,
        padding=unit([6, 2, 6, 2], "mm"),
    )
    reset_state()

    align5 = {k: [i for i, v in enumerate(split5) if v == k] for k in ("a", "b")}
    save_ht(
        Heatmap(mat, name="mat", row_split=split5)
        + rowAnnotation(
            textbox=anno_textbox(align5, {k: sent_d[k] for k in ("a", "b")}, word_wrap=True, add_new_line=True, max_width=unit(10.0, "cm"))
        ),
        "ch03_anno_textbox_align_subset.png",
        width=9.4,
        height=7.4,
        padding=unit([6, 2, 6, 2], "mm"),
    )
    reset_state()

    save_ht(
        Heatmap(
            mat,
            name="mat",
            cluster_rows=False,
            right_annotation=rowAnnotation(
                textbox=anno_textbox(
                    {"a": list(range(10)), "b": list(range(10, 30))},
                    {k: sent_d[k] for k in ("a", "b")},
                    word_wrap=True,
                    add_new_line=True,
                )
            ),
        ),
        "ch03_anno_textbox_nosplit.png",
        width=8.5,
        height=6.0,
    )
    reset_state()

    split_rep = np.tile(list("abcdefghij"), 10)
    text_rep = {k: random_text(10) for k in list("abcdefghij")}
    save_ht(
        Heatmap(
            mat,
            name="mat",
            row_split=split_rep,
            right_annotation=rowAnnotation(textbox=anno_textbox(split_rep, text_rep, by="anno_block")),
        ),
        "ch03_anno_textbox_block.png",
        width=8.5,
        height=6.8,
    )
    reset_state()

    save_ht(
        Heatmap(
            mat,
            name="mat",
            cluster_rows=False,
            right_annotation=rowAnnotation(
                textbox=anno_textbox({"a": list(range(10)), "b": list(range(49, 70))}, {k: text_rep[k] for k in ("a", "b")}, by="anno_block")
            ),
        ),
        "ch03_anno_textbox_block_subset.png",
        width=8.5,
        height=6.2,
    )
    reset_state()

    split_df = np.tile(list("abcde"), 20)
    text_df_dict = {}
    for k in list("abcde"):
        df_k = pd.DataFrame({"text": random_text(10), "fontsize": rng.uniform(6, 11, 10)})
        df_k["col"] = [rand_color(1, friendly=True)[0]] * 10 if rng.uniform() > 0.5 else list(range(1, 11))
        text_df_dict[k] = df_k
    save_ht(
        Heatmap(mat, name="mat", cluster_rows=False, row_split=split_df)
        + rowAnnotation(textbox=anno_textbox(split_df, text_df_dict, max_width=unit(11.5, "cm"))),
        "ch03_anno_textbox_df.png",
        width=9.6,
        height=7.4,
        padding=unit([6, 2, 6, 2], "mm"),
    )
    reset_state()


def render_3_20_summary() -> None:
    rng = np.random.default_rng(123)
    ha_disc = HeatmapAnnotation(summary=anno_summary(height=unit(4, "cm")))
    v_cat = rng.choice(["a", "b"], 50)
    sp = rng.choice(["a", "b"], 50)
    save_ht(
        Heatmap(v_cat, name="mat", col={"a": "red", "b": "blue"}, top_annotation=ha_disc, width=unit(2, "cm"), row_split=sp),
        "ch03_anno_summary_discrete.png",
        width=5.0,
        height=5.8,
    )
    reset_state()

    ha_cont = HeatmapAnnotation(summary=anno_summary(gp=gpar(fill=[2, 3]), height=unit(4, "cm")))
    v_num = rng.standard_normal(50)
    save_ht(
        Heatmap(v_num, name="mat", top_annotation=ha_cont, width=unit(2, "cm"), row_split=sp),
        "ch03_anno_summary_continuous.png",
        width=5.0,
        height=5.8,
    )
    reset_state()

    m = rng.standard_normal((50, 10))
    ht_list = Heatmap(m, name="main_matrix")
    ht_list = ht_list + Heatmap(v_cat, name="mat1", top_annotation=HeatmapAnnotation(summary=anno_summary(height=unit(3, "cm"))), width=unit(1, "cm"))
    ht_list = ht_list + Heatmap(v_num, name="mat2", top_annotation=HeatmapAnnotation(summary=anno_summary(gp=gpar(fill=[2, 3]), height=unit(3, "cm"))), width=unit(1, "cm"))
    lgd_bp = Legend(labels=["group a", "group b"], title="group", legend_gp=gpar(fill=["red", "blue"]))
    save_ht(
        ht_list,
        "ch03_anno_summary_combined.png",
        width=7.5,
        height=6.0,
        row_split=sp,
        ht_gap=unit(5, "mm"),
        heatmap_legend_list=[lgd_bp],
    )
    reset_state()


def render_3_21_multiple() -> None:
    rng = np.random.default_rng(123)
    mat2 = np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)])

    ha_leg = HeatmapAnnotation(
        foo=np.arange(1, 11),
        bar=mat2,
        pt=anno_points(np.arange(1, 11)),
        show_legend={"bar": False},
    )
    save_ht(Heatmap(rng.standard_normal((10, 10)), name="mat", top_annotation=ha_leg), "ch03_multi_show_legend.png", width=6.4, height=5.2)
    reset_state()

    save_anno(
        HeatmapAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)), gp=gpar(col="red")),
        "ch03_multi_gp.png",
    )
    save_anno(
        HeatmapAnnotation(
            foo=np.arange(1, 11),
            bar=mat2,
            pt=anno_points(np.arange(1, 11)),
            show_annotation_name={"bar": False},
            border={"foo": True},
        ),
        "ch03_multi_show_name_border.png",
    )
    save_anno(
        HeatmapAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)), gap=unit(2, "mm")),
        "ch03_multi_gap_scalar.png",
    )
    save_anno(
        HeatmapAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)), gap=unit([2, 10], "mm")),
        "ch03_multi_gap_vector.png",
    )

    # 3.21.2 Size of annotations
    save_anno(
        HeatmapAnnotation(foo=mat2, bar=np.arange(1, 11), pt=anno_points(np.arange(1, 11))),
        "ch03_size_default.png",
    )
    save_anno(
        HeatmapAnnotation(foo=mat2, bar=np.arange(1, 11), pt=anno_points(np.arange(1, 11)), height=unit(6, "cm")),
        "ch03_size_height.png",
    )
    save_anno(
        HeatmapAnnotation(
            foo=mat2,
            bar=np.arange(1, 11),
            pt=anno_points(np.arange(1, 11)),
            simple_anno_size=unit(1, "cm"),
            height=unit(6, "cm"),
        ),
        "ch03_size_simple_anno_size.png",
    )
    save_anno(
        HeatmapAnnotation(
            foo=mat2,
            bar=np.arange(1, 11),
            pt=anno_points(np.arange(1, 11)),
            annotation_height=unit([1, 2, 3], "cm"),
        ),
        "ch03_size_anno_height_abs.png",
    )
    save_anno(
        HeatmapAnnotation(
            foo=mat2,
            bar=np.arange(1, 11),
            pt=anno_points(np.arange(1, 11)),
            annotation_height=[1, 2, 3],
            height=unit(6, "cm"),
        ),
        "ch03_size_anno_height_rel.png",
    )
    save_anno(
        HeatmapAnnotation(foo=mat2, bar=np.arange(1, 11), height=unit(6, "cm")),
        "ch03_size_simple_only_no_adj.png",
    )
    save_anno(
        HeatmapAnnotation(foo=mat2, bar=np.arange(1, 11), height=unit(6, "cm"), simple_anno_size_adjust=True),
        "ch03_size_simple_only_adj.png",
    )

    # 3.21.3 Annotation labels
    save_anno(
        HeatmapAnnotation(
            foo=np.arange(1, 11),
            bar=mat2,
            pt=anno_points(np.arange(1, 11)),
            annotation_label=["Annotation_foo", "Annotation_bar", "Annotation_pt"],
        ),
        "ch03_anno_labels_custom.png",
        width=7.8,
        margin_cm=4.5,
    )

    ha_gt = HeatmapAnnotation(
        foo=np.arange(1, 11),
        bar=mat2,
        pt=anno_points(np.arange(1, 11)),
        annotation_label=[
            gt_render("**Anno**_<span style='color:red'>foo</span>", padding=unit([2, 2, 2, 2], "pt")),
            gt_render("**Anno**_<span style='color:blue'>bar</span>", padding=unit([2, 2, 2, 2], "pt")),
            gt_render("**Anno**_<span style='color:green'>pt</span>", padding=unit([2, 2, 2, 2], "pt")),
        ],
        annotation_name_gp=gpar(box_fill="lightgrey"),
    )
    save_ht(Heatmap(rng.standard_normal((10, 10)), name="mat", top_annotation=ha_gt), "ch03_anno_labels_gt_render.png", width=6.6, height=5.0)
    reset_state()

    ha_rot_col = HeatmapAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)), annotation_name_rot=45)
    save_ht(
        Heatmap(rng.standard_normal((10, 10)), name="mat", top_annotation=ha_rot_col),
        "ch03_anno_labels_rot_col.png",
        width=6.6,
        height=5.0,
        padding=unit([2, 2, 2, 8], "mm"),
    )
    reset_state()

    ha_rot_row = rowAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)), annotation_name_rot=45)
    save_ht(
        Heatmap(rng.standard_normal((10, 10)), name="mat", left_annotation=ha_rot_row),
        "ch03_anno_labels_rot_row.png",
        width=6.4,
        height=5.2,
        padding=unit([12, 2, 2, 2], "mm"),
    )
    reset_state()


def render_3_22_utils() -> None:
    rng = np.random.default_rng(123)
    mat2 = np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)])
    ha1 = HeatmapAnnotation(foo=np.arange(1, 11), bar=mat2, pt=anno_points(np.arange(1, 11)))
    ha2 = HeatmapAnnotation(FOO=rng.uniform(size=10), BAR=rng.choice(["a", "b"], 10), PT=anno_points(rng.standard_normal(10)))
    ha = c(ha1, ha2)
    save_anno(ha, "ch03_util_concat.png")

    ha_sub = ha[:5, ["foo", "PT"]]
    save_anno(ha_sub, "ch03_util_subset.png")


def render_3_23_custom_af() -> None:
    rng = np.random.default_rng(123)
    x = np.arange(1, 11, dtype=float)

    def _lollipop(index, k=1, n=1):
        m_len = len(index)
        pushViewport(Viewport(xscale=(0.5, m_len + 0.5), yscale=(0, 10)))
        grid_rect()
        grid_points(np.arange(1, m_len + 1), x[index], default_units="native")
        grid_segments(np.arange(1, m_len + 1), 0, np.arange(1, m_len + 1), x[index], default_units="native")
        if k == 1:
            ch.grid_yaxis()
        popViewport()

    anno1 = AnnotationFunction(
        fun=_lollipop,
        var_import={"x": x},
        n=10,
        subsettable=True,
        height=unit(2, "cm"),
    )
    m = np.vstack([np.arange(1, 11), np.arange(11, 21)])
    ht1 = Heatmap(m, name="mat1", top_annotation=HeatmapAnnotation(foo=anno1))
    ht2 = Heatmap(m, name="mat2", top_annotation=HeatmapAnnotation(foo=anno1), column_split=np.repeat(["A", "B"], 5))

    g1 = grid_grabExpr(lambda: ch.draw(ht1, newpage=False), width_in=6.4, height_in=3.4)
    g2 = grid_grabExpr(lambda: ch.draw(ht2, newpage=False), width_in=6.4, height_in=3.4)

    def _draw_stacked() -> None:
        pushViewport(Viewport(y=1.0, height=0.5, just="top"))
        grid_draw(g1)
        popViewport()
        pushViewport(Viewport(y=0.5, height=0.5, just="top"))
        grid_draw(g2)
        popViewport()

    save_Standalone(_draw_stacked, "ch03_custom_annotation_fun.png", width=6.4, height=6.8)
    reset_state()

    ht_empty = Heatmap(m, name="mat", top_annotation=HeatmapAnnotation(foo=anno_empty(height=unit(2, "cm"))), column_split=np.repeat(["A", "B"], 5))

    def _dec_lollipop(res) -> None:
        co = column_order(res)
        with decorate_annotation("foo", slice=1):
            od = co[0]
            pushViewport(Viewport(xscale=(0.5, len(od) + 0.5), yscale=(0, 10)))
            grid_points(np.arange(1, len(od) + 1), x[od], default_units="native")
            grid_segments(np.arange(1, len(od) + 1), 0, np.arange(1, len(od) + 1), x[od], default_units="native")
            ch.grid_yaxis()
            popViewport()
        with decorate_annotation("foo", slice=2):
            od = co[1]
            pushViewport(Viewport(xscale=(0.5, len(od) + 0.5), yscale=(0, 10)))
            grid_points(np.arange(1, len(od) + 1), x[od], default_units="native")
            grid_segments(np.arange(1, len(od) + 1), 0, np.arange(1, len(od) + 1), x[od], default_units="native")
            popViewport()

    save_decorated(ht_empty, "ch03_custom_anno_empty_reimpl.png", _dec_lollipop, width=6.4, height=3.8)
    reset_state()

    def anno_pct(xv: np.ndarray) -> AnnotationFunction:
        max_x = float(np.max(xv))
        text_lbl = [f"{val * 100:.2f}%" for val in xv]

        def cell_fun_pct(i: int) -> None:
            pushViewport(Viewport(xscale=(0, max_x)))
            grid_roundrect(
                x=unit(1, "npc"),
                width=unit(float(xv[i]), "native"),
                height=unit(1, "npc") - unit(4, "pt"),
                just="right",
                gp=gpar(fill="#0000FF80", col=None),
            )
            grid_text(text_lbl[i], x=unit(1, "npc"), just="right")
            popViewport()

        return AnnotationFunction(
            cell_fun=cell_fun_pct,
            var_import={"max_x": max_x, "x": xv, "text": text_lbl},
            which="row",
            n=len(xv),
            width=max_text_width(text_lbl) * 1.25,
        )

    xv = rng.uniform(size=10)
    ha_pct = rowAnnotation(foo=anno_pct(xv), annotation_name_rot=0)
    m_pct = pd.DataFrame(rng.standard_normal((10, 10)), index=[f"{val:.4f}" for val in xv])
    save_ht(ha_pct + Heatmap(m_pct, name="mat"), "ch03_custom_cell_fun_pct.png", width=6.4, height=5.0)
    reset_state()


def main() -> None:
    reset_state()
    with tempfile.TemporaryDirectory() as tmpdir:
        icon_paths = _generate_demo_icons(Path(tmpdir))
        render_intro()
        render_3_1_simple()
        render_3_2_anno_simple()
        render_3_3_empty()
        render_3_4_block()
        render_3_5_image(icon_paths)
        render_3_6_points()
        render_3_7_lines()
        render_3_8_barplot()
        render_3_9_boxplot()
        render_3_10_histogram()
        render_3_11_density()
        render_3_12_joyplot()
        render_3_13_horizon()
        render_3_14_text()
        render_3_15_numeric()
        render_3_16_customize(icon_paths[5])
        render_3_17_mark()
        render_3_18_zoom()
        render_3_19_textbox()
        render_3_20_summary()
        render_3_21_multiple()
        render_3_22_utils()
        render_3_23_custom_af()
    print("Chapter 3 figures generated successfully.")


if __name__ == "__main__":
    main()
