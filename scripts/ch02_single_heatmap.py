"""Generate all Chapter 2 figures: A Single Heatmap (Sections 2.1 - 2.13)."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any, Callable

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, optimal_leaf_ordering

import complexheatmap as ch
from complexheatmap import (
    Dendrogram,
    Heatmap,
    HeatmapAnnotation,
    Legend,
    Viewport,
    anno_barplot,
    anno_block,
    anno_points,
    as_dendrogram,
    cluster_within_group,
    colorRamp2,
    color_branches,
    column_dend,
    column_order,
    convertHeight,
    convertWidth,
    dist,
    draw,
    gpar,
    grid_circle,
    grid_draw,
    grid_grabExpr,
    grid_points,
    grid_rect,
    grid_segments,
    grid_text,
    gt_render,
    hcl_colors,
    hclust,
    ht_opt,
    ht_size,
    max_text_width,
    pindex,
    popViewport,
    pushViewport,
    rand_color,
    restore_matrix,
    rowAnnotation,
    row_dend,
    row_order,
    unit,
)

from _common import reset_state, save_Standalone, save_ht


def _build_main_matrix() -> pd.DataFrame:
    """Create the canonical 18x24 3-group block matrix used throughout Chapter 2."""
    rng = np.random.default_rng(123)
    nr1, nr2, nr3 = 4, 8, 6
    nc1, nc2, nc3 = 6, 8, 10
    nr = nr1 + nr2 + nr3
    nc = nc1 + nc2 + nc3

    b1 = np.vstack([
        rng.normal(1.0, 0.5, (nr1, nc1)),
        rng.normal(0.0, 0.5, (nr2, nc1)),
        rng.normal(0.0, 0.5, (nr3, nc1)),
    ])
    b2 = np.vstack([
        rng.normal(0.0, 0.5, (nr1, nc2)),
        rng.normal(1.0, 0.5, (nr2, nc2)),
        rng.normal(0.0, 0.5, (nr3, nc2)),
    ])
    b3 = np.vstack([
        rng.normal(0.5, 0.5, (nr1, nc3)),
        rng.normal(0.5, 0.5, (nr2, nc3)),
        rng.normal(1.0, 0.5, (nr3, nc3)),
    ])
    mat = np.hstack([b1, b2, b3])
    r_perm = rng.permutation(nr)
    c_perm = rng.permutation(nc)
    mat = mat[np.ix_(r_perm, c_perm)]
    return pd.DataFrame(
        mat,
        index=[f"row{i + 1}" for i in range(nr)],
        columns=[f"column{j + 1}" for j in range(nc)],
    )


def _save_grid_panels(
    panels: list[tuple[Callable[[], Any], float, float]],
    filename: str,
    ncol: int,
    dpi: int = 140,
) -> Path:
    """Capture multiple draw callables with grid_grabExpr and tile them in a grid."""
    nrow = (len(panels) + ncol - 1) // ncol
    col_ws = [max(panels[r * ncol + c][1] for r in range(nrow) if r * ncol + c < len(panels)) for c in range(ncol)]
    row_hs = [max(panels[r * ncol + c][2] for c in range(ncol) if r * ncol + c < len(panels)) for r in range(nrow)]
    tot_w = sum(col_ws)
    tot_h = sum(row_hs)

    grobs: list[tuple[int, int, Any]] = []
    for idx, (fn, pw, ph) in enumerate(panels):
        r = idx // ncol
        c = idx % ncol
        g = grid_grabExpr(fn, width_in=pw, height_in=ph)
        grobs.append((r, c, g))
        reset_state()

    def _compose() -> None:
        y_top = 1.0
        for r in range(nrow):
            rh_frac = row_hs[r] / tot_h
            x_left = 0.0
            for c in range(ncol):
                cw_frac = col_ws[c] / tot_w
                for gr, gc, g in grobs:
                    if gr == r and gc == c:
                        pushViewport(
                            Viewport(
                                x=x_left,
                                y=y_top,
                                width=cw_frac,
                                height=rh_frac,
                                just=("left", "top"),
                            )
                        )
                        grid_draw(g)
                        popViewport()
                x_left += cw_frac
            y_top -= rh_frac

    out = save_Standalone(_compose, filename, width=tot_w, height=tot_h, dpi=dpi)
    reset_state()
    return out


def render_intro(mat: pd.DataFrame) -> None:
    ht = Heatmap(mat)
    save_ht(ht, "ch02_default.png", width=6.0, height=5.2)
    reset_state()


def render_2_1_colors(mat: pd.DataFrame) -> None:
    col_fun = colorRamp2([-2, 0, 2], ["green", "white", "red"])
    save_ht(Heatmap(mat, name="mat", col=col_fun), "ch02_col_ramp2.png", width=6.0, height=5.2)
    reset_state()

    mat2 = mat.copy()
    mat2.iloc[0, 0] = 100000.0
    save_ht(
        Heatmap(mat2, name="mat", col=col_fun, column_title="a matrix with outliers"),
        "ch02_col_outlier.png",
        width=6.0,
        height=5.2,
    )
    reset_state()

    _save_grid_panels(
        [
            (lambda: draw(Heatmap(mat, name="mat", col=col_fun, column_title="mat"), newpage=False), 5.0, 4.4),
            (lambda: draw(Heatmap(mat / 4.0, name="mat", col=col_fun, column_title="mat / 4"), newpage=False), 5.0, 4.4),
            (lambda: draw(Heatmap(np.abs(mat), name="mat", col=col_fun, column_title="abs(mat)"), newpage=False), 5.0, 4.4),
        ],
        "ch02_col_compare.png",
        ncol=2,
    )

    rainbow_cols = [
        matplotlib.colors.to_hex(matplotlib.colormaps["rainbow"](x))
        for x in np.linspace(1.0, 0.0, 10)
    ]
    save_ht(
        Heatmap(mat, name="mat", col=rainbow_cols, column_title="set a color vector for a continuous matrix"),
        "ch02_col_vector.png",
        width=6.0,
        height=5.2,
    )
    reset_state()

    rng = np.random.default_rng(123)
    discrete_num = rng.integers(1, 5, size=(10, 10))
    colors_num = {1: "black", 2: "#DF536B", 3: "#61D04F", 4: "#2297E6"}
    save_ht(
        Heatmap(discrete_num, name="mat", col=colors_num, column_title="a discrete numeric matrix"),
        "ch02_col_discrete_num.png",
        width=5.6,
        height=4.8,
    )
    reset_state()

    discrete_char = rng.choice(["a", "b", "c", "d"], size=(10, 10))
    colors_char = {"a": "black", "b": "#DF536B", "c": "#61D04F", "d": "#2297E6"}
    save_ht(
        Heatmap(discrete_char, name="mat", col=colors_char, column_title="a discrete character matrix"),
        "ch02_col_discrete_char.png",
        width=5.6,
        height=4.8,
    )
    reset_state()

    mat_na = mat.copy()
    na_mask = rng.choice([True, False], size=mat_na.shape, p=[0.1, 0.9])
    mat_na = mat_na.mask(na_mask)
    save_ht(
        Heatmap(mat_na, name="mat", na_col="black", column_title="a matrix with NA values"),
        "ch02_col_na.png",
        width=6.0,
        height=5.2,
    )
    reset_state()

    vmin, vmax = float(mat.min().min()), float(mat.max().max())
    brks = np.linspace(vmin, vmax, 3)
    f_lab = colorRamp2(brks, ["blue", "#EEEEEE", "red"], space="LAB")
    f_rgb = colorRamp2(brks, ["blue", "#EEEEEE", "red"], space="RGB")
    _save_grid_panels(
        [
            (lambda: draw(Heatmap(mat, name="mat1", col=f_lab, column_title="LAB color space"), newpage=False), 5.2, 4.6),
            (lambda: draw(Heatmap(mat, name="mat2", col=f_rgb, column_title="RGB color space"), newpage=False), 5.2, 4.6),
        ],
        "ch02_col_space_lab_rgb.png",
        ncol=2,
    )

    spaces = ["RGB", "LAB", "XYZ", "HSV"]
    palettes = [
        ("green-black-red", ["green", "black", "red"]),
        ("blue-white-red", ["blue", "white", "red"]),
    ]

    def _draw_spaces() -> None:
        vals = np.linspace(-1.0, 1.0, 120)
        x_centers = (np.arange(len(vals)) + 0.5) / len(vals)
        w_cell = 1.0 / len(vals)
        for r_idx, (pal_name, cols) in enumerate(palettes):
            y_top = 1.0 - r_idx * 0.5
            for c_idx, sp in enumerate(spaces):
                f_sp = colorRamp2([-1.0, 0.0, 1.0], cols, space=sp)
                c_vals = list(f_sp(vals))
                pushViewport(
                    Viewport(
                        x=c_idx / len(spaces),
                        y=y_top,
                        width=1.0 / len(spaces),
                        height=0.5,
                        just=("left", "top"),
                    )
                )
                grid_text(f"{sp} ({pal_name})", x=0.5, y=0.86, gp=gpar(fontsize=10.5, fontface="bold"))
                pushViewport(Viewport(x=0.5, y=0.45, width=0.88, height=0.52))
                grid_rect(
                    x=x_centers,
                    y=0.5,
                    width=w_cell * 1.02,
                    height=1.0,
                    gp=gpar(fill=c_vals, col=None),
                )
                grid_rect(gp=gpar(fill="transparent", col="#333333", lwd=1.0))
                popViewport()
                popViewport()

    save_Standalone(_draw_spaces, "ch02_col_spaces_strip.png", width=10.5, height=3.4, dpi=140)
    reset_state()

    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", border=True, border_gp=gpar(col="black", lty=2), column_title="set heatmap borders"),
                    newpage=False,
                ),
                4.8,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", rect_gp=gpar(col="white", lwd=2), column_title="set cell borders"),
                    newpage=False,
                ),
                4.8,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", rect_gp=gpar(type="none"), column_title="rect_gp = gpar(type='none')"),
                    newpage=False,
                ),
                4.8,
                4.4,
            ),
        ],
        "ch02_col_borders.png",
        ncol=3,
    )


def render_2_2_titles(mat: pd.DataFrame) -> None:
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", column_title="I am a column title", row_title="I am a row title"),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", column_title="I am a column title at the bottom", column_title_side="bottom"),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        column_title="I am a big column title",
                        column_title_gp=gpar(fontsize=20, fontface="bold"),
                    ),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", row_title="row title", row_title_rot=0),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
        ],
        "ch02_titles_basic.png",
        ncol=2,
    )

    def _draw_boxed() -> None:
        ht_opt.TITLE_PADDING = unit([8.5, 8.5], "points")
        draw(
            Heatmap(
                mat,
                name="mat",
                column_title="I am a column title",
                column_title_gp=gpar(fill="red", col="white", border="blue"),
            ),
            newpage=False,
        )
        ht_opt(RESET=True)

    def _draw_gt_title() -> None:
        ht_opt.TITLE_PADDING = unit([10, 10], "points")
        draw(
            Heatmap(
                mat,
                name="mat",
                column_title=gt_render(
                    "<span style='color:blue'>**Blue bold**</span> / *italic*<br>"
                    "<span style='font-size:15pt; color:black'>Large</span> rich-text title",
                    r=unit(2, "pt"),
                    padding=unit([3, 4, 3, 4], "pt"),
                ),
            ),
            newpage=False,
        )
        ht_opt(RESET=True)

    _save_grid_panels(
        [
            (_draw_boxed, 4.8, 4.6),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", column_title=r"$\hat{\beta} = (X^T X)^{-1} X^T y$"),
                    newpage=False,
                ),
                4.8,
                4.6,
            ),
            (_draw_gt_title, 5.0, 4.6),
        ],
        "ch02_titles_box_math_rich.png",
        ncol=3,
    )


def render_2_3_clustering(mat: pd.DataFrame) -> None:
    _save_grid_panels(
        [
            (lambda: draw(Heatmap(mat, name="mat", cluster_rows=False, column_title="cluster_rows = False"), newpage=False), 5.0, 4.4),
            (lambda: draw(Heatmap(mat, name="mat", show_column_dend=False, column_title="show_column_dend = False"), newpage=False), 5.0, 4.4),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_dend_side="right",
                        column_dend_side="bottom",
                        column_title="dendrogram sides",
                    ),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        column_dend_height=unit(3.0, "cm"),
                        row_dend_width=unit(3.0, "cm"),
                        column_title="custom dendrogram sizes",
                    ),
                    newpage=False,
                ),
                5.0,
                4.4,
            ),
        ],
        "ch02_cluster_general.png",
        ncol=2,
    )

    mat_out = mat.copy()
    for i in range(10):
        mat_out.iloc[i, i] = 1000.0

    def robust_dist(x: np.ndarray, y: np.ndarray) -> float:
        qx = np.quantile(x, [0.1, 0.9])
        qy = np.quantile(y, [0.1, 0.9])
        mask = (x > qx[0]) & (x < qx[1]) & (y > qy[0]) & (y < qy[1])
        x_s, y_s = x[mask], y[mask]
        return float(np.sqrt(np.sum((x_s - y_s) ** 2)))

    col_fun = colorRamp2([-2, 0, 2], ["green", "white", "red"])
    _save_grid_panels(
        [
            (lambda: draw(Heatmap(mat_out, name="mat", col=col_fun, column_title="euclidean dist"), newpage=False), 5.2, 4.6),
            (
                lambda: draw(
                    Heatmap(
                        mat_out,
                        name="mat",
                        col=col_fun,
                        clustering_distance_rows=robust_dist,
                        clustering_distance_columns=robust_dist,
                        column_title="robust_dist",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
        ],
        "ch02_cluster_dist_outliers.png",
        ncol=2,
    )

    rng = np.random.default_rng(123)
    mat_letters = rng.choice(["a", "b", "c", "d"], size=(10, 10))

    def dist_letters(x: np.ndarray, y: np.ndarray) -> float:
        ax = np.array([ord(str(ch_c)) for ch_c in x], dtype=float)
        ay = np.array([ord(str(ch_c)) for ch_c in y], dtype=float)
        return float(np.sqrt(np.sum((ax - ay) ** 2)))

    save_ht(
        Heatmap(
            mat_letters,
            name="letters",
            col={"a": 2, "b": 3, "c": 4, "d": 5},
            clustering_distance_rows=dist_letters,
            clustering_distance_columns=dist_letters,
            cell_fun=lambda j, i, x, y, w, h, col: grid_text(mat_letters[i, j], x, y, gp=gpar(fontsize=10)),
            column_title="clustering a character matrix",
        ),
        "ch02_cluster_character.png",
        width=5.8,
        height=5.0,
    )
    reset_state()

    hc_r = hclust(dist(mat), method="ward.D2")
    hc_c = hclust(dist(mat.T), method="average")
    grp = [f"G{k}" for k in fcluster(hc_c.linkage, 3, criterion="maxclust")]
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", clustering_method_rows="single", column_title="single linkage"),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", cluster_rows=hc_r, cluster_columns=hc_c, column_title="pre-computed Hclust"),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        cluster_columns=cluster_within_group(mat.to_numpy(), grp),
                        column_title="cluster_within_group",
                    ),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
        ],
        "ch02_cluster_methods.png",
        ncol=3,
    )

    row_dend_colored = color_branches(as_dendrogram(hclust(dist(mat))), k=2)
    row_dend_nodes = row_dend_colored.copy()
    rng_n = np.random.default_rng(7)
    pchs = [16, 17, 15, 18]
    node_cols = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"]

    def _style_nodes(node: Dendrogram) -> None:
        node.node_par = {
            "pch": int(pchs[int(rng_n.integers(0, len(pchs)))]),
            "col": str(node_cols[int(rng_n.integers(0, len(node_cols)))]),
            "size": unit(2.2, "mm"),
        }
        for ch_node in node.children:
            _style_nodes(ch_node)

    _style_nodes(row_dend_nodes)

    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", cluster_rows=row_dend_colored, column_title="color_branches (k = 2)"),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        cluster_rows=row_dend_colored,
                        row_dend_gp=gpar(col="red"),
                        column_title="row_dend_gp overrides",
                    ),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        cluster_rows=row_dend_nodes,
                        row_dend_width=unit(2.0, "cm"),
                        column_title="custom node symbols",
                    ),
                    newpage=False,
                ),
                5.0,
                4.5,
            ),
        ],
        "ch02_cluster_render_dend.png",
        ncol=3,
    )

    m2 = np.arange(1, 101).reshape(10, 10)
    _save_grid_panels(
        [
            (lambda: draw(Heatmap(m2, name="mat", row_dend_reorder=False, column_title="no reordering"), newpage=False), 5.0, 4.5),
            (lambda: draw(Heatmap(m2, name="mat", row_dend_reorder=True, column_title="applied reordering"), newpage=False), 5.0, 4.5),
        ],
        "ch02_cluster_reorder_dend.png",
        ncol=2,
    )


def render_2_4_manual_order(mat: pd.DataFrame) -> None:
    r_num = np.argsort([int(r.replace("row", "")) for r in mat.index])
    c_num = np.argsort([int(c.replace("column", "")) for c in mat.columns])
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", row_order=r_num, column_order=c_num, column_title="reorder by numeric index"),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_order=sorted(mat.index),
                        column_order=sorted(mat.columns),
                        column_title="reorder by row/column names",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
        ],
        "ch02_manual_order.png",
        ncol=2,
    )


def render_2_5_seriation(mat: pd.DataFrame) -> None:
    d_r = dist(mat)
    d_c = dist(mat.T)
    hc_r = hclust(d_r)
    hc_c = hclust(d_c)
    z_opt_r = optimal_leaf_ordering(hc_r.linkage, d_r.condensed)
    z_opt_c = optimal_leaf_ordering(hc_c.linkage, d_c.condensed)
    hc_opt_r = ch.Hclust(z_opt_r, n=mat.shape[0], labels=list(mat.index), method=hc_r.method)
    hc_opt_c = ch.Hclust(z_opt_c, n=mat.shape[1], labels=list(mat.columns), method=hc_c.method)

    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_order=hc_opt_r.order,
                        column_order=hc_opt_c.order,
                        column_title="seriation order (optimal leaf ordering)",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        cluster_rows=as_dendrogram(hc_opt_r),
                        cluster_columns=as_dendrogram(hc_opt_c),
                        column_title="seriation with dendrograms",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
        ],
        "ch02_seriation.png",
        ncol=2,
    )


def render_2_6_dimnames(mat: pd.DataFrame) -> None:
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_names_side="left",
                        row_dend_side="right",
                        column_names_side="top",
                        column_dend_side="bottom",
                        column_title="swapped sides",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_names_gp=gpar(col=["red"] * 10 + ["blue"] * 8, fontsize=11),
                        column_title="colored row names",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", column_names_rot=45, column_title="column_names_rot = 45 (bottom)"),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        column_names_rot=45,
                        column_names_side="top",
                        column_dend_side="bottom",
                        column_title="column_names_rot = 45 (top)",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
        ],
        "ch02_dimnames_style.png",
        ncol=2,
    )

    mat2 = mat.copy()
    new_idx = list(mat2.index)
    new_idx[0] = "very_long_feature_identifier_example_01"
    mat2.index = new_idx
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat2, name="mat", row_title="default row_names_max_width (6 cm)"),
                    newpage=False,
                ),
                10.2,
                4.4,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        name="mat",
                        row_title="row_names_max_width = max_text_width(...)",
                        row_names_max_width=max_text_width(mat2.index, gp=gpar(fontsize=12)),
                    ),
                    newpage=False,
                ),
                10.2,
                4.4,
            ),
        ],
        "ch02_dimnames_max_width.png",
        ncol=1,
    )

    row_lbl_map = {f"row{i}": f"{chr(96 + i)}{i}" for i in range(1, 25)}
    col_lbl_map = {f"column{j}": f"{chr(64 + j)}{j}" for j in range(1, 25)}
    greek_labels = [
        r"$\alpha$", r"$\beta$", r"$\gamma$", r"$\delta$", r"$\epsilon$", r"$\zeta$",
        r"$\eta$", r"$\theta$", r"$\iota$", r"$\kappa$", r"$\lambda$", r"$\mu$",
        r"$\nu$", r"$\xi$", r"$o$", r"$\pi$", r"$\rho$", r"$\sigma$",
    ]
    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_labels=[row_lbl_map[r] for r in mat.index],
                        column_labels=[col_lbl_map[c] for c in mat.columns],
                        column_title="custom row/column labels",
                    ),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_labels=greek_labels,
                        column_title="math labels",
                    ),
                    newpage=False,
                ),
                4.8,
                4.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_labels=gt_render(list("abcdefghijklmnopqr"), padding=unit([2, 10, 2, 10], "pt")),
                        row_names_gp=gpar(box_col=["red", "green"] * 9),
                        column_title="gt_render boxed labels",
                    ),
                    newpage=False,
                ),
                5.0,
                4.5,
            ),
        ],
        "ch02_dimnames_custom_labels.png",
        ncol=3,
    )


def render_2_7_split(mat: pd.DataFrame) -> None:
    def _km_2() -> None:
        ht_opt.kmeans_seed = 123
        draw(Heatmap(mat, name="mat", row_km=2, column_title="row_km = 2"), newpage=False)

    def _km_2_3() -> None:
        ht_opt.kmeans_seed = 123
        draw(Heatmap(mat, name="mat", row_km=2, column_km=3, column_title="row_km = 2, column_km = 3"), newpage=False)

    def _km_no_parent() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(mat, name="mat", row_km=2, column_km=3, show_parent_dend_line=False, column_title="show_parent_dend_line = False"),
            newpage=False,
        )

    _save_grid_panels(
        [
            (_km_2, 4.9, 4.5),
            (_km_2_3, 4.9, 4.5),
            (_km_no_parent, 4.9, 4.5),
        ],
        "ch02_split_kmeans.png",
        ncol=3,
    )

    split_df = pd.DataFrame({"s1": ["A", "B"] * 9, "s2": ["C"] * 9 + ["D"] * 9})
    r_cat = pd.Categorical(["A", "B", "C"] * 6, categories=["C", "B", "A"], ordered=True)
    c_cat = pd.Categorical(["a", "b", "c", "d", "e", "f"] * 4, categories=["f", "e", "d", "c", "b", "a"], ordered=True)
    dend_r = color_branches(as_dendrogram(hclust(dist(mat))), k=2)

    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(mat, name="mat", row_split=["A", "B"] * 9, column_split=["C", "D"] * 12, column_title="vector split"),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", row_split=split_df, column_title="DataFrame multi-column split"),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat,
                        name="mat",
                        row_split=r_cat,
                        column_split=c_cat,
                        cluster_row_slices=False,
                        cluster_column_slices=False,
                        column_title="explicit factor order (no slice cluster)",
                    ),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
            (
                lambda: draw(
                    Heatmap(mat, name="mat", cluster_rows=dend_r, row_split=2, column_title="split by dendrogram (cutree)"),
                    newpage=False,
                ),
                5.2,
                4.6,
            ),
        ],
        "ch02_split_categorical.png",
        ncol=2,
    )

    rng = np.random.default_rng(123)

    def _draw_gp_split() -> None:
        ht_opt.kmeans_seed = 123
        ht_opt.TITLE_PADDING = unit([4, 4], "points")
        draw(
            Heatmap(
                mat,
                name="mat",
                row_km=2,
                row_title_gp=gpar(col=["red", "blue"], fontface=[1, 2]),
                row_names_gp=gpar(col=["green", "orange"], fontsize=[10, 14]),
                column_km=3,
                column_title_gp=gpar(fill=["red", "blue", "green"], fontface=[1, 2, 3]),
                column_names_gp=gpar(col=["green", "orange", "purple"], fontsize=[10, 14, 8]),
            ),
            newpage=False,
        )
        ht_opt(RESET=True)

    def _draw_split_global_title() -> None:
        ht = Heatmap(
            mat,
            name="mat",
            row_split=split_df,
            row_title=["top_slice", "middle_top", "middle_bot", "bottom_slice"],
            row_title_rot=0,
        )
        draw(ht, row_title="Global Row Title", newpage=False)

    def _draw_gaps_border() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(
                mat,
                name="mat",
                row_km=3,
                row_gap=unit([2, 4], "mm"),
                column_km=3,
                column_gap=unit([2, 4], "mm"),
                border=True,
                column_title="custom gaps + border=True",
            ),
            newpage=False,
        )

    def _draw_split_anno_block() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(
                mat,
                name="mat",
                row_km=2,
                column_km=3,
                top_annotation=HeatmapAnnotation(
                    Group=anno_block(gp=gpar(fill=["#E41A1C", "#377EB8", "#4DAF4A"]), labels=["C1", "C2", "C3"], labels_gp=gpar(col="white")),
                    foo1=np.arange(1, 25),
                    bar1=anno_points(rng.uniform(size=24)),
                ),
                right_annotation=rowAnnotation(
                    foo2=np.arange(18, 0, -1),
                    bar2=anno_barplot(rng.uniform(size=18)),
                ),
                border=True,
            ),
            newpage=False,
        )

    _save_grid_panels(
        [
            (_draw_split_global_title, 6.2, 4.9),
            (_draw_gp_split, 6.2, 4.9),
            (_draw_gaps_border, 6.2, 4.9),
            (_draw_split_anno_block, 6.2, 4.9),
        ],
        "ch02_split_titles_gp_gaps.png",
        ncol=2,
    )


def render_2_8_raster() -> None:
    rng = np.random.default_rng(123)
    mat2 = np.empty((5000, 50), dtype=float)
    for i in range(5000):
        half = np.sort(np.abs(rng.standard_normal(50)))[:25]
        mat2[i, :] = rng.uniform(size=50) + np.concatenate([half, half[::-1]]) * ((i + 1) / 1000.0)
    mat2 = mat2[::-1, :]
    brks = np.linspace(np.quantile(mat2, 0.01), np.quantile(mat2, 0.99), 11)
    col_fun2 = colorRamp2(brks, hcl_colors(11, "Spectral", rev=True))

    _save_grid_panels(
        [
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        col=col_fun2,
                        cluster_rows=False,
                        cluster_columns=False,
                        show_heatmap_legend=False,
                        use_raster=True,
                        row_title="use_raster = True",
                        row_title_gp=gpar(fontsize=9),
                    ),
                    newpage=False,
                ),
                2.2,
                2.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        col=col_fun2,
                        cluster_rows=False,
                        cluster_columns=False,
                        show_heatmap_legend=False,
                        use_raster=True,
                        raster_resize_mat=np.nanmean,
                        row_title="resize = mean",
                        row_title_gp=gpar(fontsize=9),
                    ),
                    newpage=False,
                ),
                2.2,
                2.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        col=col_fun2,
                        cluster_rows=False,
                        cluster_columns=False,
                        show_heatmap_legend=False,
                        use_raster=True,
                        raster_resize_mat=np.nanmax,
                        row_title="resize = max",
                        row_title_gp=gpar(fontsize=9),
                    ),
                    newpage=False,
                ),
                2.2,
                2.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        col=col_fun2,
                        cluster_rows=False,
                        cluster_columns=False,
                        show_heatmap_legend=False,
                        use_raster=True,
                        raster_resize_mat=np.nanmin,
                        row_title="resize = min",
                        row_title_gp=gpar(fontsize=9),
                    ),
                    newpage=False,
                ),
                2.2,
                2.5,
            ),
            (
                lambda: draw(
                    Heatmap(
                        mat2,
                        col=col_fun2,
                        cluster_rows=False,
                        cluster_columns=False,
                        show_heatmap_legend=False,
                        use_raster=True,
                        raster_resize_mat=lambda sub: float(sub.flat[0]),
                        row_title="resize = first",
                        row_title_gp=gpar(fontsize=9),
                    ),
                    newpage=False,
                ),
                2.2,
                2.5,
            ),
        ],
        "ch02_raster_compare.png",
        ncol=5,
    )


def render_2_9_customize_body(mat: pd.DataFrame) -> None:
    small_mat = mat.iloc[:9, :9]
    small_arr = small_mat.to_numpy()
    col_fun = colorRamp2([-2, 0, 2], ["green", "white", "red"])

    cor_mat = small_mat.corr()
    od = hclust(dist(cor_mat)).order
    cor_mat = cor_mat.iloc[od, od]
    cor_arr = cor_mat.to_numpy()
    nm = list(cor_mat.index)
    col_cor = colorRamp2([-1, 0, 1], ["green", "white", "red"])

    def _cell_cor(j: int, i: int, x: Any, y: Any, w: Any, h: Any, fill: Any) -> None:
        grid_rect(x=x, y=y, width=w, height=h, gp=gpar(col="grey", fill=None))
        if i == j:
            grid_text(nm[i], x=x, y=y, gp=gpar(fontsize=8.5))
        elif i > j:
            r_mm = min(convertWidth(w, "mm", valueOnly=True), convertHeight(h, "mm", valueOnly=True)) * 0.45 * abs(cor_arr[i, j])
            grid_circle(x=x, y=y, r=unit(r_mm, "mm"), gp=gpar(fill=col_cor(cor_arr[i, j]), col=None))
        else:
            grid_text(f"{cor_arr[i, j]:.1f}", x=x, y=y, gp=gpar(fontsize=8.5))

    def _draw_cell_all() -> None:
        draw(
            Heatmap(
                small_mat,
                name="mat",
                col=col_fun,
                cell_fun=lambda j, i, x, y, w, h, fill: grid_text(f"{small_arr[i, j]:.1f}", x, y, gp=gpar(fontsize=9)),
                column_title="cell_fun: all values",
            ),
            newpage=False,
        )

    def _draw_cell_split() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(
                small_mat,
                name="mat",
                col=col_fun,
                row_km=2,
                column_km=2,
                cell_fun=lambda j, i, x, y, w, h, fill: (
                    grid_text(f"{small_arr[i, j]:.1f}", x, y, gp=gpar(fontsize=9)) if small_arr[i, j] > 0 else None
                ),
                column_title="cell_fun: positive + split",
            ),
            newpage=False,
        )

    def _draw_corrplot() -> None:
        draw(
            Heatmap(
                cor_mat,
                name="correlation",
                col=col_cor,
                rect_gp=gpar(type="none"),
                cell_fun=_cell_cor,
                cluster_rows=False,
                cluster_columns=False,
                show_row_names=False,
                show_column_names=False,
                column_title="corrplot style via cell_fun",
            ),
            newpage=False,
        )

    _save_grid_panels(
        [
            (_draw_cell_all, 4.8, 4.6),
            (_draw_cell_split, 4.8, 4.6),
            (_draw_corrplot, 5.2, 4.6),
        ],
        "ch02_cell_fun_corrplot.png",
        ncol=3,
    )

    # Go game (19x19)
    moves_str = (
        "B[cp];W[pq];B[dc];W[qd];B[eq];W[od];B[de];W[jc];B[qk];W[qn]"
        ";B[qh];W[ck];B[ci];W[cn];B[hc];W[je];B[jq];W[df];B[ee];W[cf]"
        ";B[ei];W[bc];B[ce];W[be];B[bd];W[cd];B[bf];W[ad];B[bg];W[cc]"
        ";B[eb];W[db];B[ec];W[lq];B[nq];W[jp];B[iq];W[kq];B[pp];W[op]"
        ";B[po];W[oq];B[rp];W[ql];B[oo];W[no];B[pl];W[pm];B[np];W[qq]"
        ";B[om];W[ol];B[pk];W[qp];B[on];W[rm];B[mo];W[nr];B[rl];W[rk]"
        ";B[qm];W[dp];B[dq];W[ql];B[or];W[mp];B[nn];W[mq];B[qm];W[bp]"
        ";B[co];W[ql];B[no];W[pr];B[qm];W[dd];B[pn];W[ed];B[bo];W[eg]"
        ";B[ef];W[dg];B[ge];W[gh];B[gf];W[gg];B[ek];W[ig];B[fd];W[en]"
        ";B[bn];W[ip];B[dm];W[ff];B[cb];W[fe];B[hp];W[ho];B[hq];W[el]"
        ";B[dl];W[fk];B[ej];W[fp];B[go];W[hn];B[fo];W[em];B[dn];W[eo]"
        ";B[gp];W[ib];B[gc];W[pg];B[qg];W[ng];B[qc];W[re];B[pf];W[of]"
        ";B[rc];W[ob];B[ph];W[qo];B[rn];W[mi];B[og];W[oe];B[qe];W[rd]"
        ";B[rf];W[pd];B[gm];W[gl];B[fm];W[fl];B[lj];W[mj];B[lk];W[ro]"
        ";B[hl];W[hk];B[ik];W[dk];B[bi];W[di];B[dj];W[dh];B[hj];W[gj]"
        ";B[li];W[lh];B[kh];W[lg];B[jn];W[do];B[cl];W[ij];B[gk];W[bl]"
        ";B[cm];W[hk];B[jk];W[lo];B[hi];W[hm];B[gk];W[bm];B[cn];W[hk]"
        ";B[il];W[cq];B[bq];W[ii];B[sm];W[jo];B[kn];W[fq];B[ep];W[cj]"
        ";B[bk];W[er];B[cr];W[gr];B[gk];W[fj];B[ko];W[kp];B[hr];W[jr]"
        ";B[nh];W[mh];B[mk];W[bb];B[da];W[jh];B[ic];W[id];B[hb];W[jb]"
        ";B[oj];W[fn];B[fs];W[fr];B[gs];W[es];B[hs];W[gn];B[kr];W[is]"
        ";B[dr];W[fi];B[bj];W[hd];B[gd];W[ln];B[lm];W[oi];B[oh];W[ni]"
        ";B[pi];W[ki];B[kj];W[ji];B[so];W[rq];B[if];W[jf];B[hh];W[hf]"
        ";B[he];W[ie];B[hg];W[ba];B[ca];W[sp];B[im];W[sn];B[rm];W[pe]"
        ";B[qf];W[if];B[hk];W[nj];B[nk];W[lr];B[mn];W[af];B[ag];W[ch]"
        ";B[bh];W[lp];B[ia];W[ja];B[ha];W[sf];B[sg];W[se];B[eh];W[fh]"
        ";B[in];W[ih];B[ae];W[so];B[af]"
    )
    coords = [chr(97 + k) for k in range(19)]
    go_mat = pd.DataFrame(np.full((19, 19), None, dtype=object), index=coords, columns=coords)
    for step in moves_str.split(";"):
        m = re.match(r"^(B|W)\[(.)(.)\]$", step.strip())
        if m:
            stone, r_c, c_c = m.group(1), m.group(2), m.group(3)
            if r_c in go_mat.index and c_c in go_mat.columns:
                go_mat.loc[r_c, c_c] = stone
    go_arr = go_mat.to_numpy()
    nr_g, nc_g = go_arr.shape

    def _go_cell(j: int, i: int, x: Any, y: Any, w: Any, h: Any, fill: Any) -> None:
        grid_rect(x, y, w, h, gp=gpar(fill="#dcb35c", col=None))
        if i == 0:
            grid_segments(x, y - h * 0.5, x, y, gp=gpar(col="#222222", lwd=0.8))
        elif i == nr_g - 1:
            grid_segments(x, y, x, y + h * 0.5, gp=gpar(col="#222222", lwd=0.8))
        else:
            grid_segments(x, y - h * 0.5, x, y + h * 0.5, gp=gpar(col="#222222", lwd=0.8))

        if j == 0:
            grid_segments(x, y, x + w * 0.5, y, gp=gpar(col="#222222", lwd=0.8))
        elif j == nc_g - 1:
            grid_segments(x - w * 0.5, y, x, y, gp=gpar(col="#222222", lwd=0.8))
        else:
            grid_segments(x - w * 0.5, y, x + w * 0.5, y, gp=gpar(col="#222222", lwd=0.8))

        if i in (3, 9, 15) and j in (3, 9, 15):
            grid_points(x, y, pch=16, size=unit(2.0, "mm"))

        r_mm = min(convertWidth(w, "mm", valueOnly=True), convertHeight(h, "mm", valueOnly=True)) * 0.43
        val = go_arr[i, j]
        if val == "W":
            grid_circle(x, y, r=unit(r_mm, "mm"), gp=gpar(fill="white", col="#444444", lwd=0.6))
        elif val == "B":
            grid_circle(x, y, r=unit(r_mm, "mm"), gp=gpar(fill="black", col="black"))

    ht_go = Heatmap(
        go_mat,
        name="go",
        col={"B": "black", "W": "white"},
        rect_gp=gpar(type="none"),
        cell_fun=_go_cell,
        cluster_rows=False,
        cluster_columns=False,
        show_row_names=False,
        show_column_names=False,
        column_title="One famous Go game",
        show_heatmap_legend=False,
    )
    lgd_go = Legend(
        title="Players",
        labels=["player1 (B)", "player2 (W)"],
        type="points",
        pch=16,
        legend_gp=gpar(col=["black", "white"]),
        size=unit(4.5, "mm"),
        background="#dcb35c",
    )
    save_ht(ht_go, "ch02_cell_fun_go_game.png", width=6.2, height=5.8, heatmap_legend_list=[lgd_go])
    reset_state()

    # layer_fun showcase
    def _layer_slice_border(j: Any, i: Any, x: Any, y: Any, w: Any, h: Any, fill: Any, slice_r: int, slice_c: int) -> None:
        v = pindex(small_arr, i, j)
        grid_text([f"{val:.1f}" for val in v], x, y, gp=gpar(fontsize=9))
        if slice_r != slice_c:
            grid_rect(gp=gpar(lwd=2.0, col="black", fill="transparent"))

    def _layer_neighbor(j: Any, i: Any, x: Any, y: Any, w: Any, h: Any, fill: Any) -> None:
        ind_mat = restore_matrix(j, i, x, y)
        for ir in range(ind_mat.shape[0]):
            for ic in range(1, ind_mat.shape[1]):
                ind1 = int(ind_mat[ir, ic - 1])
                ind2 = int(ind_mat[ir, ic])
                v1 = small_arr[i[ind1], j[ind1]]
                v2 = small_arr[i[ind2], j[ind2]]
                if v1 * v2 > 0:
                    line_col = "darkred" if v1 > 0 else "darkgreen"
                    grid_segments(x[ind1], y[ind1], x[ind2], y[ind2], gp=gpar(col=line_col, lwd=2.0))
                    grid_points(x[[ind1, ind2]], y[[ind1, ind2]], pch=16, gp=gpar(col=line_col), size=unit(3.2, "mm"))

    def _layer_crosshair(j: Any, i: Any, x: Any, y: Any, w: Any, h: Any, fill: Any) -> None:
        ind_mat = restore_matrix(j, i, x, y)
        inds = np.unique(np.concatenate([ind_mat[1, :], ind_mat[:, 2]]))
        grid_points(x[inds], y[inds], pch=16, size=unit(3.6, "mm"), gp=gpar(col="black"))

    def _draw_l1() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(small_mat, name="mat", col=col_fun, row_km=2, column_km=2, layer_fun=_layer_slice_border, column_title="slice_r != slice_c border"),
            newpage=False,
        )

    def _draw_l2() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(small_mat, name="mat", col=col_fun, row_km=2, column_km=2, layer_fun=_layer_neighbor, column_title="restore_matrix neighbor links"),
            newpage=False,
        )

    def _draw_l3() -> None:
        ht_opt.kmeans_seed = 123
        draw(
            Heatmap(small_mat, name="mat", col=col_fun, row_km=2, column_km=2, layer_fun=_layer_crosshair, column_title="row 1 & col 2 per slice"),
            newpage=False,
        )

    _save_grid_panels(
        [
            (_draw_l1, 4.9, 4.6),
            (_draw_l2, 4.9, 4.6),
            (_draw_l3, 4.9, 4.6),
        ],
        "ch02_layer_fun.png",
        ncol=3,
    )


def render_2_10_size(mat: pd.DataFrame) -> None:
    def _draw_body_8cm() -> None:
        draw(
            Heatmap(
                mat,
                name="mat",
                width=unit(8, "cm"),
                height=unit(8, "cm"),
                column_title="width = 8cm, height = 8cm (body)",
                row_names_gp=gpar(fontsize=9),
                column_names_gp=gpar(fontsize=9),
            ),
            newpage=False,
        )
        grid_rect(gp=gpar(col="#CCCCCC", lty="dashed", fill="transparent"))

    def _draw_total_8cm() -> None:
        ht_opt.TITLE_PADDING = unit([6, 6], "points")
        draw(
            Heatmap(
                mat,
                name="mat",
                heatmap_width=unit(8, "cm"),
                heatmap_height=unit(8, "cm"),
                column_title="heatmap_width = 8cm, heatmap_height = 8cm",
                column_title_gp=gpar(fontsize=9.5),
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8),
            ),
            newpage=False,
        )
        ht_opt(RESET=True)
        grid_rect(gp=gpar(col="#CCCCCC", lty="dashed", fill="transparent"))

    _save_grid_panels(
        [
            (_draw_body_8cm, 5.6, 5.2),
            (_draw_total_8cm, 5.6, 5.2),
        ],
        "ch02_size.png",
        ncol=2,
    )


def render_2_13_subset(mat: pd.DataFrame) -> None:
    rng = np.random.default_rng(123)
    col_fun = colorRamp2([-2, 0, 2], ["green", "white", "red"])
    ht = Heatmap(
        mat,
        name="mat",
        row_km=2,
        column_km=3,
        col=col_fun,
        top_annotation=HeatmapAnnotation(foo1=np.arange(1, 25), bar1=anno_points(rng.uniform(size=24))),
        right_annotation=rowAnnotation(foo2=np.arange(18, 0, -1), bar2=anno_barplot(rng.uniform(size=18))),
    )

    def _sub1() -> None:
        ht_opt.kmeans_seed = 123
        draw(ht[:10, :10], newpage=False)

    def _sub2() -> None:
        ht_opt.kmeans_seed = 123
        draw(ht[::2, 1::2], newpage=False)

    _save_grid_panels(
        [
            (_sub1, 5.4, 4.8),
            (_sub2, 5.4, 4.8),
        ],
        "ch02_subset.png",
        ncol=2,
    )


def main() -> None:
    reset_state()
    mat = _build_main_matrix()
    render_intro(mat)
    render_2_1_colors(mat)
    render_2_2_titles(mat)
    render_2_3_clustering(mat)
    render_2_4_manual_order(mat)
    render_2_5_seriation(mat)
    render_2_6_dimnames(mat)
    render_2_7_split(mat)
    render_2_8_raster()
    render_2_9_customize_body(mat)
    render_2_10_size(mat)
    render_2_13_subset(mat)


if __name__ == "__main__":
    main()
