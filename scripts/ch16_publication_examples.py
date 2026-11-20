"""Generate Chapter 16 figures reproducing examples from the two ComplexHeatmap publications:

1. Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns and
   correlations in multidimensional genomic data. Bioinformatics, 32(18), 2847-2849.
2. Gu, Z. (2022). Complex heatmap visualization. iMeta, 1(3), e43.
"""

from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path
from typing import Any

# Ensure _common is imported first so MPLCONFIGDIR and GridDevice.draw_raster zorder patch apply.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import REPO_ROOT, reset_state, save_Standalone, save_decorated, save_ht

import matplotlib
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy.stats import gaussian_kde, spearmanr

import complexheatmap as ch
from complexheatmap import (
    AnnotationFunction,
    CombMat,
    Heatmap,
    HeatmapAnnotation,
    Legend,
    UpSet,
    Viewport,
    anno_barplot,
    anno_boxplot,
    anno_density,
    anno_histogram,
    anno_horizon,
    anno_image,
    anno_joyplot,
    anno_lines,
    anno_link,
    anno_mark,
    anno_numeric,
    anno_oncoprint_barplot,
    anno_points,
    anno_simple,
    anno_summary,
    anno_text,
    anno_textbox,
    colorRamp2,
    color_branches,
    convertHeight,
    decorate_annotation,
    decorate_column_dend,
    decorate_heatmap_body,
    densityHeatmap,
    frequencyHeatmap,
    gpar,
    grid_boxplot,
    grid_draw,
    grid_grabExpr,
    grid_lines,
    grid_points,
    grid_rect,
    grid_segments,
    grid_text,
    grid_xaxis,
    grid_yaxis,
    grobHeight,
    gt_render,
    oncoPrint,
    popViewport,
    pushViewport,
    restore_matrix,
    rowAnnotation,
    textGrob,
    unit,
    upset_right_annotation,
    upset_top_annotation,
    v,
)

DATA_DIR = REPO_ROOT / "data"


# ==============================================================================
# Publication 1: Gu, Eils & Schlesner (2016), Bioinformatics 32(18), 2847-2849
# ==============================================================================

def render_2016_fig1a_oncoprint_go() -> None:
    """Bioinformatics (2016) Figure 1A / Supplementary S1: TCGA LUAD OncoPrint + Clinical + GO BP."""
    reset_state()
    mut_df = pd.read_csv(DATA_DIR / "pub2016_fig1a_mutsig.csv", index_col=0, keep_default_na=False)
    clin_df = pd.read_csv(DATA_DIR / "pub2016_fig1a_clinical.csv", index_col=0)
    go_df = pd.read_csv(DATA_DIR / "pub2016_fig1a_gobp.csv", index_col=0)

    def _get_type(val: Any) -> list[str]:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return []
        s = str(val).strip()
        if not s:
            return []
        out: list[str] = []
        for part in s.split(";"):
            tp = re.sub(r":.*$", "", part).strip()
            if tp and tp not in out:
                out.append(tp)
        return out

    col_alt = {"MUT": "#008000", "AMP": "red", "HOMDEL": "blue"}
    alter_fun = {
        "background": lambda x, y, w, h: grid_rect(
            x, y, w - unit(0.4, "mm"), h - unit(0.4, "mm"), gp=gpar(fill="#CCCCCC", col=None)
        ),
        "HOMDEL": lambda x, y, w, h: grid_rect(
            x, y, w - unit(0.4, "mm"), h - unit(0.4, "mm"), gp=gpar(fill=col_alt["HOMDEL"], col=None)
        ),
        "AMP": lambda x, y, w, h: grid_rect(
            x, y, w - unit(0.4, "mm"), h - unit(0.4, "mm"), gp=gpar(fill=col_alt["AMP"], col=None)
        ),
        "MUT": lambda x, y, w, h: grid_rect(
            x, y, w - unit(0.4, "mm"), h * 0.35, gp=gpar(fill=col_alt["MUT"], col=None)
        ),
    }

    # Split genes into low vs high amplification rate across the 134 patients
    amp_rates = mut_df.apply(lambda row: np.mean(["AMP" in str(cell) for cell in row]), axis=1)
    amp_labels = np.where(amp_rates > 0.10, "high AMP events", "low AMP events")
    amp_split = pd.Categorical(amp_labels, categories=["low AMP events", "high AMP events"], ordered=True)

    gender = [None if pd.isna(x) else str(x) for x in clin_df["gender"]]
    stage = [None if pd.isna(x) else str(x) for x in clin_df["pathologicstage"]]
    age_vals = pd.to_numeric(clin_df["yearstobirth"], errors="coerce").to_numpy(dtype=float)

    stage_col = {
        "stage i": "#FF0000",
        "stage ia": "#FF6060",
        "stage ib": "#FFB0B0",
        "stage iia": "#60FF60",
        "stage iib": "#B0FFB0",
        "stage iiia": "#6060FF",
        "stage iiib": "#B0B0FF",
        "stage iv": "#FFFF00",
    }
    ha_clin = HeatmapAnnotation(
        gender=gender,
        stage=stage,
        age=anno_points(
            age_vals,
            ylim=(0.0, float(np.nanmax(age_vals))),
            axis=True,
            axis_param={"side": "left", "gp": gpar(fontsize=7.5)},
            size=unit(1.2, "mm"),
            height=unit(14, "mm"),
        ),
        col={"gender": {"male": "red", "female": "blue"}, "stage": stage_col},
        na_col="#E6E6E6",
        annotation_height=unit([4.5, 4.5, 14.0], "mm"),
        show_annotation_name=False,
        annotation_legend_param={"gender": {"title": "Gender"}, "stage": {"title": "Stage"}},
    )

    ht_onco = oncoPrint(
        mut_df,
        get_type=_get_type,
        alter_fun=alter_fun,
        col=col_alt,
        row_split=amp_split,
        column_title="OncoPrint for recurrently mutated genes in Lung Adenocarcinoma",
        column_title_gp=ch.gpar(fontsize=11, fontface="bold"),
        row_title_gp=ch.gpar(fontsize=10),
        row_names_gp=ch.gpar(fontsize=8),
        pct_gp=ch.gpar(fontsize=8),
        heatmap_legend_param={
            "title": "Alterations",
            "at": ["AMP", "HOMDEL", "MUT"],
            "labels": ["Amplification", "Deep deletion", "Mutation"],
        },
        bottom_annotation=ha_clin,
    )

    ha_go_top = HeatmapAnnotation(
        go_count=anno_barplot(
            go_df.sum(axis=0).to_numpy(dtype=float),
            border=False,
            gp=gpar(fill="#6A51A3", col=None),
            height=unit(16, "mm"),
            axis_param={"side": "left", "gp": gpar(fontsize=7.5)},
        ),
        show_annotation_name=False,
    )
    go_labels = [c.lower().replace("_", " ") for c in go_df.columns]
    ht_go = Heatmap(
        go_df.astype(str),
        name="go_bp",
        col={"0": "white", "1": "purple"},
        rect_gp=gpar(col="#B0B0B0", lwd=0.4),
        show_row_names=False,
        cluster_columns=True,
        show_column_dend=False,
        top_annotation=ha_go_top,
        show_column_names=True,
        column_labels=go_labels,
        column_names_side="bottom",
        column_names_rot=45,
        column_names_gp=gpar(fontsize=7.5),
        show_heatmap_legend=False,
        width=unit(5.2, "cm"),
        column_title="Map to Gene Ontology (BP)",
        column_title_gp=gpar(fontsize=10.5, fontface="bold"),
        use_raster=False,
    )

    ha_gene_names = rowAnnotation(
        genes=anno_text(mut_df.index.tolist(), location=unit(0.0, "npc"), just="left", gp=gpar(fontsize=7.5)),
        show_annotation_name=False,
    )

    def _decorate() -> None:
        with decorate_annotation("gender"):
            grid_text("Gender", x=unit(-2, "mm"), just="right", gp=gpar(fontsize=8.5))
        with decorate_annotation("stage"):
            grid_text("Stage", x=unit(-2, "mm"), just="right", gp=gpar(fontsize=8.5))
        with decorate_annotation("age"):
            grid_text("Age", x=unit(-9, "mm"), just="right", gp=gpar(fontsize=8.5))

    save_decorated(
        ht_onco + ha_gene_names + ht_go,
        "ch16_2016_fig1a_oncoprint_go.png",
        _decorate,
        width=14.4,
        height=9.4,
        row_sub_title_side="left",
        ht_gap=unit(3, "mm"),
        merge_legends=True,
        padding=unit([8, 18, 6, 6], "mm"),
    )


def render_2016_fig1b_scrnaseq_721() -> None:
    """Bioinformatics (2016) Figure 1B / Supplementary S2: Full 721-gene mouse T-cell scRNA-seq."""
    reset_state()
    mat_df = pd.read_csv(DATA_DIR / "pub2016_fig1b_scrna721.csv", index_col=0)
    mat_arr = mat_df.to_numpy(dtype=float)

    # Winsorize per gene at [q10, q90] and z-scale across cells
    mat2 = np.empty_like(mat_arr)
    for i in range(mat_arr.shape[0]):
        x = mat_arr[i].copy()
        q10, q90 = np.quantile(x, [0.1, 0.9])
        x = np.clip(x, q10, q90)
        sd = np.std(x, ddof=1)
        mat2[i] = (x - np.mean(x)) / sd if sd > 0 else 0.0
    mat2_df = pd.DataFrame(mat2, index=mat_df.index, columns=mat_df.columns)

    cc_set = set((DATA_DIR / "mouse_cell_cycle_genes.txt").read_text().splitlines())
    rp_set = set((DATA_DIR / "mouse_ribonucleoprotein_genes.txt").read_text().splitlines())
    ccl = np.array([g in cc_set for g in mat_df.index])
    rpl = np.array([g in rp_set for g in mat_df.index])
    base_mean = mat_arr.mean(axis=1)

    q25 = np.quantile(base_mean, 0.25)
    mark_idx = np.where(ccl & (base_mean > q25))[0].tolist()
    mark_labels = [str(mat_df.index[i]) for i in mark_idx]
    cor_mat = np.corrcoef(mat2)

    ht_list = (
        Heatmap(
            mat2_df,
            col=colorRamp2([-1.5, 0, 1.5], ["blue", "white", "red"]),
            name="scaled_expr",
            column_title=f"relative expression for {mat_df.shape[0]} genes",
            column_title_gp=gpar(fontsize=10.5, fontface="bold"),
            show_column_names=False,
            show_row_names=False,
            width=unit(7.6, "cm"),
            heatmap_legend_param={"title": "Scaled expr"},
            use_raster=True,
        )
        + Heatmap(
            base_mean,
            name="base_expr",
            width=unit(4.5, "mm"),
            column_names_gp=gpar(fontsize=9),
            heatmap_legend_param={"title": "Base expr"},
            use_raster=False,
        )
        + Heatmap(
            rpl.astype(int),
            name="ribonucleoprotein",
            col={"0": "white", "1": "purple"},
            show_heatmap_legend=False,
            column_names_gp=gpar(fontsize=9),
            width=unit(4.5, "mm"),
            use_raster=False,
        )
        + Heatmap(
            ccl.astype(int),
            name="cell_cycle",
            col={"0": "white", "1": "red"},
            show_heatmap_legend=False,
            column_names_gp=gpar(fontsize=9),
            width=unit(4.5, "mm"),
            use_raster=False,
        )
        + rowAnnotation(
            link=anno_mark(
                at=mark_idx,
                labels=mark_labels,
                labels_gp=gpar(fontsize=7.5),
                padding=unit(0.8, "mm"),
                link_width=unit(6, "mm"),
            )
        )
        + Heatmap(
            cor_mat,
            name="cor",
            col=colorRamp2([-1, 0, 1], ["green", "white", "red"]),
            show_row_names=False,
            show_column_names=False,
            row_dend_side="right",
            show_column_dend=False,
            column_title="pairwise correlation between genes",
            column_title_gp=gpar(fontsize=10.5, fontface="bold"),
            heatmap_legend_param={"title": "Correlation"},
            use_raster=True,
        )
    )

    def _shade_col_dend(ht_drawn: ch.HeatmapList) -> None:
        tree = ch.column_dend(ht_drawn, name="scaled_expr")
        ind = ch.cutree(tree, k=2)[ch.order_dendrogram(tree)]
        n_cells = len(ind)
        with decorate_column_dend("scaled_expr"):
            for grp_id, fill_col in zip((1, 2), ("#FF000040", "#00FF0040")):
                pos = np.where(ind == grp_id)[0]
                if len(pos) > 0:
                    x_start = float(pos[0]) / n_cells
                    w_frac = float(pos[-1] - pos[0] + 1) / n_cells
                    grid_rect(
                        x=unit(x_start, "npc"),
                        width=unit(w_frac, "npc"),
                        just="left",
                        gp=gpar(fill=fill_col, col=None),
                    )

    save_decorated(
        ht_list,
        "ch16_2016_fig1b_scrnaseq_721.png",
        _shade_col_dend,
        width=13.5,
        height=8.4,
        main_heatmap="cor",
        padding=unit([10, 4, 6, 4], "mm"),
    )


# ==============================================================================
# Publication 2: Gu, Z. (2022), iMeta 1(3), e43
# ==============================================================================

def _draw_single_heatmap_design() -> None:
    th_mm = convertHeight(grobHeight(textGrob("a", gp=gpar(fontsize=10))) * 2.6, "mm", valueOnly=True)
    band = unit(th_mm * 3.8, "mm")

    pushViewport(Viewport(width=unit(1.0, "npc") - unit(4, "mm"), height=unit(1.0, "npc") - unit(4, "mm")))
    try:
        # Top components
        pushViewport(Viewport(x=band, y=unit(1, "npc"), width=unit(1, "npc") - band * 2.0, height=band, just=("left", "top")))
        for i, lbl in enumerate(["column annotations", "column labels", "dendrogram", "title"], start=1):
            pushViewport(Viewport(y=i / 4.0, height=0.25, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()

        # Bottom components
        pushViewport(Viewport(x=band, y=unit(0, "npc"), width=unit(1, "npc") - band * 2.0, height=band, just=("left", "bottom")))
        for i, lbl in enumerate(["title", "dendrogram", "column labels", "column annotations"], start=1):
            pushViewport(Viewport(y=i / 4.0, height=0.25, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()

        # Left components
        pushViewport(Viewport(x=unit(0, "npc"), y=band, width=band, height=unit(1, "npc") - band * 2.0, just=("left", "bottom")))
        for i, lbl in enumerate(["title", "dendrogram", "row labels", "row annotations"], start=1):
            pushViewport(Viewport(x=i / 4.0, width=0.25, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()

        # Right components
        pushViewport(Viewport(x=unit(1, "npc"), y=band, width=band, height=unit(1, "npc") - band * 2.0, just=("right", "bottom")))
        for i, lbl in enumerate(["row annotations", "row labels", "dendrogram", "title"], start=1):
            pushViewport(Viewport(x=i / 4.0, width=0.25, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()

        # Center 2x2 heatmap body slices
        pushViewport(Viewport(x=band, y=band, width=unit(1, "npc") - band * 2.0, height=unit(1, "npc") - band * 2.0, just=("left", "bottom")))
        labels_2x2 = [
            ["matrix, row slice 2\ncolumn slice 1", "matrix, row slice 1\ncolumn slice 1"],
            ["matrix, row slice 2\ncolumn slice 2", "matrix, row slice 1\ncolumn slice 2"],
        ]
        for ci in range(1, 3):
            for rj in range(1, 3):
                pushViewport(Viewport(x=ci / 2.0, y=rj / 2.0, width=0.5, height=0.5, just=("right", "top")))
                pushViewport(Viewport(width=unit(1, "npc") - unit(3, "mm"), height=unit(1, "npc") - unit(3, "mm")))
                grid_rect(gp=gpar(col="red", fill="#FF000020"))
                grid_text(labels_2x2[ci - 1][rj - 1], gp=gpar(fontsize=8.5))
                popViewport()
                popViewport()
        popViewport()
    finally:
        popViewport()


def _draw_annotation_design() -> None:
    pushViewport(Viewport(width=unit(1.0, "npc") - unit(8, "mm"), height=unit(1.0, "npc") - unit(4, "mm")))
    try:
        grid_rect()
        for i, lbl in enumerate(["single annotation 1", "single annotation 2", "single annotation 3"], start=1):
            pushViewport(Viewport(x=i / 3.0, width=1.0 / 3.0, just="right"))
            pushViewport(Viewport(width=unit(1, "npc") - unit(4, "mm"), height=unit(1, "npc") - unit(4, "mm")))
            grid_rect(gp=gpar(col="blue", fill="#0000FF20"))
            grid_text(lbl, rot=90, gp=gpar(fontsize=9))
            popViewport()
            popViewport()
    finally:
        popViewport()


def _draw_list_design() -> None:
    th_mm = convertHeight(grobHeight(textGrob("a", gp=gpar(fontsize=10))) * 2.6, "mm", valueOnly=True)
    band = unit(th_mm * 1.9, "mm")

    pushViewport(Viewport(width=unit(1.0, "npc") - unit(4, "mm"), height=unit(1.0, "npc") - unit(4, "mm")))
    try:
        # Top
        pushViewport(Viewport(x=band, y=unit(1, "npc"), width=unit(1, "npc") - band * 2.0, height=band, just=("left", "top")))
        for i, lbl in enumerate(["title", "legends"], start=1):
            pushViewport(Viewport(y=i / 2.0, height=0.5, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()
        # Bottom
        pushViewport(Viewport(x=band, y=unit(0, "npc"), width=unit(1, "npc") - band * 2.0, height=band, just=("left", "bottom")))
        for i, lbl in enumerate(["legends", "title"], start=1):
            pushViewport(Viewport(y=i / 2.0, height=0.5, just="top"))
            grid_rect()
            grid_text(lbl, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()
        # Left
        pushViewport(Viewport(x=unit(0, "npc"), y=band, width=band, height=unit(1, "npc") - band * 2.0, just=("left", "bottom")))
        for i, lbl in enumerate(["legends", "title"], start=1):
            pushViewport(Viewport(x=i / 2.0, width=0.5, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()
        # Right
        pushViewport(Viewport(x=unit(1, "npc"), y=band, width=band, height=unit(1, "npc") - band * 2.0, just=("right", "bottom")))
        for i, lbl in enumerate(["title", "legends"], start=1):
            pushViewport(Viewport(x=i / 2.0, width=0.5, just="right"))
            grid_rect()
            grid_text(lbl, rot=90, gp=gpar(fontsize=8.5))
            popViewport()
        popViewport()
        # Center list items
        pushViewport(Viewport(x=band, y=band, width=unit(1, "npc") - band * 2.0, height=unit(1, "npc") - band * 2.0, just=("left", "bottom")))
        for i, lbl in enumerate(["heatmap 1", "heatmap 2", "row annotation"], start=1):
            pushViewport(Viewport(x=i / 3.0, width=1.0 / 3.0, just="right"))
            pushViewport(Viewport(width=unit(1, "npc") - unit(3, "mm"), height=unit(1, "npc") - unit(3, "mm")))
            grid_rect(gp=gpar(col="green", fill="#00FF0020"))
            grid_text(lbl, rot=90, gp=gpar(fontsize=8.5))
            popViewport()
            popViewport()
        popViewport()
    finally:
        popViewport()


def render_2022_fig1_design() -> None:
    """iMeta (2022) Figure 1: Three-panel modular design of ComplexHeatmap."""
    reset_state()

    def _draw() -> None:
        # Column proportions matching figure1.R widths = (2.5, 1, 1.2)
        total = 2.5 + 1.0 + 1.2
        w1, w2, w3 = 2.5 / total, 1.0 / total, 1.2 / total
        x1, x2, x3 = 0.0, w1, w1 + w2
        hdr_h = 0.09
        body_h = 1.0 - hdr_h

        titles = [
            (x1, w1, "A) A single heatmap", _draw_single_heatmap_design),
            (x2, w2, "B) Heatmap annotation", _draw_annotation_design),
            (x3, w3, "C) Heatmap list", _draw_list_design),
        ]
        for x0, wk, title_txt, fn in titles:
            pushViewport(Viewport(x=x0, y=1.0, width=wk, height=hdr_h, just=("left", "top")))
            grid_text(title_txt, gp=gpar(fontsize=11.5, fontface="bold"))
            popViewport()

            pushViewport(Viewport(x=x0, y=0.0, width=wk, height=body_h, just=("left", "bottom")))
            fn()
            popViewport()

    save_Standalone(_draw, "ch16_2022_fig1_design.png", width=11.8, height=4.8, dpi=150)


def _decorate_dend_nodes(dend: ch.Dendrogram, rng: np.random.Generator) -> None:
    """Recursively attach random colored glyphs onto internal dendrogram nodes."""
    if rng.random() > 0.45:
        pchs = [1, 2, 15, 16, 17, 18]
        cols = ch.brewer_pal(6, "Dark2")
        c = cols[int(rng.integers(0, len(cols)))]
        dend.node_par = {
            "pch": int(pchs[int(rng.integers(0, len(pchs)))]),
            "col": c,
            "fill": c,
            "size": unit(2.2, "mm"),
        }
    for ch_node in dend.children:
        _decorate_dend_nodes(ch_node, rng)


def render_2022_fig2ab_single() -> None:
    """iMeta (2022) Figure 2A & 2B: Multi-annotation split heatmap + layer_fun customization."""
    reset_state()
    rng = np.random.default_rng(123)
    nr1, nr2, nr3 = 4, 8, 6
    nc1, nc2, nc3 = 6, 8, 10
    nr = nr1 + nr2 + nr3
    nc = nc1 + nc2 + nc3

    b1 = np.vstack([rng.normal(1.0, 0.5, (nr1, nc1)), rng.normal(0.0, 0.5, (nr2, nc1)), rng.normal(0.0, 0.5, (nr3, nc1))])
    b2 = np.vstack([rng.normal(0.0, 0.5, (nr1, nc2)), rng.normal(1.0, 0.5, (nr2, nc2)), rng.normal(0.0, 0.5, (nr3, nc2))])
    b3 = np.vstack([rng.normal(0.5, 0.5, (nr1, nc3)), rng.normal(0.5, 0.5, (nr2, nc3)), rng.normal(1.0, 0.5, (nr3, nc3))])
    mat = np.hstack([b1, b2, b3])
    mat = mat[rng.permutation(nr), :][:, rng.permutation(nc)]
    mat_df = pd.DataFrame(mat, index=[f"row{i}" for i in range(1, nr + 1)], columns=[f"column{j}" for j in range(1, nc + 1)])

    col_dend = ch.as_dendrogram(ch.hclust(ch.dist(mat_df.to_numpy().T)))
    col_dend = color_branches(col_dend, k=3)
    _decorate_dend_nodes(col_dend, np.random.default_rng(77))

    row_ab = ["A", "B"] * 9
    ha_top = HeatmapAnnotation(
        foo1=np.arange(1, nc + 1),
        bar1=anno_points(rng.uniform(size=nc), size=unit(1.5, "mm"), axis_param={"gp": gpar(fontsize=7)}),
        annotation_name_gp=gpar(fontsize=8.5),
    )
    bar2_mat = np.column_stack([rng.uniform(size=nr), rng.uniform(size=nr)])
    ha_right = rowAnnotation(
        foo2=np.arange(nr, 0, -1),
        bar2=anno_barplot(
            bar2_mat,
            gp=gpar(fill=["#DF536B", "#61D04F"], col="white"),
            width=unit(2.0, "cm"),
            axis_param={"gp": gpar(fontsize=7)},
        ),
        annotation_name_gp=gpar(fontsize=8.5),
    )
    ht_a = Heatmap(
        mat_df,
        name="mat",
        row_dend_width=unit(1.8, "cm"),
        cluster_columns=col_dend,
        column_title="A) A heatmap with various annotations",
        column_title_gp=gpar(fontsize=11, fontface="bold"),
        show_column_names=False,
        row_split=row_ab,
        row_km=2,
        column_split=3,
        row_title_gp=gpar(fontsize=9),
        row_names_gp=gpar(fontsize=8),
        top_annotation=ha_top,
        right_annotation=ha_right,
    )
    lgd_bar2 = Legend(title="bar2", labels=["group1", "group2"], legend_gp=gpar(fill=["#DF536B", "#61D04F"]))

    def _draw_a() -> None:
        ch.draw(ht_a, newpage=False, annotation_legend_list=[lgd_bar2], merge_legends=True, padding=unit([6, 4, 4, 4], "mm"))
        for sl in (1, 2, 3):
            with decorate_annotation("bar1", slice=sl):
                grid_lines([0.0, 1.0], [0.5, 0.5], gp=gpar(lty="dashed", col="#888888"), default_units="npc")

    # Panel B: 9x9 split heatmap with layer_fun overlays
    small_mat = mat[:9, :9].copy()
    small_df = pd.DataFrame(small_mat, index=[f"row{i}" for i in range(1, 10)], columns=[f"column{j}" for j in range(1, 10)])
    col_b = colorRamp2([-2, 0, 2], ["green", "white", "red"])

    def _layer_b(j: Any, i: Any, x: Any, y: Any, w: Any, h: Any, fill: Any, kr: int, kc: int) -> None:
        ind_mat = restore_matrix(j, i, x, y)
        nr_s, nc_s = ind_mat.shape
        for ir in range(nr_s):
            for ic in range(1, nc_s):
                p1 = int(ind_mat[ir, ic - 1])
                p2 = int(ind_mat[ir, ic])
                v1 = float(small_mat[int(i[p1]), int(j[p1])])
                v2 = float(small_mat[int(i[p2]), int(j[p2])])
                if v1 * v2 > 0:
                    c_line = "violet" if v1 > 0 else "blue"
                    grid_segments(x[p1], y[p1], x[p2], y[p2], gp=gpar(col=c_line, lwd=1.8))
                    grid_points(x[[p1, p2]], y[[p1, p2]], pch=16, size=unit(2.6, "mm"), gp=gpar(col=c_line))
        if kr != kc:
            grid_rect(gp=gpar(lwd=2.2, fill="transparent", col="black"))

    ht_b = Heatmap(
        small_df,
        name="mat_b",
        col=col_b,
        column_title="B) Customize heatmap body",
        column_title_gp=gpar(fontsize=11, fontface="bold"),
        row_km=2,
        column_km=2,
        layer_fun=_layer_b,
        row_names_gp=gpar(fontsize=8.5),
        column_names_gp=gpar(fontsize=8.5),
        heatmap_legend_param={"title": "mat"},
    )

    gb_a = grid_grabExpr(_draw_a, width_in=6.8, height_in=5.4)
    gb_b = grid_grabExpr(lambda: ch.draw(ht_b, newpage=False, padding=unit([6, 4, 4, 4], "mm")), width_in=5.6, height_in=5.4)

    def _draw_ab() -> None:
        pushViewport(Viewport(x=0.0, width=0.55, just="left"))
        grid_draw(gb_a)
        popViewport()
        pushViewport(Viewport(x=0.55, width=0.45, just="left"))
        grid_draw(gb_b)
        popViewport()

    save_Standalone(_draw_ab, "ch16_2022_fig2ab_single.png", width=12.4, height=5.5, dpi=140)


def render_2022_fig2cd_presplit() -> None:
    """iMeta (2022) Figure 2C & 2D: TCGA GBM global clustering vs pre-split by consensus class."""
    reset_state()
    expr_df = pd.read_csv(DATA_DIR / "pub2022_fig2cd_expr.csv", index_col=0)
    meta_df = pd.read_csv(DATA_DIR / "pub2022_fig2cd_meta.csv")

    cls_vec = [str(c) for c in meta_df["class"]]
    sub_vec = [str(s) for s in meta_df["subtype"]]
    cls_col = {"1": "#66C2A5", "2": "#FC8D62", "3": "#8DA0CB", "4": "#E78AC3"}
    sub_col = {"Proneural": "#1B9E77", "Neural": "#D95F02", "Classical": "#7570B3", "Mesenchymal": "#E7298A"}

    col_fun = colorRamp2([-2, 0, 2], ["green", "white", "red"])
    ha_top = HeatmapAnnotation(
        **{"class": cls_vec},
        col={"class": cls_col},
        annotation_name_gp=gpar(fontsize=9),
    )
    ha_bot = HeatmapAnnotation(
        subtype=sub_vec,
        col={"subtype": sub_col},
        annotation_name_gp=gpar(fontsize=9),
    )

    ht_c = Heatmap(
        expr_df,
        name="mat",
        col=col_fun,
        show_row_names=False,
        show_row_dend=False,
        show_column_names=False,
        column_title="C) Cluster all columns together",
        column_title_gp=gpar(fontsize=11, fontface="bold"),
        top_annotation=ha_top,
        bottom_annotation=ha_bot,
        use_raster=True,
    )

    ha_top_d = HeatmapAnnotation(
        **{"class": cls_vec},
        col={"class": cls_col},
        annotation_name_gp=gpar(fontsize=9),
    )
    ha_bot_d = HeatmapAnnotation(
        subtype=sub_vec,
        col={"subtype": sub_col},
        annotation_name_gp=gpar(fontsize=9),
    )
    ht_d = Heatmap(
        expr_df,
        name="mat",
        col=col_fun,
        row_km=4,
        column_split=pd.Categorical(cls_vec, categories=["1", "2", "3", "4"], ordered=True),
        show_row_names=False,
        show_row_dend=False,
        show_column_names=False,
        row_title_gp=gpar(fontsize=9),
        column_title="D) Columns are pre-split by 'class'",
        column_title_gp=gpar(fontsize=11, fontface="bold"),
        top_annotation=ha_top_d,
        bottom_annotation=ha_bot_d,
        use_raster=True,
    )

    gb_c = grid_grabExpr(lambda: ch.draw(ht_c, newpage=False, merge_legends=True, padding=unit([4, 4, 4, 4], "mm")), width_in=5.8, height_in=5.6)
    gb_d = grid_grabExpr(lambda: ch.draw(ht_d, newpage=False, merge_legends=True, padding=unit([4, 4, 4, 4], "mm")), width_in=6.2, height_in=5.6)

    def _draw_cd() -> None:
        pushViewport(Viewport(x=0.0, width=0.49, just="left"))
        grid_draw(gb_c)
        popViewport()
        pushViewport(Viewport(x=0.49, width=0.51, just="left"))
        grid_draw(gb_d)
        popViewport()

    save_Standalone(_draw_cd, "ch16_2022_fig2cd_presplit.png", width=12.0, height=5.7, dpi=140)


def _create_gallery_icons(out_dir: Path, n: int = 7) -> list[str]:
    """Generate clean dark silhouette icons for the Figure 3A anno_image track."""
    size = 96
    pad = 14
    mid = size // 2
    ink = "#2C3E50"
    paths: list[str] = []
    for i in range(n):
        im = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        dr = ImageDraw.Draw(im)
        kind = i % 7
        if kind == 0:
            dr.ellipse([pad, pad, size - pad, size - pad], fill=ink)
        elif kind == 1:
            dr.polygon([(mid, pad), (size - pad, mid), (mid, size - pad), (pad, mid)], fill=ink)
        elif kind == 2:
            dr.rounded_rectangle([pad, pad, size - pad, size - pad], radius=16, fill=ink)
        elif kind == 3:
            dr.polygon([(mid, pad), (size - pad, size - pad), (pad, size - pad)], fill=ink)
        elif kind == 4:
            dr.rectangle([mid - 11, pad, mid + 11, size - pad], fill=ink)
            dr.rectangle([pad, mid - 11, size - pad, mid + 11], fill=ink)
        elif kind == 5:
            dr.rectangle([18, 52, 34, 80], fill=ink)
            dr.rectangle([40, 26, 56, 80], fill=ink)
            dr.rectangle([62, 38, 78, 80], fill=ink)
        else:
            dr.ellipse([pad, pad, size - pad, size - pad], outline=ink, width=10)
            dr.line([(28, 52), (44, 68), (70, 32)], fill=ink, width=10)
        p = out_dir / f"gallery_icon_{i}.png"
        im.save(p)
        paths.append(str(p))
    return paths


def render_2022_fig3a_gallery() -> None:
    """iMeta (2022) Figure 3A: 13-track gallery of built-in annotation graphics."""
    reset_state()
    rng = np.random.default_rng(123)
    n = 7
    m = rng.standard_normal((n, 100))
    m2 = rng.standard_normal((n, 100))
    pct_vals = rng.uniform(0.18, 0.96, n)

    m_joy = rng.standard_normal((100, 15))
    gx = np.linspace(-4, 4, 150)
    lt1 = [pd.DataFrame({"x": gx, "y": gaussian_kde(m_joy[:, j])(gx)}) for j in range(m_joy.shape[1])]
    lt2 = [np.cumprod(1.0 + rng.uniform(-k / 100.0, k / 100.0, 400)) - 1.0 for k in range(1, 16)]

    rich_txt = gt_render(
        [f"row$_{k}$" for k in range(1, n + 1)],
        align_widths=True,
        halign=0.5,
        r=unit(3, "pt"),
        padding=unit([2, 6, 2, 6], "pt"),
    )

    with tempfile.TemporaryDirectory() as td:
        icons = _create_gallery_icons(Path(td), n=n)
        ha1 = rowAnnotation(
            anno_simple=anno_simple(rng.uniform(size=n), width=unit(8, "mm")),
            anno_image=anno_image(icons, gp=gpar(col="black", fill="#F7F7F7"), space=unit(1.5, "mm"), width=unit(9, "mm")),
            anno_points=anno_points(rng.uniform(size=(n, 2)), pch=[16, 17], gp=gpar(col=["#DF536B", "#61D04F"]), size=unit(1.5, "mm"), width=unit(1.5, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_lines=anno_lines(rng.uniform(size=(n, 2)), gp=gpar(col=["#DF536B", "#61D04F"], lty=["solid", "dashed"]), width=unit(1.5, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_loess=anno_lines(rng.uniform(size=n), smooth=True, size=unit(1.5, "mm"), width=unit(1.5, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_barplot=anno_barplot(rng.uniform(size=(n, 2)), gp=gpar(fill=["#DF536B", "#61D04F"], col="white"), width=unit(1.6, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_pct=anno_numeric(pct_vals, labels_format=lambda v: [f"{val * 100:.2f}%" for val in v], bg_gp=gpar(fill="#8080FF", col="#8080FF"), labels_gp=gpar(fontsize=8), align_to="right", width=unit(1.8, "cm")),
            anno_boxplot=anno_boxplot(m, gp=gpar(fill=list(range(1, n + 1))), width=unit(2.0, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_text=anno_text(rich_txt, gp=gpar(box_col="blue", box_lwd=1.5, col=list(range(1, n + 1)), fontsize=10), just="center", location=unit(0.5, "npc"), width=unit(1.5, "cm")),
            anno_histogram=anno_histogram(m2, gp=gpar(fill=list(range(1, n + 1))), border=True, width=unit(2.0, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            anno_violin=anno_density(m, type="violin", gp=gpar(fill=list(range(1, n + 1))), width=unit(2.0, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            show_legend=False,
            gap=unit(1.8, "mm"),
            show_annotation_name=False,
        )
        ha3 = rowAnnotation(
            anno_joyplot=anno_joyplot(lt1, width=unit(2.4, "cm"), gp=gpar(fill=list(range(1, 16))), transparency=0.75, scale=1.8, axis_param={"gp": gpar(fontsize=7)}),
            anno_horizon=anno_horizon(lt2, width=unit(2.2, "cm"), axis_param={"gp": gpar(fontsize=7)}),
            gap=unit(2.0, "mm"),
            show_annotation_name=False,
        )

        w1_in = ch.width(ha1).to_mm_scalar() / 25.4
        w3_in = ch.width(ha3).to_mm_scalar() / 25.4
        gap_in = 0.15
        tot_w = max(13.8, w1_in + w3_in + gap_in + 0.5)

        def _draw_gallery() -> None:
            grid_text("A) Different annotation graphics supported in ComplexHeatmap-py", x=0.5, y=unit(1, "npc") - unit(2.5, "mm"), just="top", gp=gpar(fontsize=11.5, fontface="bold"))

            frac1 = w1_in / (w1_in + w3_in + gap_in)
            frac3 = w3_in / (w1_in + w3_in + gap_in)

            def _d1() -> None:
                ha1.draw(y=unit(0.45, "npc"), height=unit(0.72, "npc"), width=ch.width(ha1))
                labels_1 = {
                    "anno_simple": "Simple",
                    "anno_image": "Images",
                    "anno_points": "Points",
                    "anno_lines": "Lines",
                    "anno_loess": "Smoothed\nlines",
                    "anno_barplot": "Barplot",
                    "anno_pct": "Percent\nwith bars",
                    "anno_boxplot": "Boxplot",
                    "anno_text": "Customized\ntext",
                    "anno_histogram": "Histogram",
                    "anno_violin": "Violin plot",
                }
                for an, lbl in labels_1.items():
                    with decorate_annotation(an):
                        grid_text(lbl, y=unit(1, "npc") + unit(3.0, "mm"), just="bottom", gp=gpar(fontsize=8.5))

            def _d3() -> None:
                ha3.draw(y=unit(0.45, "npc"), height=unit(0.72, "npc"), width=ch.width(ha3))
                for an, lbl in {"anno_joyplot": "Joy plot", "anno_horizon": "Horizon chart"}.items():
                    with decorate_annotation(an):
                        grid_text(lbl, y=unit(1, "npc") + unit(5.5, "mm"), just="bottom", gp=gpar(fontsize=8.5))

            gb1 = grid_grabExpr(_d1, width_in=w1_in + 0.3, height_in=4.5)
            gb3 = grid_grabExpr(_d3, width_in=w3_in + 0.3, height_in=4.5)

            pushViewport(Viewport(x=0.01, y=0.02, width=0.98, height=0.88, just=("left", "bottom")))
            pushViewport(Viewport(x=0.0, width=frac1, just="left"))
            grid_draw(gb1)
            popViewport()
            pushViewport(Viewport(x=1.0 - frac3, width=frac3, just="left"))
            grid_draw(gb3)
            popViewport()
            popViewport()

        save_Standalone(_draw_gallery, "ch16_2022_fig3a_gallery.png", width=tot_w, height=5.2, dpi=140)


def render_2022_fig3bcd_complex_annos() -> None:
    """iMeta (2022) Figure 3B, 3C, 3D: anno_mark, anno_link boxplots, anno_textbox word clouds."""
    reset_state()
    rng = np.random.default_rng(123)
    nr1, nr2, nr3 = 40, 80, 60
    nc1, nc2, nc3 = 6, 8, 10
    nr = nr1 + nr2 + nr3

    b1 = np.vstack([rng.normal(1.0, 0.5, (nr1, nc1)), rng.normal(0.0, 0.5, (nr2, nc1)), rng.normal(0.0, 0.5, (nr3, nc1))])
    b2 = np.vstack([rng.normal(0.0, 0.5, (nr1, nc2)), rng.normal(1.0, 0.5, (nr2, nc2)), rng.normal(0.0, 0.5, (nr3, nc2))])
    b3 = np.vstack([rng.normal(0.5, 0.5, (nr1, nc3)), rng.normal(0.5, 0.5, (nr2, nc3)), rng.normal(1.0, 0.5, (nr3, nc3))])
    mat = np.hstack([b1, b2, b3])[rng.permutation(nr), :]

    # Panel B: Mark annotation
    genes = sorted(rng.choice(nr, size=20, replace=False).tolist())
    labels_b = [f"gene{k}" for k in range(1, 21)]
    ht_b = Heatmap(
        mat,
        name="mat_b",
        column_title="B) Mark annotation",
        column_title_gp=gpar(fontsize=11, fontface="bold"),
        right_annotation=rowAnnotation(mark=anno_mark(at=genes, labels=labels_b, labels_gp=gpar(fontsize=8))),
        show_row_names=False,
        heatmap_legend_param={"title": "mat"},
        use_raster=False,
    )

    # Panel C: Link annotation connecting 3 row slices to grouped boxplots
    from scipy.cluster.vq import kmeans2
    _, km_raw = kmeans2(mat, 3, minit="points", seed=123)
    km = [str(k + 1) for k in km_raw]
    ind_list = {str(k + 1): np.where(km_raw == k)[0].tolist() for k in range(3)}
    groups = ["group1"] * 6 + ["group2"] * 8 + ["group3"] * 10
    g_pal = ["#66C2A5", "#FC8D62", "#8DA0CB"]
    yrg = (float(mat[:, :6].mean(1).min()) - 0.5, float(mat.max()) * 0.45 + 0.8)

    def _panel_boxplot(index: list[int], nm: str) -> None:
        mm = mat[index, :]
        vals = [mm[:, :6].mean(axis=1), mm[:, 6:14].mean(axis=1), mm[:, 14:24].mean(axis=1)]
        all_v = np.concatenate(vals)
        ylo, yhi = float(all_v.min()) - 0.15, float(all_v.max()) + 0.15
        pushViewport(Viewport(xscale=(0.4, 3.6), yscale=(ylo, yhi), clip=False))
        try:
            grid_rect(gp=gpar(fill="#F7F7F7", col="black"))
            for gi, arr_g in enumerate(vals, start=1):
                grid_boxplot(arr_g, pos=float(gi), box_width=0.5, gp=gpar(fill=g_pal[gi - 1], col="#333333"), direction="vertical")
            grid_xaxis(at=[1.0, 2.0, 3.0], label=["group1", "group2", "group3"], gp=gpar(fontsize=7.5))
            grid_yaxis(gp=gpar(fontsize=7))
        finally:
            popViewport()

    ht_c = (
        Heatmap(
            mat,
            name="mat_c",
            row_split=pd.Categorical(km, categories=["1", "2", "3"], ordered=True),
            column_split=pd.Categorical(groups, categories=["group1", "group2", "group3"], ordered=True),
            cluster_column_slices=False,
            column_title=["group1", "group2", "group3"],
            column_title_gp=gpar(fontsize=9),
            show_row_names=False,
            show_column_names=False,
            heatmap_legend_param={"title": "mat"},
            use_raster=False,
        )
        + rowAnnotation(
            link=anno_link(
                align_to=ind_list,
                which="row",
                panel_fun=_panel_boxplot,
                size=unit(2.1, "cm"),
                width=unit(4.6, "cm"),
                gap=unit(4, "mm"),
                internal_line=False,
            )
        )
    )

    # Panel D: Textbox word clouds across 8 row slices
    mat_d = rng.standard_normal((100, 10))
    split_d = rng.choice(list("abcdefgh"), size=100, replace=True)
    levels_d = sorted(set(split_d))
    wc_pal = ch.brewer_pal(8, "Dark2")
    text_d = {}
    for lvl in levels_d:
        words = [f"function{k}" for k in range(1, 9)]
        cols = [wc_pal[int(rng.integers(0, len(wc_pal)))] for _ in range(8)]
        fsize = rng.uniform(5.5, 10.5, 8)
        text_d[lvl] = pd.DataFrame({"word": words, "col": cols, "fontsize": fsize})

    ht_d = (
        Heatmap(
            mat_d,
            name="mat_d",
            cluster_rows=False,
            row_split=pd.Categorical(split_d, categories=levels_d, ordered=True),
            show_row_names=False,
            column_title="D) Textbox annotation",
            column_title_gp=gpar(fontsize=11, fontface="bold"),
            heatmap_legend_param={"title": "mat"},
            use_raster=False,
        )
        + rowAnnotation(wc=anno_textbox(split_d, text_d, max_width=unit(6.6, "cm")))
    )

    gb_b = grid_grabExpr(lambda: ch.draw(ht_b, newpage=False, padding=unit([6, 4, 4, 4], "mm")), width_in=3.8, height_in=7.0)
    gb_c = grid_grabExpr(lambda: ch.draw(ht_c, newpage=False, column_title="C) Link annotation", column_title_gp=gpar(fontsize=11, fontface="bold"), padding=unit([8, 4, 4, 4], "mm")), width_in=5.6, height_in=7.0)
    gb_d = grid_grabExpr(lambda: ch.draw(ht_d, newpage=False, padding=unit([6, 4, 4, 4], "mm")), width_in=5.4, height_in=7.0)

    def _draw_bcd() -> None:
        w_tot = 3.8 + 5.6 + 5.4
        f1, f2, f3 = 3.8 / w_tot, 5.6 / w_tot, 5.4 / w_tot
        pushViewport(Viewport(x=0.0, width=f1, just="left"))
        grid_draw(gb_b)
        popViewport()
        pushViewport(Viewport(x=f1, width=f2, just="left"))
        grid_draw(gb_c)
        popViewport()
        pushViewport(Viewport(x=f1 + f2, width=f3, just="left"))
        grid_draw(gb_d)
        popViewport()

    save_Standalone(_draw_bcd, "ch16_2022_fig3bcd_complex_annos.png", width=14.8, height=7.1, dpi=140)


def render_2022_fig4a_meth_list() -> None:
    """iMeta (2022) Figure 4A: DNA methylation, expression, TSS distance, and enhancers."""
    reset_state()
    ds = ch.load_dataset("meth")
    mat_meth = ds["meth"].to_numpy(dtype=float)
    mat_expr = ds["expr"].to_numpy(dtype=float)
    anno_enhancer = ds["enhancer"].to_numpy(dtype=float)
    rm = ds["row_meta"]

    direction = rm["direction"].tolist()
    cor_pvalue = rm["cor_pvalue"].to_numpy(dtype=float)
    gene_type = rm["gene_type"].tolist()
    anno_gene = ["promoter" if str(x) == "TSS" else str(x) for x in rm["anno"]]

    rng = np.random.default_rng(42)
    dist_tss = np.array([
        rng.uniform(0, 1000) if row.mean() < 0.2 else (rng.uniform(1000, 5000) if row.mean() < 0.4 else rng.uniform(1000, 50000))
        for row in mat_meth
    ])

    col_names = [f"sample{i}" for i in range(1, 21)]
    sample_type = ["Tumor"] * 10 + ["Control"] * 10
    col_order = ch.hclust(ch.dist(mat_meth.T)).order

    gt_levels = list(dict.fromkeys(gene_type))
    gene_type_col = dict(zip(gt_levels, ch.brewer_pal(len(gt_levels), "Set3")))
    ag_levels = list(dict.fromkeys(anno_gene))
    anno_gene_col = dict(zip(ag_levels, ch.brewer_pal(len(ag_levels), "Set1")))

    ha1 = HeatmapAnnotation(
        Type=sample_type,
        col={"Type": {"Tumor": "pink", "Control": "royalblue"}},
        annotation_name_side="left",
        annotation_legend_param={"Type": {"legend_direction": "horizontal"}},
    )
    ha2 = HeatmapAnnotation(
        **{"type": sample_type},
        col={"type": {"Tumor": "pink", "Control": "royalblue"}},
        show_legend=False,
        show_annotation_name=False,
    )

    ht_list = (
        Heatmap(
            pd.DataFrame(mat_meth, columns=col_names),
            name="Methylation",
            col=colorRamp2([0, 0.5, 1], ["blue", "white", "red"]),
            column_order=col_order,
            top_annotation=ha1,
            column_title="Methylation",
            show_row_names=False,
            column_names_gp=gpar(fontsize=8),
            heatmap_legend_param={"legend_direction": "horizontal", "legend_width": unit(3, "cm")},
            use_raster=True,
        )
        + Heatmap(
            direction,
            name="Direction",
            col={"hyper": "red", "hypo": "blue"},
            show_row_names=False,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"legend_direction": "horizontal"},
        )
        + Heatmap(
            pd.DataFrame(mat_expr, columns=col_names),
            name="Expression",
            col=colorRamp2([-2, 0, 2], ["green", "white", "red"]),
            column_order=col_order,
            top_annotation=ha2,
            column_title="Expression",
            show_row_names=False,
            column_names_gp=gpar(fontsize=8),
            heatmap_legend_param={"legend_direction": "horizontal", "legend_width": unit(3, "cm")},
            use_raster=True,
        )
        + Heatmap(
            cor_pvalue,
            name="p-value",
            col=colorRamp2([0, 2, 4], ["white", "white", "red"]),
            show_row_names=False,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"legend_direction": "horizontal", "at": [0, 2, 4], "labels": ["1", "0.01", "<1e-4"]},
        )
        + Heatmap(
            gene_type,
            name="Gene type",
            col=gene_type_col,
            show_row_names=False,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"legend_direction": "horizontal", "nrow": 2},
        )
        + Heatmap(
            anno_gene,
            name="Gene annotation",
            col=anno_gene_col,
            show_row_names=False,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"legend_direction": "horizontal", "nrow": 2},
        )
        + rowAnnotation(
            **{
                "Distance to TSS": anno_points(
                    dist_tss,
                    size=unit(1.1, "mm"),
                    width=unit(3.6, "cm"),
                    axis_param={"at": [0, 20000, 40000], "labels": ["0kb", "20kb", "40kb"], "gp": gpar(fontsize=8)},
                )
            }
        )
        + Heatmap(
            pd.DataFrame(anno_enhancer, columns=["enhancer_1", "enhancer_2", "enhancer_3"]),
            name="Enhancer",
            col=colorRamp2([0, 1], ["white", "orange"]),
            cluster_columns=False,
            column_title="Enhancer",
            show_row_names=False,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"legend_direction": "horizontal", "legend_width": unit(3, "cm")},
            use_raster=False,
        )
    )

    save_ht(
        ht_list,
        "ch16_2022_fig4a_meth_list.png",
        width=13.5,
        height=8.4,
        row_km=2,
        row_split=direction,
        column_title="A) Visualize the associations between methylation, expression and other genomic features",
        column_title_gp=gpar(fontsize=11.5, fontface="bold"),
        padding=unit([4, 6, 8, 4], "mm"),
        merge_legends=True,
        heatmap_legend_side="bottom",
    )


def render_2022_fig4b_dmr_summary() -> None:
    """iMeta (2022) Figure 4B: Global summary statistics in a multi-omics epigenomics study."""
    reset_state()
    dmr = ch.load_dataset("dmr_summary")
    label = dmr["label"]["value"].astype(str).tolist()
    n_gr = dmr["n_gr"]["value"].to_numpy(dtype=float)
    mean_meth = dmr["mean_meth"]
    n_corr = dmr["n_corr"]
    dist_tss = dmr["dist_tss"]
    gene_anno = dmr["gene_anno"]
    cgi_anno = dmr["cgi_anno"]
    mat_enrich_gf = dmr["mat_enrich_gf"]
    mat_pct_st = dmr["mat_pct_st"]
    mat_enrich_st = dmr["mat_enrich_st"]

    meth_col_fun = colorRamp2([0, 0.5, 1], ["blue", "white", "red"])
    corr_col = ["green", "red"]
    dist_tss_col = ["#FF0000", "#FF7352", "#FFB299", "#FFD9CB"]
    gene_anno_col = ["green", "blue"]
    cgi_anno_col = ["#FFA500", "#FFD191"]
    z_score_col_fun = colorRamp2([-200, 0, 200], ["green", "white", "red"])
    state_col = ["#FF0000", "#008000", "#C2E105", "#8A91D0", "#CD5C5C", "#808080", "#000000"]

    anno_w = unit(2.8, "cm")
    ht_list = (
        rowAnnotation(text=anno_text(label, location=unit(1, "npc"), just="right", gp=gpar(fontsize=10.5)))
        + Heatmap(
            mean_meth,
            name="mean_meth",
            col=meth_col_fun,
            cluster_rows=False,
            row_title=None,
            cluster_columns=False,
            show_row_names=False,
            column_names_rot=45,
            column_names_gp=gpar(fontsize=8.5),
            heatmap_legend_param={"title": "Methylation", "legend_direction": "horizontal", "legend_width": unit(2.8, "cm")},
            width=unit(mean_meth.shape[1] * 4.5, "mm"),
        )
        + rowAnnotation(n_gr=anno_barplot(n_gr, bar_width=1, width=anno_w), show_annotation_name=False)
        + rowAnnotation(n_corr=anno_barplot(n_corr.to_numpy(), bar_width=1, gp=gpar(fill=corr_col), width=anno_w), show_annotation_name=False)
        + rowAnnotation(dist_tss=anno_barplot(dist_tss.to_numpy(), bar_width=1, gp=gpar(fill=dist_tss_col), width=anno_w), show_annotation_name=False)
        + rowAnnotation(gene_anno=anno_barplot(gene_anno.to_numpy(), bar_width=1, gp=gpar(fill=gene_anno_col), width=anno_w), show_annotation_name=False)
        + rowAnnotation(cgi_anno=anno_barplot(cgi_anno.to_numpy(), bar_width=1, gp=gpar(fill=cgi_anno_col), width=anno_w), show_annotation_name=False)
        + Heatmap(
            mat_enrich_gf,
            name="enrich_gf",
            col=z_score_col_fun,
            cluster_columns=False,
            width=unit(mat_enrich_gf.shape[1] * 4.2, "mm"),
            column_title="",
            column_names_rot=45,
            column_names_gp=gpar(fontsize=8.5),
            show_heatmap_legend=False,
        )
        + rowAnnotation(pct_st=anno_barplot(mat_pct_st.to_numpy(), bar_width=1, gp=gpar(fill=state_col), width=anno_w), show_annotation_name=False)
        + Heatmap(
            mat_enrich_st,
            name="enrich_st",
            col=z_score_col_fun,
            cluster_columns=False,
            width=unit(mat_enrich_st.shape[1] * 4.2, "mm"),
            column_title="",
            show_heatmap_legend=False,
            column_names_gp=gpar(col=state_col, fontsize=8.5),
            show_row_names=False,
            column_names_rot=45,
        )
    )

    lgd_list = [
        Legend(labels=["negative", "positive"], title="Correlation", legend_gp=gpar(fill=corr_col)),
        Legend(labels=["gene", "intergenic"], title="Gene annotation", legend_gp=gpar(fill=gene_anno_col)),
        Legend(labels=["<1kb", "1kb~5kb", "5kb~10kb", ">10kb"], title="Distance to TSS", legend_gp=gpar(fill=dist_tss_col), nrow=2),
        Legend(labels=["CGI", "CGI shore"], title="CGI annotation", legend_gp=gpar(fill=cgi_anno_col)),
        Legend(col_fun=z_score_col_fun, title="Z-score", direction="horizontal", at=[-200, 0, 200], legend_width=unit(2.8, "cm")),
        Legend(labels=list(mat_enrich_st.columns), title="Chromatin states", legend_gp=gpar(fill=state_col), nrow=2),
    ]

    split_grp = [re.sub(r"\d+$", "", s) for s in label]

    def _decorate_titles() -> None:
        anno_titles = {
            "n_gr": "Number of\nDMRs",
            "n_corr": "Significantly\ncorrelated genes",
            "dist_tss": "Distance to TSS",
            "gene_anno": "Gene annotation",
            "cgi_anno": "CGI annotation",
            "pct_st": "Overlap to\nChromatin states",
        }
        for an, ttl in anno_titles.items():
            with decorate_annotation(an, slice=1):
                grid_text(ttl, y=unit(1, "npc") + unit(3, "mm"), just="bottom", gp=gpar(fontsize=9))

        ht_titles = {
            "mean_meth": "Mean\nmethylation",
            "enrich_gf": "Enrichment to\ngenomic features",
            "enrich_st": "Enrichment to\nchromatin states",
        }
        for hn, ttl in ht_titles.items():
            with decorate_heatmap_body(hn, slice=1):
                grid_text(ttl, y=unit(1, "npc") + unit(3, "mm"), just="bottom", gp=gpar(fontsize=9))

        dev = ch.current_device(False)
        while len(dev.vp_stack) > 1:
            popViewport()
        grid_text(
            "B) Visualize the global summary statistics in a multi-omics study",
            x=unit(0.5, "npc"),
            y=unit(1, "npc") - unit(3.5, "mm"),
            just="top",
            gp=gpar(fontsize=11.5, fontface="bold"),
        )

    save_decorated(
        ht_list,
        "ch16_2022_fig4b_dmr_summary.png",
        _decorate_titles,
        width=15.0,
        height=6.6,
        padding=unit([4, 4, 24, 4], "mm"),
        row_split=split_grp,
        heatmap_legend_list=lgd_list,
        heatmap_legend_side="bottom",
    )


def render_2022_fig5ab_density_freq3d() -> None:
    """iMeta (2022) Figure 5A & 5B: densityHeatmap and 3D frequencyHeatmap."""
    reset_state()
    rng = np.random.default_rng(123)
    m = np.hstack([rng.standard_normal((100, 10)), rng.uniform(-2.0, 2.0, (100, 10)) + 0.5])
    df = pd.DataFrame(m, columns=[f"C{i}" for i in range(1, 21)])

    def _draw_a() -> None:
        ch.ht_opt.HEATMAP_LEGEND_PADDING = unit(12, "mm")
        ha_a = HeatmapAnnotation(
            Distribution=["rnorm"] * 10 + ["runif"] * 10,
            col={"Distribution": {"rnorm": "#DF536B", "runif": "#61D04F"}},
        )
        dh = densityHeatmap(df, ylab="Value", top_annotation=ha_a, title="A) Density heatmap")
        ch.draw(dh, newpage=False, merge_legends=True, padding=unit([4, 16, 4, 4], "mm"))

    def _draw_b() -> None:
        ch.ht_opt.HEATMAP_LEGEND_PADDING = unit(8, "mm")
        ch.ht_opt.TITLE_PADDING = unit([10, 2], "mm")
        ha_b = HeatmapAnnotation(
            Distribution=["rnorm"] * 10 + ["runif"] * 10,
            col={"Distribution": {"rnorm": "#DF536B", "runif": "#61D04F"}},
        )
        fh3d = frequencyHeatmap(
            df,
            ylab="Value",
            top_annotation=ha_b,
            use_3d=True,
            title="B) 3D frequency heatmap",
            bar_max_length=unit(7, "mm"),
        )
        ch.draw(fh3d, newpage=False, merge_legends=True, padding=unit([4, 8, 6, 4], "mm"))

    gb_a = grid_grabExpr(_draw_a, width_in=6.8, height_in=5.2)
    reset_state()
    gb_b = grid_grabExpr(_draw_b, width_in=6.8, height_in=5.2)
    reset_state()

    def _draw_ab() -> None:
        pushViewport(Viewport(x=0.0, width=0.5, just="left"))
        grid_draw(gb_a)
        popViewport()
        pushViewport(Viewport(x=0.5, width=0.5, just="left"))
        grid_draw(gb_b)
        popViewport()

    save_Standalone(_draw_ab, "ch16_2022_fig5ab_density_freq3d.png", width=13.6, height=5.3, dpi=140)


def render_2022_fig5c_oncoprint() -> None:
    """iMeta (2022) Figure 5C: OncoPrint of 22 pathway genes across 80 TCGA LUAD patients."""
    reset_state()
    mat = ch.load_dataset("tcga_lung_adenocarcinoma")
    sample_order = ch.load_dataset("sample_order")
    if isinstance(sample_order, pd.DataFrame):
        sample_order = sample_order.iloc[:, 0].astype(str).tolist()
    elif isinstance(sample_order, pd.Series):
        sample_order = sample_order.astype(str).tolist()

    rng = np.random.default_rng(123)
    kept_cols = rng.choice(mat.columns.tolist(), size=80, replace=False).tolist()
    sub_mat = mat[kept_cols]
    non_empty_r = sub_mat.apply(lambda r: any(str(x).strip() != "" for x in r), axis=1)
    non_empty_c = sub_mat.apply(lambda c: any(str(x).strip() != "" for x in c), axis=0)
    sub_mat = sub_mat.loc[non_empty_r, non_empty_c]
    kept_set = set(sub_mat.columns.tolist())
    sub_order = [s for s in sample_order if s in kept_set]

    col = {"HOMDEL": "blue", "AMP": "red", "MUT": "#008000"}
    alter_fun = {
        "background": lambda x, y, w, h: grid_rect(x, y, w - unit(1.5, "pt"), h - unit(1.5, "pt"), gp=gpar(fill="#CCCCCC", col=None)),
        "HOMDEL": lambda x, y, w, h: grid_rect(x, y, w - unit(1.5, "pt"), h - unit(1.5, "pt"), gp=gpar(fill=col["HOMDEL"], col=None)),
        "AMP": lambda x, y, w, h: grid_rect(x, y, w - unit(1.5, "pt"), h - unit(1.5, "pt"), gp=gpar(fill=col["AMP"], col=None)),
        "MUT": lambda x, y, w, h: grid_rect(x, y, w - unit(1.5, "pt"), h * 0.35, gp=gpar(fill=col["MUT"], col=None)),
    }

    op = oncoPrint(
        sub_mat,
        alter_fun=alter_fun,
        col=col,
        row_order=list(range(sub_mat.shape[0])),
        column_order=sub_order if len(sub_order) == sub_mat.shape[1] else None,
        remove_empty_columns=False,
        remove_empty_rows=False,
        column_title="C) OncoPrint",
        column_title_gp=gpar(fontsize=11.5, fontface="bold"),
        row_names_gp=gpar(fontsize=8.5),
        pct_gp=gpar(fontsize=8.5),
        heatmap_legend_param={
            "title": "Alterations",
            "at": ["HOMDEL", "AMP", "MUT"],
            "labels": ["Deep deletion", "Amplification", "Mutation"],
            "title_position": "leftcenter",
            "nrow": 1,
        },
    )
    save_ht(op, "ch16_2022_fig5c_oncoprint.png", width=8.6, height=6.2, heatmap_legend_side="bottom", padding=unit([4, 6, 4, 6], "mm"))


def render_2022_fig5d_upset() -> None:
    """iMeta (2022) Figure 5D: Genomic UpSet plot of Roadmap H3K4me3 peaks across 6 human tissues."""
    reset_state()
    sets_df = pd.read_csv(DATA_DIR / "pub2022_fig5d_h3k4me3_sets.csv")
    combs_df = pd.read_csv(DATA_DIR / "pub2022_fig5d_h3k4me3_combs.csv", dtype={"comb_code": str})

    set_names = sets_df["set_name"].tolist()
    set_sizes = sets_df["set_size"].to_numpy(dtype=float)
    comb_codes = combs_df["comb_code"].tolist()
    comb_sizes = combs_df["comb_size"].to_numpy(dtype=float)
    mat = combs_df[set_names].to_numpy(dtype=int).T

    cm = CombMat(
        mat=mat,
        set_names=set_names,
        comb_codes=comb_codes,
        set_sizes=set_sizes,
        comb_sizes=comb_sizes,
        mode="distinct",
    )

    set_ord = np.argsort(-set_sizes).tolist()
    comb_ord = np.argsort(-comb_sizes).tolist()

    top_anno = upset_top_annotation(
        cm,
        height=unit(3.6, "cm"),
        axis_param={"at": [0, 1e7, 2e7], "labels": ["0 MB", "10 MB", "20 MB"], "side": "left", "gp": gpar(fontsize=8)},
    )
    right_anno = upset_right_annotation(
        cm,
        width=unit(2.8, "cm"),
        annotation_name_offset=unit(8, "mm"),
        axis_param={"at": [0, 2e7, 4e7, 6e7], "labels": ["0", "20", "40", "60 MB"], "side": "bottom", "labels_rot": 0, "gp": gpar(fontsize=8)},
    )

    up = UpSet(
        cm,
        set_order=set_ord,
        comb_order=comb_ord,
        top_annotation=top_anno,
        right_annotation=right_anno,
        column_title="D) UpSet plot (Roadmap H3K4me3 ChIP-seq peaks)",
        column_title_gp=gpar(fontsize=11.5, fontface="bold"),
        row_names_gp=gpar(fontsize=9.5),
    )
    save_ht(up, "ch16_2022_fig5d_upset.png", width=9.6, height=5.4, padding=unit([12, 8, 4, 8], "mm"))


def render_2022_fig5e_genome() -> None:
    """iMeta (2022) Figure 5E: Genome-level multi-track vertical heatmap across chr1..chr22."""
    reset_state()
    chrs = [f"chr{i}" for i in range(1, 23)]
    chr_window = ch.bin_genome("hg19", bin_size=1_000_000, chromosomes=chrs)

    bed1 = ch.generateRandomBed(nr=1000, nc=10, species="hg19", seed=123)
    bed1 = bed1[bed1["chr"].isin(chrs)].reset_index(drop=True)
    gr1 = ch.GRanges(
        bed1["chr"].tolist(),
        ch.IRanges(bed1["start"].to_numpy(), bed1["end"].to_numpy()),
        **{c: bed1[c].to_numpy() for c in bed1.columns[3:]},
    )
    num_mat = ch.normalize_genomic_signals_to_bins(
        gr1, list(bed1.columns[3:]), bins=chr_window, method="weighted", empty_value=np.nan
    ).mcols.to_numpy(dtype=float)

    rng = np.random.default_rng(42)
    char_cols = []
    for j in range(10):
        bed_c = ch.generateRandomBed(nr=1000, nc=0, species="hg19", seed=100 + j)
        bed_c = bed_c[bed_c["chr"].isin(chrs)].sample(n=20, random_state=100 + j).reset_index(drop=True)
        cnv_vals = rng.choice(["gain", "loss"], size=len(bed_c))
        gr_c = ch.GRanges(bed_c["chr"].tolist(), ch.IRanges(bed_c["start"].to_numpy(), bed_c["end"].to_numpy()))
        out = np.full(len(chr_window), np.nan, dtype=object)
        qh, sh = ch.findOverlaps(chr_window, gr_c)
        for qi, si in zip(qh, sh):
            out[int(qi)] = str(cnv_vals[int(si)])
        char_cols.append(out)
    char_mat = np.column_stack(char_cols)

    bed2 = ch.generateRandomBed(nr=120, nc=2, species="hg19", seed=456)
    bed2 = bed2[bed2["chr"].isin(chrs)].reset_index(drop=True)
    gr2 = ch.GRanges(bed2["chr"].tolist(), ch.IRanges(bed2["start"].to_numpy(), bed2["end"].to_numpy()), v1=bed2["value1"].to_numpy(), v2=bed2["value2"].to_numpy())
    v_mat = ch.normalize_genomic_signals_to_bins(gr2, ["v1", "v2"], bins=chr_window, method="weighted", empty_value=np.nan).mcols.to_numpy(dtype=float)

    bed3 = ch.generateRandomBed(nr=42, nc=0, species="hg19", seed=789)
    bed3 = bed3[bed3["chr"].isin(chrs)].reset_index(drop=True)
    gr3 = ch.GRanges(bed3["chr"].tolist(), ch.IRanges(bed3["start"].to_numpy(), bed3["start"].to_numpy()), gene=[f"gene_{k + 1}" for k in range(len(bed3))])
    qh3, sh3 = ch.findOverlaps(chr_window, gr3)
    at = qh3.tolist()
    labels = gr3.mcols["gene"].iloc[sh3].tolist()

    chr_cat = pd.Categorical(chr_window.seqnames, categories=chrs, ordered=True)
    subgroup = ["A"] * 5 + ["B"] * 5
    bar_vals = np.nan_to_num(v_mat[:, 0], nan=0.0)
    bar_cols = ["#DF536B" if val > 0 else "#61D04F" for val in bar_vals]
    alt_titles = [f" \n{c}" if (idx + 1) % 2 == 0 else f"{c}\n " for idx, c in enumerate(chrs)]

    with ch.ht_opt.local(TITLE_PADDING=unit([6, 6], "points")):
        ht_vert = (
            Heatmap(
                num_mat.T,
                name="mat",
                col=colorRamp2([-1, 0, 1], ["green", "white", "red"]),
                na_col="#F2F2F2",
                column_split=chr_cat,
                cluster_columns=False,
                show_row_dend=False,
                row_split=subgroup,
                cluster_row_slices=False,
                row_title="%s",
                left_annotation=rowAnnotation(
                    Subgroup=subgroup,
                    col={"Subgroup": {"A": "#DF536B", "B": "#61D04F"}},
                    show_annotation_name=False,
                    annotation_legend_param={"Subgroup": {"legend_direction": "horizontal", "title_position": "lefttop", "nrow": 1}},
                ),
                column_title_gp=gpar(fontsize=7.5),
                border=True,
                column_gap=unit(0, "points"),
                column_title=alt_titles,
                heatmap_legend_param={"legend_direction": "horizontal", "title_position": "lefttop"},
                use_raster=False,
            )
            % v %
            Heatmap(
                char_mat.T,
                name="CNV",
                col={"gain": "red", "loss": "blue"},
                na_col="#EFEFEF",
                border=True,
                cluster_rows=False,
                row_title="CNV",
                heatmap_legend_param={"legend_direction": "horizontal", "title_position": "lefttop", "nrow": 1},
                use_raster=False,
            )
            % v %
            HeatmapAnnotation(label=anno_mark(at=at, labels=labels, side="bottom", labels_gp=gpar(fontsize=7)))
            % v %
            HeatmapAnnotation(
                points=anno_points(v_mat, gp=gpar(col=["#2297E6", "#28E2E5"]), pch=[1, 16], size=unit(1.4, "mm")),
                annotation_name_side="left",
                height=unit(1.8, "cm"),
            )
            % v %
            HeatmapAnnotation(
                bars=anno_barplot(bar_vals, gp=gpar(col=bar_cols, fill=bar_cols)),
                annotation_name_side="left",
                height=unit(1.8, "cm"),
            )
        )
        save_ht(
            ht_vert,
            "ch16_2022_fig5e_genome.png",
            width=12.0,
            height=7.8,
            padding=unit([4, 4, 8, 8], "mm"),
            heatmap_legend_side="bottom",
            merge_legends=True,
            column_title="E) Genome-level multi-track heatmap",
            column_title_gp=gpar(fontsize=11.5, fontface="bold"),
        )


def render_2022_fig6_tss_enriched() -> None:
    """iMeta (2022) Figure 6 (ComplexHeatmap-py TSS spatial enrichment reconstruction):

    Displays normalized chromatin states, DNA methylation (-5kb to +5kb around TSS),
    and Roadmap E003 gene expression across active, bivalent, and inactive promoters,
    together with top enrichment profile lines across the three groups.
    """
    reset_state()
    npz = np.load(DATA_DIR / "pub2022_fig6_tss_roadmap.npz")
    states_idx = npz["states_idx"]
    state_names = [str(s) for s in npz["state_names"]]
    meth = npz["meth"]
    expr = npz["expr"]
    groups = [str(g) for g in npz["group"]]

    # Reconstruct 2D string state array matching R ChromHMM simplified categories
    states_str = np.asarray(state_names, dtype=object)[states_idx]

    states_col = {
        "TssActive": "#E41A1C",
        "TssBivalent": "#FF7F00",
        "Transcript": "#4DAF4A",
        "Enhancer": "#FFFF33",
        "Heterochromatin": "#AFEEEE",
        "Repressive": "#999999",
        "Quiescent": "#1A1A1A",
    }

    split_levels = ["active", "bivalent", "inactive"]
    row_split = pd.Categorical(groups, categories=split_levels, ordered=True)

    # Sort within each group by center methylation + TSS activity for clean visual banding
    order_within = []
    for lvl in split_levels:
        idx_l = np.where(np.asarray(groups) == lvl)[0]
        score = np.mean(states_idx[idx_l, 45:55] == 0, axis=1) * 10.0 - np.mean(meth[idx_l, 45:55], axis=1)
        order_within.extend(idx_l[np.argsort(-score)].tolist())

    states_ord = states_str[order_within]
    meth_ord = meth[order_within]
    expr_ord = expr[order_within]
    split_ord = pd.Categorical([groups[i] for i in order_within], categories=split_levels, ordered=True)

    # Precompute group-level enrichment profiles over 100 bins (-5000..+5000)
    # States enrichment: fraction of active + bivalent promoter chromatin states at each bin
    ltys = {"active": "solid", "bivalent": "dashed", "inactive": "dotted"}
    bins_x = np.arange(states_ord.shape[1], dtype=float)
    state_profiles = {}
    meth_profiles = {}
    for lvl in split_levels:
        mask = np.asarray(split_ord) == lvl
        prom_frac = np.mean(np.isin(states_ord[mask], ["TssActive", "TssBivalent"]), axis=0)
        state_profiles[lvl] = prom_frac
        meth_profiles[lvl] = np.mean(meth_ord[mask], axis=0)

    def _make_profile_anno(profiles: dict[str, np.ndarray], ymax: float, ylab: str = "") -> AnnotationFunction:
        def _draw(index, k=1, n_slices=1, env=None):
            pushViewport(Viewport(xscale=(0.0, 99.0), yscale=(0.0, ymax), clip=False))
            try:
                grid_rect(gp=gpar(fill="white", col="black", lwd=0.8))
                grid_lines([49.5, 49.5], [0.0, ymax], gp=gpar(lty="dotted", col="#999999", lwd=0.8), default_units="native")
                for lvl in split_levels:
                    ch.grid_lines(bins_x, profiles[lvl], gp=gpar(lty=ltys[lvl], col="black", lwd=1.3), default_units="native")
                grid_yaxis(main=False, gp=gpar(fontsize=7.5))
            finally:
                popViewport()

        return AnnotationFunction(
            fun=_draw,
            fun_name="tss_profile",
            which="column",
            width=unit(1.0, "npc"),
            height=unit(2.4, "cm"),
            n=100,
            subsettable=False,
            _which_explicit=True,
        )

    pos_labels = ["-5000" if j == 0 else ("start" if j == 49 else ("5000" if j == 99 else "")) for j in range(100)]
    ha_tss_bot = HeatmapAnnotation(
        pos=anno_text(pos_labels, rot=0, just="center", location=unit(0.75, "npc"), gp=gpar(fontsize=8.5)),
        height=unit(6, "mm"),
        show_annotation_name=False,
    )
    ha_tss_bot2 = HeatmapAnnotation(
        pos=anno_text(pos_labels, rot=0, just="center", location=unit(0.75, "npc"), gp=gpar(fontsize=8.5)),
        height=unit(6, "mm"),
        show_annotation_name=False,
    )

    ht_states = Heatmap(
        states_ord,
        name="States",
        col=states_col,
        cluster_rows=False,
        cluster_columns=False,
        row_split=split_ord,
        row_title="%s",
        row_title_gp=gpar(fontsize=10, fontface="bold"),
        show_row_names=False,
        show_column_names=False,
        top_annotation=HeatmapAnnotation(enrich=_make_profile_anno(state_profiles, 1.0), show_annotation_name=False),
        bottom_annotation=ha_tss_bot,
        width=unit(4.4, "cm"),
        use_raster=True,
    )
    ht_meth = Heatmap(
        meth_ord,
        name="Methylation",
        col=colorRamp2([0.0, 0.5, 1.0], ["blue", "white", "red"]),
        cluster_rows=False,
        cluster_columns=False,
        show_row_names=False,
        show_column_names=False,
        top_annotation=HeatmapAnnotation(enrich=_make_profile_anno(meth_profiles, 0.9), show_annotation_name=False),
        bottom_annotation=ha_tss_bot2,
        width=unit(4.4, "cm"),
        use_raster=True,
    )
    ht_expr = Heatmap(
        expr_ord,
        name="Expression",
        col=colorRamp2([0.0, 4.0, 8.0], ["white", "#9E9AC8", "#E41A1C"]),
        cluster_rows=False,
        cluster_columns=False,
        show_row_names=False,
        show_column_names=False,
        width=unit(1.0, "cm"),
        top_annotation=HeatmapAnnotation(
            summary=anno_summary(gp=gpar(fill="#E8E8E8", col="black"), axis_param={"side": "right", "gp": gpar(fontsize=7.5)}, height=unit(2.4, "cm")),
            show_annotation_name=False,
        ),
        use_raster=False,
    )

    lgd_cat = Legend(
        title="Category",
        labels=split_levels,
        legend_gp=gpar(lty=["solid", "dashed", "dotted"], col="black", lwd=1.3),
        type="lines",
        grid_width=unit(9, "mm"),
    )

    def _decorate_tss() -> None:
        for hn in ("States", "Methylation"):
            for sl in (1, 2, 3):
                with decorate_heatmap_body(hn, slice=sl):
                    grid_lines([0.5, 0.5], [0.0, 1.0], gp=gpar(lty="dotted", col="#222222", lwd=0.8), default_units="npc")

    save_decorated(
        ht_states + ht_meth + ht_expr,
        "ch16_2022_fig6_tss_enriched.png",
        _decorate_tss,
        width=9.6,
        height=9.2,
        ht_gap=unit(10, "mm"),
        column_title="Visualize chromatin states, methylation and expression",
        column_title_gp=gpar(fontsize=11.5, fontface="bold"),
        heatmap_legend_list=[lgd_cat],
        padding=unit([6, 12, 8, 6], "mm"),
    )


def main() -> None:
    print("Rendering 2016 Bioinformatics Fig 1A...")
    render_2016_fig1a_oncoprint_go()
    print("Rendering 2016 Bioinformatics Fig 1B...")
    render_2016_fig1b_scrnaseq_721()
    print("Rendering 2022 iMeta Fig 1...")
    render_2022_fig1_design()
    print("Rendering 2022 iMeta Fig 2A-B...")
    render_2022_fig2ab_single()
    print("Rendering 2022 iMeta Fig 2C-D...")
    render_2022_fig2cd_presplit()
    print("Rendering 2022 iMeta Fig 3A...")
    render_2022_fig3a_gallery()
    print("Rendering 2022 iMeta Fig 3B-D...")
    render_2022_fig3bcd_complex_annos()
    print("Rendering 2022 iMeta Fig 4A...")
    render_2022_fig4a_meth_list()
    print("Rendering 2022 iMeta Fig 4B...")
    render_2022_fig4b_dmr_summary()
    print("Rendering 2022 iMeta Fig 5A-B...")
    render_2022_fig5ab_density_freq3d()
    print("Rendering 2022 iMeta Fig 5C...")
    render_2022_fig5c_oncoprint()
    print("Rendering 2022 iMeta Fig 5D...")
    render_2022_fig5d_upset()
    print("Rendering 2022 iMeta Fig 5E...")
    render_2022_fig5e_genome()
    print("Rendering 2022 iMeta Fig 6 TSS analog...")
    render_2022_fig6_tss_enriched()
    reset_state()
    print("Chapter 16 publication figures generated successfully.")


if __name__ == "__main__":
    main()
