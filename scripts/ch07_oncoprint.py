"""Generate all figures for Chapter 7: OncoPrint."""

from __future__ import annotations

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import complexheatmap as ch
from complexheatmap import (
    Heatmap,
    HeatmapAnnotation,
    alter_graphic,
    anno_barplot,
    anno_oncoprint_barplot,
    anno_points,
    default_get_type,
    gpar,
    grid,
    oncoPrint,
    rowAnnotation,
    test_alter_fun,
    unify_mat_list,
    unit,
    unit_c,
    v,
)

from _common import IMAGES_DIR, reset_state, save_ht


def main() -> None:
    # ---- 7.1 General settings ----
    csv_txt = (
        "gene,s1,s2,s3\n"
        "g1,snv;indel,snv,indel\n"
        "g2,,snv;indel,snv\n"
        "g3,snv,,indel;snv\n"
    )
    mat = pd.read_csv(io.StringIO(csv_txt), index_col=0, keep_default_na=False)
    print("default_get_type(g1, s1):", default_get_type(mat.iloc[0, 0]))
    print("default_get_type(g1, s2):", default_get_type(mat.iloc[0, 1]))

    col = {"snv": "red", "indel": "blue"}
    alter_fun_basic = {
        "snv": lambda x, y, w, h: grid.rect(x, y, w * 0.9, h * 0.9, gp=gpar(fill=col["snv"], col=None)),
        "indel": lambda x, y, w, h: grid.rect(x, y, w * 0.9, h * 0.4, gp=gpar(fill=col["indel"], col=None)),
    }

    reset_state()
    save_ht(oncoPrint(mat, alter_fun=alter_fun_basic, col=col), "ch07_01_basic_str_mat.png", width=5.0, height=4.4)

    # Dict of binary matrices + unify_mat_list
    reset_state()
    m_snv = pd.DataFrame([[1, 1, 0], [0, 1, 1], [1, 0, 1]], index=["g1", "g2", "g3"], columns=["s1", "s2", "s3"])
    m_indel = pd.DataFrame([[1, 0, 1], [0, 1, 0], [0, 0, 1]], index=["g1", "g2", "g3"], columns=["s1", "s2", "s3"])
    mat_list = unify_mat_list({"snv": m_snv, "indel": m_indel}, default=0)
    save_ht(oncoPrint(mat_list, alter_fun=alter_fun_basic, col=col), "ch07_02_mat_list_unify.png", width=5.0, height=4.4)

    # Single-function alter_fun(x, y, w, h, v)
    reset_state()
    def _single_af(x, y, w, h, v):
        types = [] if v == "background" else ([v] if isinstance(v, str) else list(v))
        if "snv" in types:
            grid.rect(x, y, w * 0.9, h * 0.9, gp=gpar(fill=col["snv"], col=None))
        if "indel" in types:
            grid.rect(x, y, w * 0.9, h * 0.4, gp=gpar(fill=col["indel"], col=None))
    save_ht(oncoPrint(mat, alter_fun=_single_af, col=col), "ch07_03_single_alter_fun.png", width=5.0, height=4.4)

    # Dynamic-height split rectangle per cell
    reset_state()
    def _stacked_af(x, y, w, h, v):
        types = [] if v == "background" else ([v] if isinstance(v, str) else list(v))
        n = len(types)
        h_box = h * 0.9
        if n > 0:
            for idx, tp in enumerate(types, start=1):
                grid.rect(
                    x,
                    y - h_box * 0.5 + (idx / n) * h_box,
                    w * 0.9,
                    (1.0 / n) * h_box,
                    gp=gpar(fill=col[tp], col=None),
                    just="top",
                )
    save_ht(oncoPrint(mat, alter_fun=_stacked_af, col=col), "ch07_04_dynamic_height_alter_fun.png", width=5.0, height=4.4)

    # Triangle alter_fun
    reset_state()
    alter_fun_tri = {
        "background": lambda x, y, w, h: (
            grid.polygon(unit_c(x - 0.5 * w, x - 0.5 * w, x + 0.5 * w), unit_c(y - 0.5 * h, y + 0.5 * h, y - 0.5 * h), gp=gpar(fill="grey", col="white")),
            grid.polygon(unit_c(x + 0.5 * w, x + 0.5 * w, x - 0.5 * w), unit_c(y + 0.5 * h, y - 0.5 * h, y + 0.5 * h), gp=gpar(fill="grey", col="white")),
        ),
        "snv": lambda x, y, w, h: grid.polygon(
            unit_c(x - 0.5 * w, x - 0.5 * w, x + 0.5 * w),
            unit_c(y - 0.5 * h, y + 0.5 * h, y - 0.5 * h),
            gp=gpar(fill=col["snv"], col="white"),
        ),
        "indel": lambda x, y, w, h: grid.polygon(
            unit_c(x + 0.5 * w, x + 0.5 * w, x - 0.5 * w),
            unit_c(y + 0.5 * h, y - 0.5 * h, y + 0.5 * h),
            gp=gpar(fill=col["indel"], col="white"),
        ),
    }
    save_ht(oncoPrint(mat, alter_fun=alter_fun_tri, col=col), "ch07_05_triangle_alter_fun.png", width=5.0, height=4.0)

    # test_alter_fun preview
    reset_state()
    alter_fun_7 = {
        "mut1": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill="red", col=None)),
        "mut2": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill="blue", col=None)),
        "mut3": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill="yellow", col=None)),
        "mut4": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill="purple", col=None)),
        "mut5": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill=None, lwd=2)),
        "mut6": lambda x, y, w, h: grid.points(x, y, pch=16),
        "mut7": lambda x, y, w, h: grid.segments(x - w * 0.5, y - h * 0.5, x + w * 0.5, y + h * 0.5, gp=gpar(lwd=2)),
    }
    dev_test = test_alter_fun(alter_fun_7, width=4.0, height=3.4)
    dev_test.fig.savefig(IMAGES_DIR / "ch07_06_test_alter_fun.png", dpi=140, facecolor="white")
    plt.close("all")

    # Background variants
    reset_state()
    af_green_bg = {
        "background": lambda x, y, w, h: grid.rect(x, y, w, h, gp=gpar(fill="#00FF0020")),
        **alter_fun_basic,
    }
    save_ht(oncoPrint(mat, alter_fun=af_green_bg, col=col), "ch07_07_bg_green.png", width=5.0, height=4.4)

    reset_state()
    af_no_bg = {
        "background": lambda *args: None,
        **alter_fun_basic,
    }
    save_ht(oncoPrint(mat, alter_fun=af_no_bg, col=col), "ch07_08_bg_none.png", width=5.0, height=4.4)

    # Complex alteration types (primary + secondary classes)
    reset_state()
    rng = np.random.default_rng(123)
    x1 = rng.choice(["", "snv"], 100, replace=True, p=[0.8, 0.2])
    x2 = rng.choice(["", "indel"], 100, replace=True, p=[0.8, 0.2])
    x2[x1 == "snv"] = ""
    x3 = rng.choice(["", "intronic"], 100, replace=True, p=[0.5, 0.5])
    x4 = rng.choice(["", "exonic"], 100, replace=True, p=[0.5, 0.5])
    x3[(x1 == "") & (x2 == "")] = ""
    x4[(x1 == "") & (x2 == "")] = ""
    x4[x3 == "intronic"] = ""
    cells = [";".join([t for t in row if t != ""]) for row in zip(x1, x2, x3, x4)]
    m_complex = pd.DataFrame(
        np.array(cells, dtype=object).reshape(10, 10),
        index=[f"g{i}" for i in range(1, 11)],
        columns=[f"s{i}" for i in range(1, 11)],
    )
    af_complex = {
        "background": lambda x, y, w, h: grid.rect(x, y, w * 0.9, h * 0.9, gp=gpar(fill="#CCCCCC", col=None)),
        "snv": lambda x, y, w, h: grid.rect(x, y, w * 0.9, h * 0.9, gp=gpar(fill="red", col=None)),
        "indel": lambda x, y, w, h: grid.rect(x, y, w * 0.9, h * 0.9, gp=gpar(fill="blue", col=None)),
        "intronic": lambda x, y, w, h: grid.points(x, y, pch=16),
        "exonic": lambda x, y, w, h: (
            grid.segments(x - w * 0.4, y - h * 0.4, x + w * 0.4, y + h * 0.4, gp=gpar(lwd=2)),
            grid.segments(x + w * 0.4, y - h * 0.4, x - w * 0.4, y + h * 0.4, gp=gpar(lwd=2)),
        ),
    }
    op_complex = oncoPrint(
        m_complex,
        alter_fun=af_complex,
        col={"snv": "red", "indel": "blue"},
        top_annotation=HeatmapAnnotation(cbar=anno_oncoprint_barplot(["snv", "indel"])),
        right_annotation=rowAnnotation(rbar=anno_oncoprint_barplot(["snv", "indel"], which="row")),
        heatmap_legend_param={"title": "Alterations", "at": ["snv", "indel"], "labels": ["snv", "indel"]},
    )
    save_ht(op_complex, "ch07_09_complex_alt_types.png", width=6.4, height=5.4)

    # Simplify alter_fun with alter_graphic()
    reset_state()
    op_ag = oncoPrint(
        mat,
        alter_fun={
            "snv": alter_graphic("rect", width=0.9, height=0.9, fill=col["snv"]),
            "indel": alter_graphic("rect", width=0.9, height=0.4, fill=col["indel"]),
        },
        col=col,
    )
    save_ht(op_ag, "ch07_10_alter_graphic.png", width=5.0, height=4.4)

    # Other heatmap-related settings
    reset_state()
    save_ht(
        oncoPrint(mat, alter_fun=alter_fun_basic, col=col, show_column_names=True),
        "ch07_11_show_column_names.png",
        width=5.0,
        height=4.6,
    )

    reset_state()
    save_ht(
        oncoPrint(
            mat, alter_fun=alter_fun_basic, col=col,
            row_names_side="left", pct_side="right", pct_digits=2,
        ),
        "ch07_12_pct_and_row_names_side.png",
        width=5.4,
        height=4.4,
    )

    reset_state()
    save_ht(
        oncoPrint(
            mat,
            alter_fun=alter_fun_basic,
            col=col,
            top_annotation=HeatmapAnnotation(cbar=anno_oncoprint_barplot(height=unit(1, "cm"))),
            right_annotation=rowAnnotation(
                rbar=anno_oncoprint_barplot(
                    which="row",
                    width=unit(4, "cm"),
                    axis_param={"at": [0, 2, 4], "labels": ["zero", "two", "four"], "side": "top", "labels_rot": 0},
                )
            ),
        ),
        "ch07_13_custom_barplot_axis.png",
        width=6.2,
        height=4.8,
    )

    reset_state()
    save_ht(
        oncoPrint(
            mat,
            alter_fun=alter_fun_basic,
            col=col,
            left_annotation=rowAnnotation(
                rbar=anno_oncoprint_barplot(which="row", axis_param={"direction": "reverse"})
            ),
            right_annotation=None,
        ),
        "ch07_14_left_barplot.png",
        width=5.6,
        height=4.6,
    )

    # ---- 7.2 Apply to cBioPortal TCGA Lung Adenocarcinoma dataset ----
    tcga = ch.load_dataset("tcga_lung_adenocarcinoma")
    sample_order = ch.load_dataset("sample_order")
    col_tcga = {"HOMDEL": "blue", "AMP": "red", "MUT": "#008000"}
    af_tcga = {
        "background": alter_graphic("rect", fill="#CCCCCC"),
        "HOMDEL": alter_graphic("rect", fill=col_tcga["HOMDEL"]),
        "AMP": alter_graphic("rect", fill=col_tcga["AMP"]),
        "MUT": alter_graphic("rect", height=0.33, fill=col_tcga["MUT"]),
    }
    column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
    hlp_tcga = {
        "title": "Alternations",
        "at": ["HOMDEL", "AMP", "MUT"],
        "labels": ["Deep deletion", "Amplification", "Mutation"],
    }

    reset_state()
    op_tcga = oncoPrint(
        tcga, alter_fun=af_tcga, col=col_tcga,
        column_title=column_title, heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8), pct_gp=gpar(fontsize=8),
    )
    save_ht(op_tcga, "ch07_15_tcga_default.png", width=12.0, height=8.0)

    reset_state()
    op_tcga_ne = oncoPrint(
        tcga, alter_fun=af_tcga, col=col_tcga,
        remove_empty_columns=True, remove_empty_rows=True,
        column_title=column_title, heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8), pct_gp=gpar(fontsize=8),
    )
    save_ht(op_tcga_ne, "ch07_16_tcga_remove_empty.png", width=12.0, height=7.6)

    # Reorder with cBioPortal sample_order (aligned to non-empty rows/columns)
    reset_state()
    non_empty_rows = op_tcga_ne._rownames
    non_empty_cols_set = set(op_tcga_ne._colnames)
    filtered_sample_order = [s for s in sample_order if s in non_empty_cols_set]
    op_tcga_cbio = oncoPrint(
        tcga,
        alter_fun=af_tcga,
        col=col_tcga,
        row_order=non_empty_rows,
        column_order=filtered_sample_order,
        remove_empty_columns=True,
        remove_empty_rows=True,
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_tcga_cbio, "ch07_17_tcga_cbioportal_order.png", width=12.0, height=7.6)

    # Subset alteration types in top/right anno_oncoprint_barplot
    reset_state()
    op_sub_bp = oncoPrint(
        tcga,
        alter_fun=af_tcga,
        col=col_tcga,
        top_annotation=HeatmapAnnotation(
            column_barplot=anno_oncoprint_barplot(["MUT"], border=True, height=unit(4, "cm")),
        ),
        right_annotation=rowAnnotation(
            row_barplot=anno_oncoprint_barplot(
                ["AMP", "HOMDEL"],
                which="row",
                border=True,
                width=unit(4, "cm"),
                axis_param={"side": "bottom", "labels_rot": 90},
            ),
        ),
        remove_empty_columns=True,
        remove_empty_rows=True,
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_sub_bp, "ch07_18_tcga_subset_barplot.png", width=12.0, height=8.2)

    # Fraction barplots
    reset_state()
    op_frac = oncoPrint(
        tcga,
        alter_fun=af_tcga,
        col=col_tcga,
        top_annotation=HeatmapAnnotation(column_barplot=anno_oncoprint_barplot(show_fraction=True)),
        right_annotation=rowAnnotation(row_barplot=anno_oncoprint_barplot(which="row", show_fraction=True)),
        remove_empty_columns=True,
        remove_empty_rows=True,
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_frac, "ch07_19_tcga_fraction_barplot.png", width=12.0, height=7.8)

    # pct_side="right", row_names_side="left"
    reset_state()
    op_sides = oncoPrint(
        tcga,
        alter_fun=af_tcga,
        col=col_tcga,
        remove_empty_columns=True,
        remove_empty_rows=True,
        pct_side="right",
        row_names_side="left",
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_sides, "ch07_20_tcga_pct_right_names_left.png", width=12.0, height=7.6)

    # Custom annotations alongside OncoPrint barplots (pre-filtered non-empty submatrix)
    reset_state()
    tcga_ne = tcga.loc[non_empty_rows, [c for c in tcga.columns if c in non_empty_cols_set]]
    ng_ne, ns_ne = tcga_ne.shape
    op_extra = oncoPrint(
        tcga_ne,
        alter_fun=af_tcga,
        col=col_tcga,
        top_annotation=HeatmapAnnotation(
            cbar=anno_oncoprint_barplot(),
            foo1=np.arange(1, ns_ne + 1),
            bar1=anno_points(np.arange(1, ns_ne + 1)),
        ),
        left_annotation=rowAnnotation(foo2=np.arange(1, ng_ne + 1)),
        right_annotation=rowAnnotation(bar2=anno_barplot(np.arange(1, ng_ne + 1))),
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_extra, "ch07_21_tcga_extra_annotations.png", width=12.5, height=8.2, padding=unit([8, 4, 4, 4], "mm"))

    # oncoPrint as a Heatmap: horizontal, vertical, and split
    rng_expr = np.random.default_rng(123)
    reset_state()
    ht_horiz = (
        oncoPrint(
            tcga, alter_fun=af_tcga, col=col_tcga,
            column_title=column_title, heatmap_legend_param=hlp_tcga,
            row_names_gp=gpar(fontsize=8), pct_gp=gpar(fontsize=8),
        )
        + Heatmap(rng_expr.standard_normal((tcga.shape[0], 10)), name="expr", width=unit(4, "cm"))
    )
    save_ht(ht_horiz, "ch07_22_tcga_concat_horiz.png", width=12.5, height=8.0)

    reset_state()
    ht_vert = (
        oncoPrint(
            tcga, alter_fun=af_tcga, col=col_tcga,
            column_title=column_title, heatmap_legend_param=hlp_tcga,
            row_names_gp=gpar(fontsize=8), pct_gp=gpar(fontsize=8),
        )
        % v
        % Heatmap(rng_expr.standard_normal((10, tcga.shape[1])), name="expr", height=unit(4, "cm"))
    )
    save_ht(ht_vert, "ch07_23_tcga_concat_vert.png", width=12.0, height=9.0)

    reset_state()
    row_grp = rng_expr.choice(["Pathway A", "Pathway B"], tcga.shape[0], replace=True)
    save_ht(ht_horiz, "ch07_24_tcga_split.png", width=12.5, height=8.2, row_split=row_grp)

    reset_state()
    kras_status = ["KRAS altered" if str(x).strip() != "" else "KRAS WT" for x in tcga_ne.loc["KRAS"]]
    gene_grp = ["RAS/RAF" if g in {"KRAS", "HRAS", "BRAF", "RAF1"} else "MAPK cascade" for g in tcga_ne.index]
    op_rc_split = oncoPrint(
        tcga_ne,
        alter_fun=af_tcga,
        col=col_tcga,
        row_split=gene_grp,
        column_split=kras_status,
        column_title=column_title,
        heatmap_legend_param=hlp_tcga,
        row_names_gp=gpar(fontsize=8),
        pct_gp=gpar(fontsize=8),
    )
    save_ht(op_rc_split, "ch07_25_tcga_row_col_split.png", width=12.5, height=8.2)

    print("Filtered genes:", list(op_tcga_ne._rownames))
    print("Filtered sample count:", len(op_tcga_ne._colnames))
    print("Generated 25 figures for Chapter 7.")


if __name__ == "__main__":
    main()
