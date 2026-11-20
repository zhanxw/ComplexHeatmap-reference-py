"""Generate all figures for Chapter 10: Integrate with other packages."""

from __future__ import annotations

import numpy as np
import pandas as pd

from _common import IMAGES_DIR, reset_state, save_Standalone, save_ht
import complexheatmap as ch


def _build_pheatmap_test_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Create the classic pheatmap example matrix and annotations."""
    rng = np.random.default_rng(123)
    arr = rng.standard_normal((20, 10))
    arr[:10, ::2] += 3.0
    arr[10:, 1::2] += 2.0
    arr[14:, 1::2] += 4.0

    test = pd.DataFrame(
        arr,
        index=[f"Gene{i + 1}" for i in range(20)],
        columns=[f"Test{j + 1}" for j in range(10)],
    )
    annotation_col = pd.DataFrame(
        {
            "CellType": ["CT1", "CT2"] * 5,
            "Time": [1, 2, 3, 4, 5] * 2,
        },
        index=test.columns,
    )
    annotation_row = pd.DataFrame(
        {
            "GeneClass": ["Path1"] * 10 + ["Path2"] * 4 + ["Path3"] * 6,
        },
        index=test.index,
    )
    ann_colors = {
        "Time": ["white", "firebrick"],
        "CellType": {"CT1": "#1B9E77", "CT2": "#D95F02"},
        "GeneClass": {"Path1": "#7570B3", "Path2": "#E7298A", "Path3": "#66A61E"},
    }
    return test, annotation_col, annotation_row, ann_colors


def render_pheatmap_examples() -> None:
    """Render Section 10.1 pheatmap, heatmap, and heatmap.2 figures."""
    test, annotation_col, annotation_row, ann_colors = _build_pheatmap_test_data()

    save_ht(ch.pheatmap(test), "ch10_pheatmap_basic.png", width=5.6, height=6.2)
    reset_state()

    save_ht(
        ch.pheatmap(
            test,
            annotation_col=annotation_col,
            annotation_row=annotation_row,
            annotation_colors=ann_colors,
        ),
        "ch10_pheatmap_anno.png",
        width=6.6,
        height=6.7,
        padding=ch.unit([18, 4, 4, 4], "mm"),
    )
    reset_state()

    p1 = ch.pheatmap(test, name="mat1")
    p2 = ch.rowAnnotation(foo=ch.anno_barplot(np.arange(1, len(test) + 1), width=ch.unit(2, "cm")))
    p3 = ch.pheatmap(
        test,
        name="mat2",
        color=ch.hcl_colors(50, "RdBu", rev=True),
    )
    save_ht(p1 + p2 + p3, "ch10_pheatmap_concat.png", width=9.6, height=6.2)
    reset_state()

    p_orig = ch.pheatmap(
        test,
        annotation_col=annotation_col,
        annotation_row=annotation_row,
        annotation_colors=ann_colors,
    )
    save_ht(
        p_orig,
        "ch10_pheatmap_legend_grouping.png",
        width=6.8,
        height=6.7,
        padding=ch.unit([18, 4, 4, 4], "mm"),
        legend_grouping="original",
    )
    reset_state()

    # 2x2 grid of pheatmap comparison variants
    star_mat = np.where(test.to_numpy() > 4.5, "*", "")
    g1 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_pheatmap(test, scale="row", clustering_distance_rows="correlation", main="Scaled + Pearson"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=4.8,
    )
    g2 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_pheatmap(test, display_numbers=True, number_format="{:.1f}", fontsize_number=7, main="Display Numbers"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=4.8,
    )
    g3 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_pheatmap(test, annotation_col=annotation_col, cluster_rows=False, gaps_row=[10, 14], cutree_cols=2, main="Gaps + Cutree"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=4.8,
    )
    g4 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_pheatmap(test, display_numbers=star_mat, angle_col=45, annotation_colors=ann_colors, annotation_col=annotation_col, main="Star Mask + Angle 45"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=4.8,
    )

    def _draw_variants() -> None:
        panels = [(0.0, 1.0, g1), (0.5, 1.0, g2), (0.0, 0.5, g3), (0.5, 0.5, g4)]
        for x0, y1, gr in panels:
            ch.pushViewport(ch.viewport(x=x0, y=y1, width=0.5, height=0.5, just=("left", "top")))
            ch.grid_draw(gr)
            ch.popViewport()

    save_Standalone(_draw_variants, "ch10_pheatmap_variants.png", width=10.6, height=9.8)
    reset_state()

    # Base stats::heatmap and gplots::heatmap.2 translations
    row_side = ["#1B9E77"] * 10 + ["#D95F02"] * 4 + ["#7570B3"] * 6
    col_side = ["#E7298A", "#66A61E"] * 5
    gh1 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_heatmap(test, RowSideColors=row_side, ColSideColors=col_side, main="stats::heatmap() port"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=5.2,
    )
    gh2 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.compare_heatmap_2(test, scale="row", dendrogram="row", RowSideColors=row_side, main="gplots::heatmap.2() port"),
            newpage=False,
        ),
        width_in=5.2,
        height_in=5.2,
    )

    def _draw_legacy() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.48, height=1.0, just="left"))
        ch.grid_draw(gh1)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, height=1.0, just="left"))
        ch.grid_draw(gh2)
        ch.popViewport()

    save_Standalone(_draw_legacy, "ch10_heatmap_and_heatmap2.png", width=11.0, height=5.4)
    reset_state()


def render_multipanel_and_matplotlib() -> None:
    """Render Section 10.2 multi-panel grid_grabExpr and Matplotlib integration figures."""
    rng = np.random.default_rng(1234)
    # Construct a 60x60 similarity matrix with 3 clusters for clean anno_zoom demo
    sizes = [18, 22, 20]
    n_tot = sum(sizes)
    sim = rng.uniform(0.02, 0.18, size=(n_tot, n_tot))
    offset = 0
    cl = []
    for k, sz in enumerate(sizes, start=1):
        sim[offset : offset + sz, offset : offset + sz] = rng.uniform(0.55, 0.95, size=(sz, sz))
        cl.extend([f"C{k}"] * sz)
        offset += sz
    sim = 0.5 * (sim + sim.T)
    np.fill_diagonal(sim, 1.0)
    mat_sim = pd.DataFrame(sim, index=[f"Term_{i + 1}" for i in range(n_tot)], columns=[f"Term_{i + 1}" for i in range(n_tot)])

    term_boxes = {
        "C1": pd.DataFrame({"word": ["Cell cycle", "Mitosis", "DNA repair", "Spindle"], "col": ["#E41A1C", "#B2182B", "#D6604D", "#F4A582"], "fontsize": [9.5, 8.8, 8.2, 8]}),
        "C2": pd.DataFrame({"word": ["Immune response", "Cytokine", "T-cell", "Interferon"], "col": ["#377EB8", "#2166AC", "#4393C3", "#92C5DE"], "fontsize": [9.5, 8.8, 8.2, 8]}),
        "C3": pd.DataFrame({"word": ["Metabolism", "Oxidative", "Lipids", "ATP synth."], "col": ["#4DAF4A", "#1B7837", "#5AAE61", "#A6DBA0"], "fontsize": [9.5, 8.8, 8.2, 8]}),
    }

    def ht_clusters() -> None:
        ht = ch.Heatmap(
            mat_sim,
            name="Similarity",
            col=ch.colorRamp2([0.0, 0.5, 1.0], ["#FFFFFF", "#FCAE91", "#CB181D"]),
            cluster_rows=False,
            cluster_columns=False,
            show_row_names=False,
            show_column_names=False,
            column_title="Panel B: Similarity Heatmap with Matched-Size anno_textbox",
            column_title_gp=ch.gpar(fontsize=10.5, fontface="bold"),
            border=True,
            left_annotation=ch.rowAnnotation(
                Cluster=cl,
                col={"Cluster": {"C1": "#E41A1C", "C2": "#377EB8", "C3": "#4DAF4A"}},
                show_annotation_name=False,
            ),
            right_annotation=ch.rowAnnotation(
                Terms=ch.anno_textbox(cl, term_boxes, max_width=ch.unit(5.0, "cm")),
            ),
        )
        ch.draw(ht, newpage=False)

    # Target panel occupies 88% width and 78% height of a 10.4x5.2 inch canvas
    fig_w, fig_h = 10.4, 5.2
    p2_w = fig_w * 0.88
    p2_h = fig_h * 0.78
    g_hm = ch.grid_grabExpr(ht_clusters, width_in=p2_w, height_in=p2_h)

    def _draw_multipanel() -> None:
        # Left sidebar panel p1 (width = 0.12)
        ch.pushViewport(ch.viewport(x=0.0, y=0.5, width=0.12, height=0.96, just=("left", "center")))
        ch.grid_rect(width=0.86, height=0.96, gp=ch.gpar(fill="#F3F5F8", col="#607D8B", lwd=1.2))
        ch.grid_text("Panel A\nMetadata", rot=90, gp=ch.gpar(fontsize=10, fontface="bold", col="#37474F"))
        ch.popViewport()

        # Top-right heatmap panel p2 (width = 0.88, height = 0.78)
        ch.pushViewport(ch.viewport(x=0.12, y=1.0, width=0.88, height=0.78, just=("left", "top")))
        ch.grid_draw(g_hm)
        ch.popViewport()

        # Bottom-right summary panel p3 (width = 0.88, height = 0.22)
        ch.pushViewport(ch.viewport(x=0.12, y=0.0, width=0.88, height=0.22, just=("left", "bottom")))
        ch.grid_rect(width=0.96, height=0.82, gp=ch.gpar(fill="#F8F9FA", col="#607D8B", lwd=1.2))
        ch.grid_text("Panel C: Summary Track / Matched-Size grid_grabExpr Composition", gp=ch.gpar(fontsize=10, fontface="bold", col="#37474F"))
        ch.popViewport()

    save_Standalone(_draw_multipanel, "ch10_multipanel_grid.png", width=fig_w, height=fig_h)
    reset_state()

    # Matplotlib Axes + ComplexHeatmap on a shared GridDevice.fig
    rng_m = np.random.default_rng(77)
    nr, nc = 30, 12
    expr = rng_m.normal(0.0, 0.6, size=(nr, nc))
    expr[:15, :6] += 1.9
    expr[15:, 6:] += 1.9
    df_expr = pd.DataFrame(expr, index=[f"G{i + 1:02d}" for i in range(nr)], columns=[f"S{j + 1:02d}" for j in range(nc)])
    grp = ["Control"] * 6 + ["Treated"] * 6
    gene_cls = ["Up in Ctrl"] * 15 + ["Up in Treated"] * 15

    ht_left = ch.Heatmap(
        df_expr,
        name="Z-score",
        col=ch.colorRamp2([-2.5, 0.0, 2.5], ["#2166AC", "#F7F7F7", "#B2182B"]),
        row_split=gene_cls,
        top_annotation=ch.HeatmapAnnotation(
            Group=grp,
            col={"Group": {"Control": "#1B9E77", "Treated": "#D95F02"}},
            annotation_name_gp=ch.gpar(fontsize=8),
        ),
        show_row_names=False,
        column_names_gp=ch.gpar(fontsize=8),
        column_names_rot=45,
        column_title="A. Clustered Expression Heatmap",
        column_title_gp=ch.gpar(fontsize=11, fontface="bold"),
        border=True,
    )

    g_ht_left = ch.grid_grabExpr(lambda: ch.draw(ht_left, newpage=False), width_in=11.2 * 0.56, height_in=5.3)

    dev = ch.grid_newpage(width_in=11.2, height_in=5.5, dpi=140)
    try:
        ch.pushViewport(ch.viewport(x=0.0, y=0.5, width=0.56, height=0.96, just=("left", "center")))
        ch.grid_draw(g_ht_left)
        ch.popViewport()

        # Native Matplotlib Axes on the right half of dev.fig
        ax1 = dev.fig.add_axes([0.64, 0.57, 0.32, 0.34])
        log2fc = expr[:, 6:].mean(axis=1) - expr[:, :6].mean(axis=1)
        neg_log10_p = np.abs(log2fc) * 2.2 + rng_m.uniform(0.3, 1.5, size=nr)
        colors_pts = ["#1B9E77" if c == "Up in Ctrl" else "#D95F02" for c in gene_cls]
        ax1.scatter(log2fc, neg_log10_p, c=colors_pts, s=36, edgecolors="black", linewidths=0.5, alpha=0.88)
        ax1.axvline(0.0, color="#888888", linestyle="--", linewidth=0.9)
        ax1.axhline(1.3, color="#888888", linestyle=":", linewidth=0.9)
        ax1.set_title("B. Differential Expression Volcano", fontsize=10.5, fontweight="bold")
        ax1.set_xlabel("log2 Fold Change (Treated - Control)", fontsize=8.5)
        ax1.set_ylabel("-log10(p-value)", fontsize=8.5)
        ax1.tick_params(labelsize=8)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        ax2 = dev.fig.add_axes([0.63, 0.11, 0.33, 0.32])
        x_s = np.arange(nc)
        mean_ctrl = expr[:15, :].mean(axis=0)
        mean_trt = expr[15:, :].mean(axis=0)
        ax2.plot(x_s, mean_ctrl, marker="o", markersize=4.5, color="#1B9E77", linewidth=1.8, label="Up in Ctrl")
        ax2.plot(x_s, mean_trt, marker="s", markersize=4.5, color="#D95F02", linewidth=1.8, label="Up in Treated")
        ax2.axvspan(-0.5, 5.5, color="#1B9E77", alpha=0.08)
        ax2.axvspan(5.5, 11.5, color="#D95F02", alpha=0.08)
        ax2.set_xticks(x_s)
        ax2.set_xticklabels(df_expr.columns, rotation=45, fontsize=7.5)
        ax2.set_ylabel("Mean Z-score", fontsize=8.5)
        ax2.set_ylim(-0.5, 2.75)
        ax2.set_title("C. Module Mean Profile Across Samples", fontsize=10.5, fontweight="bold")
        ax2.legend(fontsize=7.5, frameon=False, loc="upper center", ncol=2)
        ax2.tick_params(labelsize=8)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        dev.fig.savefig(IMAGES_DIR / "ch10_matplotlib_integration.png", dpi=140, facecolor="white")
    finally:
        import matplotlib.pyplot as plt
        plt.close("all")
    reset_state()


def render_gt_render_examples() -> None:
    """Render Section 10.3 gt_render rich text examples."""
    rng = np.random.default_rng(123)
    mat = pd.DataFrame(rng.standard_normal((10, 10)), index=list("abcdefghij"), columns=[f"C{j + 1}" for j in range(10)])

    with ch.ht_opt.local(TITLE_PADDING=ch.unit([26, 26], "pt")):
        save_ht(
            ch.Heatmap(
                mat,
                name="mat",
                column_title=ch.gt_render(
                    "Some <span style='color:blue'>blue text **in bold.**</span><br>And *italics text.*<br>And some <span style='font-size:15pt; color:black'>large</span> text.",
                    r=ch.unit(2, "pt"),
                    padding=ch.unit([3, 5, 3, 5], "pt"),
                ),
                column_title_gp=ch.gpar(box_fill="orange"),
            ),
            "ch10_gt_render_title.png",
            width=5.8,
            height=5.8,
        )
    reset_state()

    save_ht(
        ch.Heatmap(
            mat,
            name="mat",
            row_km=2,
            row_title=[ch.gt_render("**title1**", padding=ch.unit([2, 4, 2, 4], "pt")), ch.gt_render("*title2*", padding=ch.unit([2, 4, 2, 4], "pt"))],
            row_title_gp=ch.gpar(box_fill=["yellow", "lightblue"]),
        ),
        "ch10_gt_render_split_title.png",
        width=5.6,
        height=5.2,
    )
    reset_state()

    save_ht(
        ch.Heatmap(
            mat,
            name="mat",
            row_labels=ch.gt_render(list("abcdefghij"), padding=ch.unit([2, 10, 2, 10], "pt")),
            row_names_gp=ch.gpar(box_col=["red", "green"] * 5),
            column_labels=ch.gt_render([f"**C{j + 1}**" if j % 2 == 0 else f"*C{j + 1}*" for j in range(10)]),
        ),
        "ch10_gt_render_dimnames.png",
        width=5.6,
        height=5.2,
    )
    reset_state()

    ha = ch.HeatmapAnnotation(
        foo=list("abcdefghij"),
        annotation_label="Annotation one",
        show_legend=False,
    )
    with ch.ht_opt.local(TITLE_PADDING=ch.unit([14, 14], "pt")):
        save_ht(
            ch.Heatmap(
                mat,
                name="mat",
                top_annotation=ha,
                column_title=ch.gt_render("**Rich Title** with <span style='color:#1B9E77'>**colored span**</span>", padding=ch.unit([2, 5, 2, 5], "pt")),
                column_title_gp=ch.gpar(box_fill="#F0F4F8", box_col="#377EB8"),
            ),
            "ch10_gt_render_anno.png",
            width=6.2,
            height=5.4,
            padding=ch.unit([4, 4, 4, 16], "mm"),
        )
    reset_state()

    save_ht(
        ch.rowAnnotation(
            foo=ch.anno_text(
                ch.gt_render([c * 8 for c in "ABCDEFGHIJ"], align_widths=True, padding=ch.unit([1, 3, 1, 3], "pt")),
                gp=ch.gpar(box_col="blue", box_lwd=2),
                just="right",
                location=ch.unit(1, "npc"),
            )
        )
        + ch.Heatmap(mat, name="mat"),
        "ch10_gt_render_anno_text.png",
        width=7.2,
        height=5.2,
    )
    reset_state()

    save_ht(
        ch.Heatmap(
            mat,
            name="mat",
            heatmap_legend_param={
                "title": ch.gt_render("<span style='color:orange'>**Legend title**</span>", padding=ch.unit([2, 4, 2, 4], "pt")),
                "title_gp": ch.gpar(box_fill="grey"),
                "at": [-3, 0, 3],
                "labels": ["negative three", "zero", "positive three"],
            },
        ),
        "ch10_gt_render_legend.png",
        width=6.0,
        height=5.0,
    )
    reset_state()

    # Additional rich text word-box annotation (anno_textbox) aligned to contiguous groups
    mat_big = rng.standard_normal((30, 8))
    groups = ["Pathway A"] * 10 + ["Pathway B"] * 10 + ["Pathway C"] * 10
    tb_data = {
        "Pathway A": pd.DataFrame({"word": ["Cell Cycle", "Mitosis", "G2/M Check", "CDK1"], "col": ["#B2182B", "#D6604D", "#E41A1C", "#333333"], "fontsize": [11, 9.5, 9, 8.5]}),
        "Pathway B": pd.DataFrame({"word": ["Immune", "Cytokine", "NF-kB", "IL-6"], "col": ["#2166AC", "#4393C3", "#377EB8", "#333333"], "fontsize": [11, 9.5, 9, 8.5]}),
        "Pathway C": pd.DataFrame({"word": ["Glycolysis", "TCA Cycle", "ATP", "OXPHOS"], "col": ["#1B7837", "#4DAF4A", "#5AAE61", "#333333"], "fontsize": [11, 9.5, 9, 8.5]}),
    }
    with ch.ht_opt.local(TITLE_PADDING=ch.unit([10, 10], "pt")):
        save_ht(
            ch.Heatmap(
                mat_big,
                name="Expr",
                cluster_rows=False,
                show_row_names=False,
                column_title=ch.gt_render("**Pathway Word Boxes** via *anno_textbox()*"),
                left_annotation=ch.rowAnnotation(
                    Pathway=groups,
                    col={"Pathway": {"Pathway A": "#B2182B", "Pathway B": "#2166AC", "Pathway C": "#1B7837"}},
                    show_annotation_name=False,
                ),
                right_annotation=ch.rowAnnotation(Summary=ch.anno_textbox(groups, tb_data)),
            ),
            "ch10_gt_render_textbox.png",
            width=7.6,
            height=5.6,
            padding=ch.unit([4, 14, 4, 4], "mm"),
        )
    reset_state()


def main() -> None:
    reset_state()
    render_pheatmap_examples()
    render_multipanel_and_matplotlib()
    render_gt_render_examples()
    print(f"Chapter 10 figures written to {IMAGES_DIR}")


if __name__ == "__main__":
    main()
