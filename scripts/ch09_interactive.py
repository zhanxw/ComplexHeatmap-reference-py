"""Generate figures for Chapter 9: Interactive ComplexHeatmap."""

from __future__ import annotations

import numpy as np
import pandas as pd

from _common import IMAGES_DIR, reset_state, save_Standalone
import complexheatmap as ch


def _build_demo_matrix() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create a 24x18 expression matrix with planted co-expressed blocks."""
    rng = np.random.default_rng(2026)
    nr, nc = 24, 18
    mat = rng.normal(0.0, 0.55, size=(nr, nc))
    # Block 1: high in first 8 samples for first 9 genes
    mat[:9, :8] += 2.2
    # Block 2: low in middle samples for middle genes
    mat[8:16, 6:13] -= 2.0
    # Block 3: high in last 6 samples for bottom genes
    mat[15:, 11:] += 2.3

    df = pd.DataFrame(
        mat,
        index=[f"Gene_{i + 1:02d}" for i in range(nr)],
        columns=[f"Sample_{j + 1:02d}" for j in range(nc)],
    )
    col_meta = pd.DataFrame(
        {
            "Condition": ["Ctrl"] * 6 + ["TreatA"] * 6 + ["TreatB"] * 6,
            "Score": np.round(rng.uniform(0.2, 1.0, size=nc), 2),
        },
        index=df.columns,
    )
    row_meta = pd.DataFrame(
        {"Pathway": ["CellCycle"] * 9 + ["Metabolism"] * 7 + ["Immune"] * 8},
        index=df.index,
    )
    return df, col_meta, row_meta


def render_interactive_zoom() -> None:
    """Render side-by-side original heatmap with brushed region and zoomed sub-heatmap."""
    df, col_meta, row_meta = _build_demo_matrix()
    col_fun = ch.colorRamp2([-2.5, 0.0, 2.5], ["#2166AC", "#F7F7F7", "#B2182B"])

    top_ha = ch.HeatmapAnnotation(
        Condition=col_meta["Condition"].tolist(),
        Score=ch.anno_barplot(col_meta["Score"].to_numpy(), gp=ch.gpar(fill="#5B8DBE"), height=ch.unit(10, "mm")),
        col={"Condition": {"Ctrl": "#66C2A5", "TreatA": "#FC8D62", "TreatB": "#8DA0CB"}},
        annotation_name_gp=ch.gpar(fontsize=8),
    )
    left_ha = ch.rowAnnotation(
        Pathway=row_meta["Pathway"].tolist(),
        col={"Pathway": {"CellCycle": "#E41A1C", "Metabolism": "#377EB8", "Immune": "#4DAF4A"}},
        show_annotation_name=False,
    )

    ht = ch.Heatmap(
        df,
        name="Expr",
        col=col_fun,
        top_annotation=top_ha,
        left_annotation=left_ha,
        show_row_names=True,
        show_column_names=True,
        row_names_gp=ch.gpar(fontsize=6.5),
        column_names_gp=ch.gpar(fontsize=6.5),
        column_names_rot=45,
        column_title="Original Heatmap (Brushed Region)",
        column_title_gp=ch.gpar(fontsize=11, fontface="bold"),
        border=True,
    )

    # Step 1: Draw offscreen / prepare a layout copy to extract post-clustering order
    res = ch.draw(ht, show=False)
    r_ord = ch.row_order(res)
    c_ord = ch.column_order(res)

    # Step 2: Brush display rows [r0, r1) and display cols [c0, c1)
    r0, r1 = 2, 10   # 8 display rows from top
    c0, c1 = 1, 7    # 6 display columns from left
    r_sel = r_ord[r0:r1]
    c_sel = c_ord[c0:c1]

    # Step 3: Slice original un-initialized Heatmap (preserving annotations!) and lock display order
    ht_sub = ht[r_sel, c_sel]
    sub_df = df.iloc[r_sel, c_sel]

    def _cell_fun(j: int, i: int, x, y, w, h, fill) -> None:
        val = float(sub_df.iloc[i, j])
        txt_col = "white" if abs(val) > 1.4 else "#222222"
        ch.grid_text(f"{val:.2f}", x, y, gp=ch.gpar(fontsize=7.5, col=txt_col))

    ht_zoom = ch.Heatmap(
        sub_df,
        name="ZoomExpr",
        col=col_fun,
        top_annotation=ht_sub.top_annotation,
        left_annotation=ht_sub.left_annotation,
        cluster_rows=False,
        cluster_columns=False,
        show_row_names=True,
        show_column_names=True,
        row_names_gp=ch.gpar(fontsize=8.5),
        column_names_gp=ch.gpar(fontsize=8.5),
        column_names_rot=45,
        cell_fun=_cell_fun,
        rect_gp=ch.gpar(col="#FFFFFF", lwd=0.8),
        border=True,
        show_heatmap_legend=True,
        heatmap_legend_param={"title": "Expr"},
        column_title=f"Zoomed Sub-heatmap ({len(r_sel)} x {len(c_sel)})",
        column_title_gp=ch.gpar(fontsize=11, fontface="bold"),
    )

    nr_tot, nc_tot = df.shape

    def _draw_left() -> None:
        ch.draw(ht, newpage=False, merge_legends=True)
        with ch.decorate_heatmap_body("Expr"):
            # Native coordinates: x in [0, nc_tot], y in [0, nr_tot] bottom-to-top
            x_center = 0.5 * (c0 + c1)
            y_center = nr_tot - 0.5 * (r0 + r1)
            bw = float(c1 - c0)
            bh = float(r1 - r0)
            ch.grid_rect(
                x=x_center,
                y=y_center,
                width=bw,
                height=bh,
                default_units="native",
                gp=ch.gpar(fill="#FFD70030", col="#FFD700", lwd=2.6, lty="solid"),
            )
            ch.grid_rect(
                x=x_center,
                y=y_center,
                width=bw,
                height=bh,
                default_units="native",
                gp=ch.gpar(fill="transparent", col="#111111", lwd=1.0, lty="dashed"),
            )
            ch.grid_rect(
                x=c1 + 0.25,
                y=nr_tot - r0 - 0.1,
                width=6.3,
                height=1.1,
                just=["left", "top"],
                default_units="native",
                gp=ch.gpar(fill="#FFFFFFE6", col="#333333", lwd=0.8),
            )
            ch.grid_text(
                f"Brush [{r0}:{r1}, {c0}:{c1}]",
                x=c1 + 0.55,
                y=nr_tot - r0 - 0.65,
                just=["left", "center"],
                default_units="native",
                gp=ch.gpar(fontsize=7.5, fontface="bold", col="#111111"),
            )

    def _draw_right() -> None:
        ch.draw(ht_zoom, newpage=False, merge_legends=True)

    g_left = ch.grid_grabExpr(_draw_left, width_in=5.6, height_in=5.4)
    g_right = ch.grid_grabExpr(_draw_right, width_in=5.8, height_in=5.4)

    def _compose() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.47, height=1.0, just="left"))
        ch.grid_draw(g_left)
        ch.popViewport()

        # Connector cue in center gutter
        ch.pushViewport(ch.viewport(x=0.49, width=0.04, height=1.0, just="left"))
        ch.grid_segments(
            ch.unit(0.1, "npc"),
            ch.unit(0.55, "npc"),
            ch.unit(0.85, "npc"),
            ch.unit(0.55, "npc"),
            gp=ch.gpar(col="#444444", lwd=2.0),
        )
        ch.grid_polygon(
            ch.unit([0.85, 0.65, 0.65], "npc"),
            ch.unit([0.55, 0.57, 0.53], "npc"),
            gp=ch.gpar(fill="#444444", col="#444444"),
        )
        ch.popViewport()

        ch.pushViewport(ch.viewport(x=0.52, width=0.48, height=1.0, just="left"))
        ch.grid_draw(g_right)
        ch.popViewport()

    save_Standalone(_compose, "ch09_interactive_zoom.png", width=11.6, height=5.5)
    reset_state()


def render_split_slice_zoom() -> None:
    """Render brushing/extracting a slice from a k-means split heatmap."""
    df, col_meta, _ = _build_demo_matrix()
    col_fun = ch.colorRamp2([-2.5, 0.0, 2.5], ["#2166AC", "#F7F7F7", "#B2182B"])

    ht_split = ch.Heatmap(
        df,
        name="SplitExpr",
        col=col_fun,
        row_km=2,
        column_km=2,
        row_title="Row Cluster %s",
        column_title="Col Cluster %s",
        top_annotation=ch.HeatmapAnnotation(
            Condition=col_meta["Condition"].tolist(),
            col={"Condition": {"Ctrl": "#66C2A5", "TreatA": "#FC8D62", "TreatB": "#8DA0CB"}},
        ),
        show_row_names=True,
        show_column_names=True,
        row_names_gp=ch.gpar(fontsize=6.5),
        column_names_gp=ch.gpar(fontsize=6.5),
        column_names_rot=45,
        border=True,
    )
    res = ch.draw(ht_split, show=False)
    r_slices = ch.row_order(res)
    c_slices = ch.column_order(res)

    # Pick slice (row_slice=1, col_slice=1) -> 0-based indices r_slices[0], c_slices[0]
    r_sel = r_slices[0]
    c_sel = c_slices[0]
    sub_df = df.iloc[r_sel, c_sel]

    ht_sub_base = ht_split[r_sel, c_sel]
    ht_slice = ch.Heatmap(
        sub_df,
        name="SliceExpr",
        col=col_fun,
        top_annotation=ht_sub_base.top_annotation,
        cluster_rows=False,
        cluster_columns=False,
        show_row_names=True,
        show_column_names=True,
        row_names_gp=ch.gpar(fontsize=8),
        column_names_gp=ch.gpar(fontsize=8),
        column_names_rot=45,
        border=True,
        column_title=f"Extracted Slice (1, 1): {len(r_sel)} x {len(c_sel)}",
        column_title_gp=ch.gpar(fontsize=11, fontface="bold"),
    )

    def _draw_left() -> None:
        ch.draw(ht_split, newpage=False)
        with ch.decorate_heatmap_body("SplitExpr", row_slice=1, column_slice=1):
            ch.grid_rect(
                gp=ch.gpar(fill="#FFD70022", col="#FF9800", lwd=3.0),
            )

    def _draw_right() -> None:
        ch.draw(ht_slice, newpage=False)

    g_l = ch.grid_grabExpr(_draw_left, width_in=5.5, height_in=5.0)
    g_r = ch.grid_grabExpr(_draw_right, width_in=5.5, height_in=5.0)

    def _compose() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.49, just="left"))
        ch.grid_draw(g_l)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, just="left"))
        ch.grid_draw(g_r)
        ch.popViewport()

    save_Standalone(_compose, "ch09_split_slice_zoom.png", width=11.2, height=5.2)
    reset_state()


def main() -> None:
    reset_state()
    render_interactive_zoom()
    render_split_slice_zoom()
    print(f"Chapter 9 figures written to {IMAGES_DIR}")


if __name__ == "__main__":
    main()
