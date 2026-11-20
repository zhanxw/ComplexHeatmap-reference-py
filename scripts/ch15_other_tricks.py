"""Generate figures for Chapter 15: Other Tricks."""

from __future__ import annotations

import numpy as np
import pandas as pd

import complexheatmap as ch
from _common import reset_state, save_Standalone, save_ht


def _random_mat(nr: int, seed: int = 123) -> pd.DataFrame:
    rng = np.random.default_rng(seed + nr)
    m = rng.normal(size=(nr, 10))
    return pd.DataFrame(m, index=[f"G{i + 1}" for i in range(nr)], columns=list("abcdefghij"))


def _make_fixed_cell_heatmap(nr: int, title: str | None = None) -> ch.Heatmap:
    m = _random_mat(nr)
    return ch.Heatmap(
        m,
        name="mat",
        height=ch.unit(5, "mm") * nr,
        width=ch.unit(5, "mm") * 10,
        column_title=title or f"{nr} rows, each row is 5 mm",
        top_annotation=ch.HeatmapAnnotation(bar=np.arange(1, 11)),
    )


def render_same_cell_size() -> tuple[float, float]:
    """Section 15.1: Fit exact linear relation between row count and figure height."""
    reset_state()
    xs = np.array([10.0, 20.0])
    heights_in = []
    widths_in = []
    for nr in (10, 20):
        ht = _make_fixed_cell_heatmap(int(nr), title="foo")
        w_in, h_in = ch.ht_size(ht, "inches", valueOnly=True)
        widths_in.append(w_in)
        heights_in.append(h_in)

    slope, intercept = np.polyfit(xs, heights_in, 1)
    # Add 4 mm (2 mm top + 2 mm bottom padding) when drawing onto PNG device
    pad_mm = 4.0
    pad_in = pad_mm / 25.4
    dev_w = max(widths_in) + pad_in

    for nr, fname in ((10, "ch15_cell_size_10rows.png"), (20, "ch15_cell_size_20rows.png")):
        reset_state()
        ht = _make_fixed_cell_heatmap(nr, title=f"{nr} rows, each row is 5mm")
        fig_h = float(slope * nr + intercept + pad_in)
        save_ht(ht, fname, width=float(dev_w), height=fig_h, padding=ch.unit([2, 2, 2, 2], "mm"))

    return float(slope), float(intercept)


def render_ht_opt_demo() -> None:
    """Section 15.2: Global and scoped styling via ht_opt and ht_opt.local."""
    reset_state()
    rng = np.random.default_rng(42)
    m = pd.DataFrame(rng.normal(size=(10, 8)), index=[f"Gene{i + 1}" for i in range(10)], columns=[f"S{j + 1}" for j in range(8)])
    ha = ch.HeatmapAnnotation(
        Group=["Ctrl"] * 4 + ["Treat"] * 4,
        Score=ch.anno_barplot(rng.uniform(0.2, 1.0, 8)),
        col={"Group": {"Ctrl": "#66C2A5", "Treat": "#FC8D62"}},
    )

    # Scoped global options applied uniformly across heatmaps & annotations
    with ch.ht_opt.local(
        heatmap_border=True,
        annotation_border=True,
        heatmap_row_names_gp=ch.gpar(fontsize=9, fontface="italic", col="#222222"),
        heatmap_column_names_gp=ch.gpar(fontsize=9, fontface="bold"),
        heatmap_column_title_gp=ch.gpar(fontsize=12, fontface="bold", col="#1B4F72"),
        legend_title_gp=ch.gpar(fontsize=9, fontface="bold"),
        TITLE_PADDING=ch.unit([6, 2], "mm"),
    ):
        ht = ch.Heatmap(
            m,
            name="Z-score",
            top_annotation=ha,
            row_km=2,
            column_split=["Ctrl"] * 4 + ["Treat"] * 4,
            column_title="Styled with ht_opt.local(...)",
        )
        save_ht(ht, "ch15_ht_opt_scoped.png", width=7.5, height=5.8)

    # Verify reset restores defaults
    ch.ht_opt(RESET=True)


def render_grid_grabexpr_composition() -> None:
    """Section 15.3: Capture and compose independent heatmaps with grid_grabExpr and Viewport."""
    reset_state()
    rng = np.random.default_rng(7)
    m1 = pd.DataFrame(rng.normal(size=(12, 8)), index=[f"R{i + 1}" for i in range(12)], columns=[f"A{j + 1}" for j in range(8)])
    m2 = pd.DataFrame(rng.uniform(0, 10, size=(8, 6)), index=[f"P{i + 1}" for i in range(8)], columns=[f"B{j + 1}" for j in range(6)])

    ht_a = ch.Heatmap(
        m1,
        name="Expr",
        col=ch.colorRamp2([-2, 0, 2], ["#2166AC", "white", "#B2182B"]),
        column_title="Panel A: Z-score Expression",
        row_km=2,
        row_names_gp=ch.gpar(fontsize=8),
        column_names_gp=ch.gpar(fontsize=8),
    )
    ht_b = ch.Heatmap3D(
        m2,
        name="Score",
        col=ch.colorRamp2([0, 5, 10], ["#FFFFCC", "#FD8D3C", "#BD0026"]),
        column_title="Panel B: 3D Pathway Activity",
        top_annotation=ch.HeatmapAnnotation(_spacer=ch.anno_empty(border=False, height=ch.unit(6, "mm")), show_annotation_name=False),
        bar_max_length=ch.unit(7, "mm"),
        row_names_gp=ch.gpar(fontsize=8),
        column_names_gp=ch.gpar(fontsize=8),
    )

    grob_a = ch.grid_grabExpr(
        lambda: ch.draw(ht_a, padding=ch.unit([3, 3, 3, 3], "mm"), newpage=False),
        width_in=5.5,
        height_in=5.0,
    )
    grob_b = ch.grid_grabExpr(
        lambda: ch.draw(ht_b, padding=ch.unit([3, 3, 5, 5], "mm"), newpage=False),
        width_in=5.5,
        height_in=5.0,
    )

    def _compose() -> None:
        ch.grid_rect(gp=ch.gpar(fill="white", col=None))
        ch.pushViewport(ch.Viewport(x=ch.unit(0.0, "npc"), y=ch.unit(0.0, "npc"), width=ch.unit(0.5, "npc"), height=ch.unit(1.0, "npc"), just=("left", "bottom")))
        ch.grid_draw(grob_a)
        ch.grid_text("A", x=ch.unit(4, "mm"), y=ch.unit(1, "npc") - ch.unit(4, "mm"), just=("left", "top"), gp=ch.gpar(fontsize=14, fontface="bold"))
        ch.popViewport()

        ch.pushViewport(ch.Viewport(x=ch.unit(0.5, "npc"), y=ch.unit(0.0, "npc"), width=ch.unit(0.5, "npc"), height=ch.unit(1.0, "npc"), just=("left", "bottom")))
        ch.grid_draw(grob_b)
        ch.grid_text("B", x=ch.unit(4, "mm"), y=ch.unit(1, "npc") - ch.unit(4, "mm"), just=("left", "top"), gp=ch.gpar(fontsize=14, fontface="bold"))
        ch.popViewport()

    save_Standalone(_compose, "ch15_grid_grab_compose.png", width=11.0, height=5.0)


def main() -> None:
    slope, intercept = render_same_cell_size()
    print(f"Linear model: height_in = {slope:.4f} * nr + {intercept:.4f}")
    render_ht_opt_demo()
    render_grid_grabexpr_composition()
    reset_state()
    print("Chapter 15 figures generated successfully.")


if __name__ == "__main__":
    main()
