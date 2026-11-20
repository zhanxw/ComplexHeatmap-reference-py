"""Generate all figures for Chapter 6: Heatmap Decoration."""

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
    anno_barplot,
    anno_boxplot,
    anno_empty,
    anno_points,
    colorRamp2,
    column_dend,
    cutree,
    decorate_annotation,
    decorate_column_dend,
    decorate_column_names,
    decorate_column_title,
    decorate_dend,
    decorate_dimnames,
    decorate_heatmap_body,
    decorate_row_dend,
    decorate_row_names,
    decorate_row_title,
    decorate_title,
    generateRandomBed,
    gpar,
    grid,
    list_components,
    order_dendrogram,
    popViewport,
    pushViewport,
    rowAnnotation,
    row_order,
    seekViewport,
    unit,
    viewport,
)

from _common import reset_state, save_decorated


def _build_demo_heatmap_list() -> ch.HeatmapList:
    rng = np.random.default_rng(123)
    mat = np.vstack([rng.normal(2.0, 1.0, (8, 10)), rng.normal(-2.0, 1.0, (4, 10))])
    df = pd.DataFrame(
        mat,
        index=[f"R{i}" for i in range(1, 13)],
        columns=[f"C{j}" for j in range(1, 11)],
    )

    ha_column1 = HeatmapAnnotation(
        points=anno_points(rng.normal(0.0, 1.0, 10)),
        annotation_name_side="left",
    )
    ht1 = Heatmap(
        df,
        name="ht1",
        km=2,
        row_title=["Slice 1", "Slice 2"],
        column_title="Heatmap 1",
        top_annotation=ha_column1,
        row_names_side="left",
    )

    ha_column2 = HeatmapAnnotation(
        type=["a"] * 5 + ["b"] * 5,
        col={"type": {"a": "red", "b": "blue"}},
    )
    ht2 = Heatmap(
        df,
        name="ht2",
        row_title="Heatmap 2",
        column_title="Heatmap 2",
        bottom_annotation=ha_column2,
        column_km=2,
    )

    return ht1 + ht2 + rowAnnotation(bar=anno_barplot(df.mean(axis=1).to_numpy(), width=unit(2, "cm")))


def main() -> None:
    # ---- 6.1 Inspect and outline all component viewports ----
    reset_state()
    ht_list = _build_demo_heatmap_list()

    def _outline_components() -> None:
        comps = list_components()
        print("Components:", comps)
        for hm_name in ("ht1", "ht2"):
            for i in (1, 2):
                vp = f"{hm_name}_heatmap_body_{i}_1"
                if vp in comps:
                    seekViewport(vp)
                    grid.rect(gp=gpar(fill=None, col="red"))
                if hm_name == "ht1":
                    seekViewport(f"{hm_name}_row_names_{i}")
                    grid.rect(gp=gpar(fill=None, col="red"))
            seekViewport(f"{hm_name}_column_title_1")
            grid.rect(gp=gpar(fill=None, col="red"))
            seekViewport(f"{hm_name}_dend_column_1")
            grid.rect(gp=gpar(fill=None, col="red"))
            seekViewport(f"{hm_name}_column_names_1")
            grid.rect(gp=gpar(fill=None, col="red"))

        for vp in (
            "ht2_heatmap_body_1_2",
            "ht2_heatmap_body_2_2",
            "ht2_dend_column_2",
            "ht2_column_names_2",
            "ht1_dend_row_1",
            "ht1_dend_row_2",
            "ht1_row_title_1",
            "ht1_row_title_2",
            "annotation_points_1",
            "annotation_type_1",
            "annotation_type_2",
            "annotation_bar_1",
            "annotation_bar_2",
        ):
            seekViewport(vp)
            grid.rect(gp=gpar(fill=None, col="red"))

    save_decorated(
        ht_list,
        "ch06_01_access_components.png",
        _outline_components,
        width=10.5,
        height=7.2,
        row_title="Heatmap list",
        column_title="Heatmap list",
    )

    # ---- 6.2 Decoration functions showcase ----
    reset_state()
    ht_list2 = _build_demo_heatmap_list()

    def _decorate_showcase(res: ch.HeatmapList) -> None:
        with decorate_heatmap_body("ht1", slice=2):
            grid.text("outlier", 1.5 / 10.0, 2.5 / 4.0, default_units="npc")
            grid.lines([0.5, 0.5], [0.0, 1.0], gp=gpar(lty=2, lwd=2))

        with decorate_column_dend("ht1"):
            tree = column_dend(res)["ht1"]
            cl = cutree(tree, k=2)
            od = order_dendrogram(tree)
            ind = cl[od]
            n_col = len(ind)
            x1 = [int(np.where(ind == k_id)[0][0]) for k_id in (1, 2)]
            x2 = [int(np.where(ind == k_id)[0][-1]) + 1 for k_id in (1, 2)]
            grid.rect(
                x=[v / n_col for v in x1],
                width=[(b - a) / n_col for a, b in zip(x1, x2)],
                just="left",
                default_units="npc",
                gp=gpar(fill=["#FF000040", "#00FF0040"], col=None),
            )

        with decorate_row_names("ht1", slice=2):
            grid.rect(gp=gpar(fill="#FF000040"))

        with decorate_row_title("ht1", slice=1):
            grid.rect(gp=gpar(fill="#00FF0040"))

        with decorate_annotation("points"):
            grid.lines([0, 1], unit([0, 0], "native"), gp=gpar(col="red"))

    save_decorated(
        ht_list2,
        "ch06_02_decorate_functions.png",
        _decorate_showcase,
        width=10.5,
        height=7.2,
        row_title="Heatmap list",
        column_title="Heatmap list",
        heatmap_legend_side="right",
        annotation_legend_side="left",
    )

    # ---- 6.3 2D slice targeting & all decorate_* variants ----
    reset_state()
    rng_s = np.random.default_rng(11)
    mat_s = pd.DataFrame(
        rng_s.standard_normal((10, 8)),
        index=[f"G{i}" for i in range(1, 11)],
        columns=[f"S{j}" for j in range(1, 9)],
    )
    ht_s = Heatmap(
        mat_s,
        name="M",
        row_split=["R1"] * 5 + ["R2"] * 5,
        column_split=["C1"] * 4 + ["C2"] * 4,
        row_title=["R1", "R2"],
        column_title=["C1", "C2"],
        top_annotation=HeatmapAnnotation(grp=["A", "B"] * 4),
    )

    def _decorate_slices() -> None:
        with decorate_heatmap_body("M", row_slice=1, column_slice=2):
            grid.rect(gp=gpar(fill="#FF000030", col="red", lwd=2))
            grid.text("slice (1, 2)", gp=gpar(col="darkred", fontsize=9, fontface="bold"))
        with decorate_heatmap_body("M", row_slice=2, column_slice=1):
            grid.circle(r=unit(0.35, "snpc"), gp=gpar(fill="#0000FF30", col="blue", lwd=2))
        with decorate_row_dend("M", slice=2):
            grid.rect(gp=gpar(col="blue", lwd=1.5, fill=None))
        with decorate_dend("M", which="column", slice=1):
            grid.rect(gp=gpar(col="purple", lwd=1.5, fill=None))
        with decorate_column_names("M", slice=2):
            grid.rect(gp=gpar(fill="#00AA0030", col=None))
        with decorate_dimnames("M", which="row", slice=1):
            grid.rect(gp=gpar(fill="#FFAA0030", col=None))
        with decorate_column_title("M", slice=1):
            grid.rect(gp=gpar(col="brown", lwd=1.5, fill=None))
        with decorate_title("M", which="row", slice=2):
            grid.rect(gp=gpar(col="orange", lwd=1.5, fill=None))
        with decorate_annotation("grp", slice=1):
            grid.circle(x=unit([0.25, 0.75], "npc"), gp=gpar(fill="#FF000080"))

    save_decorated(ht_s, "ch06_03_slice_targeting.png", _decorate_slices, width=7.5, height=5.8)

    # ---- 6.4 DMR single-column heatmap with summary barplot annotation ----
    reset_state()
    rng_d = np.random.default_rng(123)
    bed = generateRandomBed(nr=1000, nc=0, seed=123)
    frac = np.concatenate([rng_d.uniform(0.4, 1.0, 400), rng_d.uniform(0.2, 0.6, len(bed) - 400)])
    col_fun = colorRamp2([0, 1], ["white", "orange"])
    split = ["group1"] * 400 + ["group2"] * (len(bed) - 400)
    ht_dmr = Heatmap(
        frac,
        name="fraction",
        col=col_fun,
        width=unit(2.2, "cm"),
        top_annotation=HeatmapAnnotation(barplot=anno_empty(height=unit(4, "cm"))),
    )
    w = (bed["end"] - bed["start"]).to_numpy()

    def _decorate_dmr(res: ch.HeatmapList) -> None:
        ro = row_order(res)
        p = [float(np.sum(w[idx] * frac[idx]) / np.sum(w[idx])) for idx in ro]
        with decorate_annotation("barplot"):
            pushViewport(viewport(xscale=(0.5, 2.5), yscale=(0.0, max(p) * 1.1)))
            grid.rect(
                x=[1, 2],
                y=0,
                width=0.8,
                height=p,
                just="bottom",
                gp=gpar(fill="orange"),
                default_units="native",
            )
            grid.yaxis()
            grid.text("mean fraction", x=unit(-1.5, "cm"), rot=90, just="bottom")
            popViewport()

    save_decorated(
        ht_dmr,
        "ch06_04_dmr_barplot.png",
        _decorate_dmr,
        width=6.2,
        height=6.2,
        row_split=split,
        padding=unit([2, 16, 2, 2], "mm"),
    )

    # ---- 6.5 Titles for row annotations ----
    reset_state()
    rng_t = np.random.default_rng(42)
    ht_titles = (
        Heatmap(rng_t.standard_normal((10, 10)), name="mat", show_column_dend=False)
        + rowAnnotation(foo=anno_barplot(np.arange(1, 11), width=unit(4, "cm")))
        + rowAnnotation(bar=anno_boxplot(rng_t.standard_normal((10, 10)), width=unit(4, "cm")))
    )

    def _add_titles() -> None:
        with decorate_annotation("foo"):
            grid.text("title for barplot", y=unit(1, "npc") + unit(2, "mm"), just="bottom")
        with decorate_annotation("bar"):
            grid.text("title for boxplot", y=unit(1, "npc") + unit(2, "mm"), just="bottom")
        with decorate_heatmap_body("mat"):
            grid.text("title for the heatmap", y=unit(1, "npc") + unit(2, "mm"), just="bottom")

    save_decorated(
        ht_titles,
        "ch06_05_row_anno_titles.png",
        _add_titles,
        width=8.6,
        height=5.5,
        padding=unit([10, 2, 10, 2], "mm"),
    )

    print("Generated 5 figures for Chapter 6.")


if __name__ == "__main__":
    main()
