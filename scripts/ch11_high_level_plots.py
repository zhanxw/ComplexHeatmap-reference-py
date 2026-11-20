"""Generate figures for Chapter 11: Other High-level Plots."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import complexheatmap as ch
from _common import reset_state, save_decorated, save_ht


def build_density_matrix(seed: int = 123) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    m1 = rng.normal(size=(100, 10))
    m2 = rng.uniform(-2.0, 2.0, size=(100, 10)) + 0.5
    m = np.hstack([m1, m2])
    cols = [f"C{i + 1}" for i in range(m.shape[1])]
    return pd.DataFrame(m, columns=cols)


def render_density_examples() -> None:
    m = build_density_matrix(123)

    # 1. Default densityHeatmap (leave room on left for y-axis & right for quantile labels)
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    dh1 = ch.densityHeatmap(m)
    save_ht(dh1, "ch11_density_default.png", width=7.5, height=5.2, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 2. Custom ylim, title, ylab
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    dh2 = ch.densityHeatmap(m, ylim=(-2, 2), title="Distribution as heatmap", ylab="some values")
    save_ht(dh2, "ch11_density_ylim_title.png", width=7.5, height=5.2, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 3. Custom column_order (0-based permutation)
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    rng = np.random.default_rng(42)
    perm = rng.permutation(m.shape[1]).tolist()
    dh3 = ch.densityHeatmap(m, column_order=perm, title="Reordered columns")
    save_ht(dh3, "ch11_density_col_order.png", width=7.5, height=5.2, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 4. Custom color palette
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    dh4 = ch.densityHeatmap(m, col=ch.hcl_colors(10, palette="viridis"), title="Custom color palette")
    save_ht(dh4, "ch11_density_colors.png", width=7.5, height=5.2, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 5. Clustering columns by KS distance
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    dh5 = ch.densityHeatmap(m, cluster_columns=True, clustering_distance_columns="ks", title="KS-distance column clustering")
    save_ht(dh5, "ch11_density_cluster_ks.png", width=7.5, height=5.5, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 6. Dict of unequal-length distributions
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    rng_l = np.random.default_rng(99)
    dist_dict = {
        "Normal(0, 1)": rng_l.normal(0.0, 1.0, 120),
        "Normal(1.2, 0.6)": rng_l.normal(1.2, 0.6, 80),
        "Uniform(-2, 2)": rng_l.uniform(-2.0, 2.0, 200),
        "Bimodal": np.concatenate([rng_l.normal(-1.4, 0.45, 70), rng_l.normal(1.4, 0.45, 70)]),
        "Gamma-shifted": rng_l.gamma(2.0, 0.6, 150) - 1.2,
    }
    dh6 = ch.densityHeatmap(dist_dict, cluster_columns=True, title="Unequal-length distributions", ylab="signal")
    save_ht(dh6, "ch11_density_list.png", width=7.2, height=5.2, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 7. With top and bottom annotations
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    rng_a = np.random.default_rng(7)
    ha1 = ch.HeatmapAnnotation(
        dist=["rnorm"] * 10 + ["runif"] * 10,
        col={"dist": {"rnorm": "#E41A1C", "runif": "#377EB8"}},
    )
    ha2 = ch.HeatmapAnnotation(foo=ch.anno_points(rng_a.normal(size=20)))
    dh7 = ch.densityHeatmap(m, anno=ha1, bottom_annotation=ha2, title="Density heatmap with annotations")
    save_ht(dh7, "ch11_density_anno.png", width=8.0, height=5.8, padding=ch.unit([4, 18, 4, 4], "mm"), merge_legends=True)

    # 8. Vertical concatenation with HeatmapAnnotation and Heatmap
    reset_state()
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(14, "mm")
    rng_v = np.random.default_rng(21)
    mat_v = pd.DataFrame(rng_v.normal(size=(20, 20)), columns=m.columns)
    ht_v = (
        ch.densityHeatmap(m, title="Density + Heatmap vertical stack")
        % ch.v %
        ch.HeatmapAnnotation(foo=ch.anno_barplot(np.arange(1, 21)))
        % ch.v %
        ch.Heatmap(mat_v, name="mat", height=ch.unit(6, "cm"), show_column_names=True)
    )
    save_ht(ht_v, "ch11_density_concat.png", width=8.2, height=8.8, padding=ch.unit([4, 18, 4, 4], "mm"))

    # 9. frequencyHeatmap (2D and 3D)
    reset_state()
    fh = ch.frequencyHeatmap(m, title="Frequency heatmap (count)")
    save_ht(fh, "ch11_frequency_heatmap.png", width=7.5, height=5.2, padding=ch.unit([4, 10, 4, 4], "mm"))

    reset_state()
    ch.ht_opt.TITLE_PADDING = ch.unit([10, 2], "mm")
    ch.ht_opt.HEATMAP_LEGEND_PADDING = ch.unit(8, "mm")
    fh3d = ch.frequencyHeatmap(
        m,
        use_3d=True,
        title="3D Frequency heatmap",
        top_annotation=ch.HeatmapAnnotation(_spacer=ch.anno_empty(border=False, height=ch.unit(6, "mm")), show_annotation_name=False),
        bar_max_length=ch.unit(8, "mm"),
    )
    save_ht(fh3d, "ch11_frequency_heatmap_3d.png", width=8.2, height=6.0, padding=ch.unit([4, 10, 6, 6], "mm"))


def render_dmr_summary() -> None:
    reset_state()
    dmr = ch.load_dataset("dmr_summary")
    label = dmr["label"].index.tolist()
    n_gr = dmr["n_gr"].index.to_numpy(dtype=float)
    mean_meth = dmr["mean_meth"]
    n_corr = dmr["n_corr"]
    dist_tss = dmr["dist_tss"]
    gene_anno = dmr["gene_anno"]
    cgi_anno = dmr["cgi_anno"]
    mat_enrich_gf = dmr["mat_enrich_gf"]
    mat_pct_st = dmr["mat_pct_st"]
    mat_enrich_st = dmr["mat_enrich_st"]

    meth_col_fun = ch.colorRamp2([0, 0.5, 1], ["blue", "white", "red"])
    corr_col = ["green", "red"]
    dist_tss_col = ["#FF0000", "#FF7352", "#FFB299", "#FFD9CB"]
    gene_anno_col = ["green", "blue"]
    cgi_anno_col = ["#FFA500", "#FFD191"]
    z_score_col_fun = ch.colorRamp2([-200, 0, 200], ["green", "white", "red"])
    state_col = ["#FF0000", "#008000", "#C2E105", "#8A91D0", "#CD5C5C", "#808080", "#000000"]

    anno_width = ch.unit(2.8, "cm")
    ht_list = (
        ch.rowAnnotation(
            text=ch.anno_text(label, location=ch.unit(1, "npc"), just="right", gp=ch.gpar(fontsize=11)),
        )
        + ch.Heatmap(
            mean_meth,
            name="mean_meth",
            col=meth_col_fun,
            cluster_rows=False,
            row_title=None,
            cluster_columns=False,
            show_row_names=False,
            heatmap_legend_param={"title": "Methylation"},
            width=ch.unit(mean_meth.shape[1] * 5, "mm"),
        )
        + ch.rowAnnotation(n_gr=ch.anno_barplot(n_gr, bar_width=1, width=anno_width), show_annotation_name=False)
        + ch.rowAnnotation(n_corr=ch.anno_barplot(n_corr.to_numpy(), bar_width=1, gp=ch.gpar(fill=corr_col), width=anno_width), show_annotation_name=False)
        + ch.rowAnnotation(dist_tss=ch.anno_barplot(dist_tss.to_numpy(), bar_width=1, gp=ch.gpar(fill=dist_tss_col), width=anno_width), show_annotation_name=False)
        + ch.rowAnnotation(gene_anno=ch.anno_barplot(gene_anno.to_numpy(), bar_width=1, gp=ch.gpar(fill=gene_anno_col), width=anno_width), show_annotation_name=False)
        + ch.rowAnnotation(cgi_anno=ch.anno_barplot(cgi_anno.to_numpy(), bar_width=1, gp=ch.gpar(fill=cgi_anno_col), width=anno_width), show_annotation_name=False)
        + ch.Heatmap(
            mat_enrich_gf,
            name="enrich_gf",
            col=z_score_col_fun,
            cluster_columns=False,
            width=ch.unit(mat_enrich_gf.shape[1] * 4.5, "mm"),
            column_title="",
            heatmap_legend_param={"title": "Z-score"},
        )
        + ch.rowAnnotation(pct_st=ch.anno_barplot(mat_pct_st.to_numpy(), bar_width=1, gp=ch.gpar(fill=state_col), width=anno_width), show_annotation_name=False)
        + ch.Heatmap(
            mat_enrich_st,
            name="enrich_st",
            col=z_score_col_fun,
            cluster_columns=False,
            width=ch.unit(mat_enrich_st.shape[1] * 5.5, "mm"),
            column_title="",
            show_heatmap_legend=False,
            column_names_gp=ch.gpar(col=state_col, fontsize=9),
            show_row_names=False,
        )
    )

    lgd_list = [
        ch.Legend(labels=["neg", "pos"], title="Corr. genes", legend_gp=ch.gpar(fill=corr_col)),
        ch.Legend(labels=["gene", "intergenic"], title="Gene annotation", legend_gp=ch.gpar(fill=gene_anno_col)),
        ch.Legend(labels=["<1kb", "1kb~5kb", "5kb~10kb", ">10kb"], title="Distance to TSS", legend_gp=ch.gpar(fill=dist_tss_col)),
        ch.Legend(labels=["CGI", "CGI shore"], title="CGI annotation", legend_gp=ch.gpar(fill=cgi_anno_col)),
        ch.Legend(labels=list(mat_enrich_st.columns), title="Chromatin states", legend_gp=ch.gpar(fill=state_col)),
    ]

    split_grp = [re.sub(r"\d+$", "", s) for s in label]

    def _add_titles() -> None:
        anno_titles = {
            "n_gr": "Number of\nDMRs",
            "n_corr": "Significantly\ncorrelated genes",
            "dist_tss": "Distance to TSS",
            "gene_anno": "Gene annotation",
            "cgi_anno": "CGI annotation",
            "pct_st": "Overlap to\nChromatin states",
        }
        for an, ttl in anno_titles.items():
            with ch.decorate_annotation(an, slice=1):
                ch.grid_text(ttl, y=ch.unit(1, "npc") + ch.unit(3, "mm"), just="bottom", gp=ch.gpar(fontsize=9.5))

        ht_titles = {
            "mean_meth": "Mean\nmethylation",
            "enrich_gf": "Enrichment to\ngenomic features",
            "enrich_st": "Enrichment to\nchromatin states",
        }
        for hn, ttl in ht_titles.items():
            with ch.decorate_heatmap_body(hn, slice=1):
                ch.grid_text(ttl, y=ch.unit(1, "npc") + ch.unit(3, "mm"), just="bottom", gp=ch.gpar(fontsize=9.5))

    save_decorated(
        ht_list,
        "ch11_dmr_summary.png",
        _add_titles,
        width=15.0,
        height=6.2,
        padding=ch.unit([4, 4, 18, 4], "mm"),
        row_split=split_grp,
        heatmap_legend_list=lgd_list,
        heatmap_legend_side="bottom",
    )


def render_vertical_annotations() -> None:
    reset_state()
    rng = np.random.default_rng(123)
    prop = np.concatenate([
        rng.uniform(0.1, 0.5, 10),
        rng.uniform(0.2, 0.4, 10),
        rng.uniform(0.3, 0.6, 10),
        rng.uniform(0.4, 0.8, 10),
    ])
    median_length = np.concatenate([
        rng.uniform(5000, 10000, 10),
        rng.uniform(6000, 20000, 10),
        rng.uniform(7000, 15000, 10),
        rng.uniform(6000, 30000, 10),
    ])
    group = np.repeat(["a", "b", "c", "d"], 10).tolist()

    ht_list = (
        ch.HeatmapAnnotation(
            prop=ch.anno_barplot(
                prop,
                height=ch.unit(3.5, "cm"),
                axis_param={"at": [0, 0.2, 0.4, 0.6, 0.8], "labels": ["0%", "20%", "40%", "60%", "80%"]},
            ),
            annotation_name_rot=90,
        )
        % ch.v %
        ch.HeatmapAnnotation(
            median_length=ch.anno_barplot(
                median_length,
                height=ch.unit(3.5, "cm"),
                axis_param={"at": [0, 10000, 20000], "labels": ["0kb", "10kb", "20kb"]},
            ),
            annotation_name_rot=90,
        )
        % ch.v %
        ch.HeatmapAnnotation(
            group=group,
            annotation_legend_param={"group": {"legend_direction": "horizontal", "nrow": 1, "title_position": "leftcenter"}},
        )
    )
    save_ht(
        ht_list,
        "ch11_stacked_anno_vertical.png",
        width=8.5,
        height=5.5,
        padding=ch.unit([4, 14, 4, 18], "mm"),
        heatmap_legend_side="bottom",
        annotation_legend_side="bottom",
        column_title="Statistics for a list of genomic regions",
    )


def main() -> None:
    render_density_examples()
    render_dmr_summary()
    render_vertical_annotations()
    reset_state()
    print("Chapter 11 figures generated successfully.")


if __name__ == "__main__":
    main()
