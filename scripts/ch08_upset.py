"""Generate all figures for Chapter 8: UpSet plot."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from _common import IMAGES_DIR, reset_state, save_Standalone, save_decorated, save_ht
import complexheatmap as ch
from complexheatmap import v


def _make_sample_sets() -> dict[str, list[str]]:
    """Create reproducible 3-set dictionary with all 7 distinct combinations populated."""
    rng = np.random.default_rng(2)
    letters = list("abcdefghijklmnopqrstuvwxyz")
    return {
        "a": sorted(rng.choice(letters, 6, replace=False).tolist()),
        "b": sorted(rng.choice(letters, 10, replace=False).tolist()),
        "c": sorted(rng.choice(letters, 15, replace=False).tolist()),
    }


def _build_movies_df() -> pd.DataFrame:
    """Construct deterministic synthetic movies dataset mirroring UpSetR's movies schema."""
    rng = np.random.default_rng(2026)
    n = 3883

    # Decade distribution peaking strongly in the 1990s (1919-2000)
    decades = np.array([1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000])
    dec_probs = np.array([0.008, 0.015, 0.032, 0.045, 0.060, 0.080, 0.110, 0.175, 0.425, 0.050])
    dec_probs /= dec_probs.sum()
    chosen_dec = rng.choice(decades, size=n, p=dec_probs)
    years = np.clip(chosen_dec + rng.integers(0, 10, size=n), 1919, 2000)

    genres = ["Action", "Romance", "Horror", "Children", "SciFi", "Documentary", "Comedy", "Drama"]
    mat = np.zeros((n, len(genres)), dtype=int)

    # Archetypal genre combinations matching MovieLens/UpSetR patterns
    archetypes = [
        (["Drama"], 0.24),
        (["Comedy"], 0.19),
        (["Comedy", "Drama"], 0.08),
        (["Comedy", "Romance"], 0.07),
        (["Drama", "Romance"], 0.06),
        (["Action"], 0.06),
        (["Horror"], 0.055),
        (["Action", "SciFi"], 0.045),
        (["Children", "Comedy"], 0.04),
        (["Documentary"], 0.035),
        (["Action", "Drama"], 0.035),
        (["SciFi"], 0.025),
        (["Children"], 0.02),
        (["Action", "Horror"], 0.015),
        (["Children", "SciFi"], 0.015),
        (["Action", "Children"], 0.010),
    ]
    names_a = [a[0] for a in archetypes]
    probs_a = np.array([a[1] for a in archetypes])
    probs_a /= probs_a.sum()

    picked = rng.choice(len(archetypes), size=n, p=probs_a)
    g_idx = {g: i for i, g in enumerate(genres)}
    for r_i, a_i in enumerate(picked):
        for g in names_a[a_i]:
            mat[r_i, g_idx[g]] = 1

    # Add realistic older skew for SciFi+Children & rating patterns across genres
    scifi_child = (mat[:, g_idx["SciFi"]] == 1) & (mat[:, g_idx["Children"]] == 1)
    years[scifi_child] = rng.integers(1965, 1986, size=int(scifi_child.sum()))

    base_rating = rng.normal(3.25, 0.72, size=n)
    base_rating[mat[:, g_idx["Horror"]] == 1] -= 0.75
    base_rating[mat[:, g_idx["Romance"]] == 1] += 0.35
    base_rating[mat[:, g_idx["Documentary"]] == 1] += 0.55
    base_rating[mat[:, g_idx["Drama"]] == 1] += 0.20
    act_child = (mat[:, g_idx["Action"]] == 1) & (mat[:, g_idx["Children"]] == 1)
    base_rating[act_child] -= 0.65
    avg_rating = np.round(np.clip(base_rating, 0.4, 4.95), 2)

    watches = rng.integers(30, 600, size=n)
    watches[mat[:, g_idx["Action"]] == 1] += rng.integers(120, 900, size=int((mat[:, g_idx["Action"]] == 1).sum()))
    watches[mat[:, g_idx["Comedy"]] == 1] += rng.integers(60, 400, size=int((mat[:, g_idx["Comedy"]] == 1).sum()))
    watches[mat[:, g_idx["Documentary"]] == 1] = np.maximum(
        15, watches[mat[:, g_idx["Documentary"]] == 1] - 120
    )

    df = pd.DataFrame({"Name": [f"Movie_{i + 1}" for i in range(n)], "ReleaseDate": years})
    for gi, g in enumerate(genres):
        df[g] = mat[:, gi]
    df["AvgRating"] = avg_rating
    df["Watches"] = watches
    return df


def _build_genomic_peaks() -> dict[str, ch.GRanges]:
    """Generate 6 representative H3K4me3 peak sets (ESC, ES-derived, somatic tissues)."""
    rng = np.random.default_rng(77)
    chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"]
    chr_len = 120_000_000

    def _make_intervals(n_iv: int, min_w: int, max_w: int, anchor_starts: np.ndarray | None = None) -> tuple[list[str], np.ndarray, np.ndarray]:
        c_list = [chroms[i % len(chroms)] for i in range(n_iv)]
        if anchor_starts is None:
            starts = rng.integers(10_000, chr_len - max_w - 10_000, size=n_iv)
        else:
            starts = anchor_starts + rng.integers(-150, 150, size=len(anchor_starts))
            starts = np.clip(starts, 1_000, chr_len - max_w - 1_000)
        widths = rng.integers(min_w, max_w, size=n_iv)
        return c_list, starts, starts + widths

    # Core housekeeping promoters shared by all 6 samples (~14 Mb)
    c_core, s_core, e_core = _make_intervals(4200, 2600, 4200)
    # Stem-cell-shared promoters (~8 Mb)
    c_es, s_es, e_es = _make_intervals(2400, 2400, 4000)
    # Somatic-shared promoters (~7 Mb)
    c_som, s_som, e_som = _make_intervals(2100, 2400, 4000)

    specs = {
        "ESC": [(c_core, s_core, e_core), (c_es, s_es, e_es), _make_intervals(3200, 2200, 3800)],
        "ES-deriv1": [(c_core, s_core, e_core), (c_es, s_es, e_es), _make_intervals(2600, 2200, 3800)],
        "ES-deriv2": [(c_core, s_core, e_core), (c_es, s_es, e_es), _make_intervals(2800, 2200, 3800)],
        "Brain": [(c_core, s_core, e_core), (c_som, s_som, e_som), _make_intervals(3500, 2200, 3800)],
        "Muscle": [(c_core, s_core, e_core), (c_som, s_som, e_som), _make_intervals(2400, 2200, 3800)],
        "Heart": [(c_core, s_core, e_core), (c_som, s_som, e_som), _make_intervals(2500, 2200, 3800)],
    }

    peaks: dict[str, ch.GRanges] = {}
    for name, blocks in specs.items():
        all_c = [c for blk in blocks for c in blk[0]]
        all_s = np.concatenate([blk[1] for blk in blocks])
        all_e = np.concatenate([blk[2] for blk in blocks])
        peaks[name] = ch.GRanges(all_c, ch.IRanges(all_s, all_e))
    return peaks


def render_upset_mode_diagram() -> None:
    """Port of upset_mode.R illustrating distinct, intersect, and union modes."""
    r = 1.0
    h = 0.56

    def _circle_arc(cx: float, cy: float, frm: float = 0.0, to: float = 2.0 * math.pi) -> np.ndarray:
        theta = np.linspace(frm, to, 100)
        return np.column_stack([cx + r * np.cos(theta), cy + r * np.sin(theta)])

    def top_c(frm: float = 0.0, to: float = 2.0 * math.pi) -> np.ndarray:
        return _circle_arc(0.0, h, frm, to)

    def left_c(frm: float = 0.0, to: float = 2.0 * math.pi) -> np.ndarray:
        return _circle_arc(-h * math.sqrt(3.0) / 2.0, -h / 2.0, frm, to)

    def right_c(frm: float = 0.0, to: float = 2.0 * math.pi) -> np.ndarray:
        return _circle_arc(h * math.sqrt(3.0) / 2.0, -h / 2.0, frm, to)

    def draw_poly(pts: np.ndarray, fill: str, col: str) -> None:
        ch.grid_polygon(pts[:, 0], pts[:, 1], default_units="native", gp=ch.gpar(fill=fill, col=col))

    def draw_outline(pts: np.ndarray, col: str) -> None:
        ch.grid_lines(pts[:, 0], pts[:, 1], default_units="native", gp=ch.gpar(col=col, lwd=1.1))

    def make_venn(bits: tuple[int, int, int], mode: str, fill: str = "#111111", col: str = "#808080") -> None:
        top, left, right = bool(bits[0]), bool(bits[1]), bool(bits[2])
        xs = np.array([-1.3, 0.0, 1.3])
        if mode == "upset":
            ch.pushViewport(ch.viewport(xscale=(-2.0, 2.0), yscale=(-2.0, 2.0)))
        else:
            ch.pushViewport(ch.viewport(xscale=(-2.0, 2.0), yscale=(-2.0, 2.0), width=ch.unit(1, "snpc"), height=ch.unit(1, "snpc")))
        try:
            pi = math.pi
            if top and left and right:
                if mode in ("distinct", "intersect"):
                    draw_poly(np.vstack([top_c(-pi / 3, -2 * pi / 3), right_c(pi, 2 * pi / 3), left_c(pi / 3, 0)]), fill, col)
                elif mode == "union":
                    draw_poly(np.vstack([top_c(pi, 0), right_c(pi / 3, -2 * pi / 3), left_c(-pi / 3, -4 * pi / 3)]), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col="black"), default_units="native")
                    ch.grid_segments(xs[0], 0.0, xs[2], 0.0, gp=ch.gpar(lwd=2.5, col="black"), default_units="native")
            elif top and left and not right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(pi, 4 * pi / 3), right_c(pi, 2 * pi / 3), left_c(pi / 3, 2 * pi / 3)]), fill, col)
                elif mode == "intersect":
                    draw_poly(np.vstack([top_c(pi, 5 * pi / 3), left_c(0, 2 * pi / 3)]), fill, col)
                elif mode == "union":
                    draw_poly(np.vstack([top_c(-pi / 3, pi), left_c(2 * pi / 3, 2 * pi)]), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["black", "black", "#CCCCCC"]), default_units="native")
                    ch.grid_segments(xs[0], 0.0, xs[1], 0.0, gp=ch.gpar(lwd=2.5, col="black"), default_units="native")
            elif top and not left and right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(0, -pi / 3), left_c(0, pi / 3), right_c(2 * pi / 3, pi / 3)]), fill, col)
                elif mode == "intersect":
                    draw_poly(np.vstack([top_c(0, -2 * pi / 3), right_c(pi, pi / 3)]), fill, col)
                elif mode == "union":
                    draw_poly(np.vstack([top_c(4 * pi / 3, 0), right_c(pi / 3, -pi)]), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["black", "#CCCCCC", "black"]), default_units="native")
                    ch.grid_segments(xs[0], 0.0, xs[2], 0.0, gp=ch.gpar(lwd=2.5, col="black"), default_units="native")
            elif not top and left and right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(-pi / 3, -2 * pi / 3), right_c(pi, 4 * pi / 3), left_c(-pi / 3, 0)]), fill, col)
                elif mode == "intersect":
                    draw_poly(np.vstack([left_c(pi / 3, -pi / 3), right_c(-2 * pi / 3, -4 * pi / 3)]), fill, col)
                elif mode == "union":
                    draw_poly(np.vstack([left_c(pi / 3, 5 * pi / 3), right_c(-2 * pi / 3, 2 * pi / 3)]), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["#CCCCCC", "black", "black"]), default_units="native")
                    ch.grid_segments(xs[1], 0.0, xs[2], 0.0, gp=ch.gpar(lwd=2.5, col="black"), default_units="native")
            elif top and not left and not right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(0, pi), left_c(2 * pi / 3, pi / 3), right_c(2 * pi / 3, pi / 3)]), fill, col)
                elif mode in ("intersect", "union"):
                    draw_poly(top_c(), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["black", "#CCCCCC", "#CCCCCC"]), default_units="native")
            elif not top and left and not right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(pi, 4 * pi / 3), right_c(pi, 4 * pi / 3), left_c(-pi / 3, -4 * pi / 3)]), fill, col)
                elif mode in ("intersect", "union"):
                    draw_poly(left_c(), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["#CCCCCC", "black", "#CCCCCC"]), default_units="native")
            elif not top and not left and right:
                if mode == "distinct":
                    draw_poly(np.vstack([top_c(-pi / 3, 0), right_c(pi / 3, -2 * pi / 3), left_c(-pi / 3, 0)]), fill, col)
                elif mode in ("intersect", "union"):
                    draw_poly(right_c(), fill, col)
                else:
                    ch.grid_points(xs, np.zeros(3), pch=16, size=ch.unit(3.2, "mm"), gp=ch.gpar(col=["#CCCCCC", "#CCCCCC", "black"]), default_units="native")

            if mode != "upset":
                draw_outline(top_c(), col)
                draw_outline(left_c(), col)
                draw_outline(right_c(), col)
        finally:
            ch.popViewport()

    def _draw() -> None:
        combos = [
            (1, 1, 1),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ]
        modes = ["upset", "distinct", "intersect", "union"]
        header_h = 0.08
        row_h = (1.0 - header_h) / len(combos)
        col_w = 1.0 / 4.0

        # Column headers
        ch.pushViewport(ch.viewport(x=0.0, y=1.0 - header_h / 2.0, width=col_w, height=header_h, just="left", xscale=(-2.0, 2.0), yscale=(-1.0, 1.0)))
        ch.grid_text(["A", "B", "C"], x=np.array([-1.3, 0.0, 1.3]), y=np.zeros(3), default_units="native", gp=ch.gpar(fontsize=11, fontface="bold"))
        ch.popViewport()

        for ci in range(1, 4):
            ch.pushViewport(ch.viewport(x=ci * col_w, y=1.0 - header_h / 2.0, width=col_w, height=header_h, just="left"))
            ch.grid_text(modes[ci], gp=ch.gpar(fontsize=12, fontface="bold"))
            ch.popViewport()

        for ri, bits in enumerate(combos):
            yc = 1.0 - header_h - (ri + 0.5) * row_h
            for ci, md in enumerate(modes):
                xc = (ci + 0.5) * col_w
                ch.pushViewport(ch.viewport(x=xc, y=yc, width=col_w * 0.92, height=row_h * 0.92))
                make_venn(bits, mode=md, fill="#1B1B1B", col="#777777")
                ch.popViewport()

    save_Standalone(_draw, "ch08_upset_mode.png", width=5.4, height=6.8)


def render_basic_and_styling_plots() -> None:
    """Render Sections 8.4 and 8.5 figures."""
    lt = _make_sample_sets()
    letters = list("abcdefghijklmnopqrstuvwxyz")
    m = ch.make_comb_mat(lt)

    save_ht(ch.UpSet(m), "ch08_upset_basic.png", width=5.5, height=3.3)
    reset_state()

    save_ht(
        ch.UpSet(m, set_order=["a", "b", "c"], comb_order=np.argsort(ch.comb_size(m)).tolist()),
        "ch08_upset_order.png",
        width=5.5,
        height=3.3,
    )
    reset_state()

    pal = {1: "red", 2: "blue", 3: "black"}
    save_ht(
        ch.UpSet(
            m,
            pt_size=ch.unit(5, "mm"),
            lwd=3,
            comb_col=[pal[int(d)] for d in ch.comb_degree(m)],
        ),
        "ch08_upset_style.png",
        width=5.5,
        height=3.3,
    )
    reset_state()

    g_bg1 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", column_title="Single bg_col"),
            newpage=False,
            padding=ch.unit([8, 18, 4, 4], "mm"),
        ),
        width_in=5.5,
        height_in=3.4,
    )
    g_bg2 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.UpSet(m, comb_col="#0000FF", bg_col=["#F0F0FF", "#FFF0F0"], bg_pt_col="#CCCCFF", column_title="Alternating bg_col"),
            newpage=False,
            padding=ch.unit([8, 18, 4, 4], "mm"),
        ),
        width_in=5.5,
        height_in=3.4,
    )

    def _draw_bg() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.49, just="left"))
        ch.grid_draw(g_bg1)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, just="left"))
        ch.grid_draw(g_bg2)
        ch.popViewport()

    save_Standalone(_draw_bg, "ch08_upset_bg.png", width=11.2, height=3.4)
    reset_state()

    save_ht(ch.UpSet(ch.t(m)), "ch08_upset_transpose.png", width=3.8, height=5.2)
    reset_state()

    g_sub1 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.UpSet(m[:, ch.comb_size(m) >= 4], column_title="comb_size(m) >= 4"),
            newpage=False,
            padding=ch.unit([8, 18, 4, 4], "mm"),
        ),
        width_in=4.7,
        height_in=3.4,
    )
    g_sub2 = ch.grid_grabExpr(
        lambda: ch.draw(
            ch.UpSet(m[:, ch.comb_degree(m) == 2], column_title="comb_degree(m) == 2"),
            newpage=False,
            padding=ch.unit([8, 18, 4, 4], "mm"),
        ),
        width_in=4.7,
        height_in=3.4,
    )

    def _draw_sub() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.49, just="left"))
        ch.grid_draw(g_sub1)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, just="left"))
        ch.grid_draw(g_sub2)
        ch.popViewport()

    save_Standalone(_draw_sub, "ch08_upset_subset.png", width=9.6, height=3.4)
    reset_state()

    m1 = ch.make_comb_mat(lt, mode="distinct")
    m2 = ch.make_comb_mat(lt, mode="intersect")
    m3 = ch.make_comb_mat(lt, mode="union")

    g_modes = [
        ch.grid_grabExpr(
            lambda mat_i=mat_i, lbl=lbl: ch.draw(
                ch.UpSet(mat_i, column_title=lbl),
                newpage=False,
                padding=ch.unit([8, 18, 4, 4], "mm"),
            ),
            width_in=6.2,
            height_in=2.8,
        )
        for mat_i, lbl in [(m1, "mode = 'distinct'"), (m2, "mode = 'intersect'"), (m3, "mode = 'union'")]
    ]

    def _draw_modes() -> None:
        for idx, g_m in enumerate(g_modes):
            y_top = 1.0 - idx / 3.0
            ch.pushViewport(ch.viewport(y=y_top, height=1.0 / 3.0, just="top"))
            ch.grid_draw(g_m)
            ch.popViewport()

    save_Standalone(_draw_modes, "ch08_upset_modes.png", width=6.2, height=8.6)
    reset_state()

    m_univ = ch.make_comb_mat(lt, universal_set=letters)
    save_ht(ch.UpSet(m_univ), "ch08_upset_universal.png", width=5.8, height=3.3)
    reset_state()

    m_comp = ch.make_comb_mat(lt, complement_size=10)
    save_ht(ch.UpSet(m_comp), "ch08_upset_complement_size.png", width=5.8, height=3.3)
    reset_state()

    m_small = ch.make_comb_mat(lt, universal_set=letters[:10])
    save_ht(ch.UpSet(m_small), "ch08_upset_smaller_universal.png", width=5.6, height=3.3)
    reset_state()

    x_bm = ch.list_to_matrix(lt, universal_set=letters)
    m_rem = ch.make_comb_mat(x_bm)
    m_rem = m_rem[:, ch.comb_degree(m_rem) > 0]
    save_ht(ch.UpSet(m_rem), "ch08_upset_remove_complement.png", width=5.5, height=3.3)
    reset_state()

    # Section 8.5: UpSet plots as heatmaps
    save_ht(
        ch.UpSet(m, top_annotation=ch.upset_top_annotation(m, gp=ch.gpar(col=ch.comb_degree(m), fill="white"))),
        "ch08_upset_top_gp.png",
        width=5.5,
        height=3.4,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.upset_top_annotation(
                m,
                ylim=(0, 15),
                bar_width=1.0,
                axis_param={"side": "right", "at": [0, 5, 10, 15], "labels": ["zero", "five", "ten", "fifteen"]},
            ),
        ),
        "ch08_upset_top_axis.png",
        width=5.8,
        height=3.4,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.upset_top_annotation(
                m,
                annotation_name_rot=90,
                annotation_name_side="right",
                axis_param={"side": "right"},
            ),
        ),
        "ch08_upset_top_name.png",
        width=5.8,
        height=3.4,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            right_annotation=ch.upset_right_annotation(
                m,
                ylim=(0, 30),
                gp=ch.gpar(fill="green"),
                annotation_name_side="top",
                axis_param={"side": "top"},
            ),
        ),
        "ch08_upset_right_custom.png",
        width=5.6,
        height=3.6,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.HeatmapAnnotation(
                degree=[str(int(d)) for d in ch.comb_degree(m)],
                **{
                    "Intersection\nsize": ch.anno_barplot(
                        ch.comb_size(m), border=False, gp=ch.gpar(fill="black"), height=ch.unit(2, "cm")
                    )
                },
                annotation_name_side="left",
                annotation_name_rot=0,
            ),
        ),
        "ch08_upset_top_multi.png",
        width=6.0,
        height=3.8,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            right_annotation=ch.rowAnnotation(
                **{
                    "Set size": ch.anno_barplot(
                        ch.set_size(m), border=False, gp=ch.gpar(fill="black"), width=ch.unit(2, "cm")
                    )
                },
                group=["group1", "group1", "group2"],
            ),
        ),
        "ch08_upset_right_multi.png",
        width=6.0,
        height=3.4,
    )
    reset_state()

    save_ht(
        ch.UpSet(m, left_annotation=ch.upset_left_annotation(m), right_annotation=None),
        "ch08_upset_left.png",
        width=5.6,
        height=3.4,
    )
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.upset_top_annotation(m, add_numbers=True),
            right_annotation=ch.upset_right_annotation(m, add_numbers=True),
        ),
        "ch08_upset_numbers.png",
        width=5.8,
        height=3.6,
    )
    reset_state()

    ht = ch.UpSet(m)
    save_ht(
        ht + ch.Heatmap(np.arange(1, 4), name="foo", width=ch.unit(5, "mm")) + ch.rowAnnotation(bar=ch.anno_points(np.arange(1, 4))),
        "ch08_upset_concat_h.png",
        width=7.0,
        height=3.5,
    )
    reset_state()

    nc = m.ncol
    save_ht(
        ht
        % v
        % ch.Heatmap(np.arange(1, nc + 1).reshape(1, -1), name="foo", row_names_side="left", height=ch.unit(5, "mm"))
        % v
        % ch.HeatmapAnnotation(bar=ch.anno_points(np.arange(1, nc + 1)), annotation_name_side="left"),
        "ch08_upset_concat_v.png",
        width=6.0,
        height=4.8,
    )
    reset_state()

    m1n, m2n, m3n = ch.normalize_comb_mat(m1, m2, m3)

    save_ht(
        ch.UpSet(m1n, row_title="distinct mode") % v % ch.UpSet(m2n, row_title="intersect mode") % v % ch.UpSet(m3n, row_title="union mode"),
        "ch08_upset_stack_v.png",
        width=6.2,
        height=7.5,
        padding=ch.unit([4, 22, 4, 4], "mm"),
    )
    reset_state()

    save_ht(
        ch.UpSet(ch.t(m1n), column_title="distinct mode") + ch.UpSet(ch.t(m2n), column_title="intersect mode") + ch.UpSet(ch.t(m3n), column_title="union mode"),
        "ch08_upset_stack_h.png",
        width=7.8,
        height=5.0,
    )
    reset_state()

    # Align combination sizes across modes onto m1's combination names
    codes1 = ch.comb_name(m1)
    cs_dist = ch.comb_size(m1)
    cs_int = np.array([dict(zip(ch.comb_name(m2), ch.comb_size(m2))).get(c, 0.0) for c in codes1])
    cs_uni = np.array([dict(zip(ch.comb_name(m3), ch.comb_size(m3))).get(c, 0.0) for c in codes1])

    top_ha = ch.HeatmapAnnotation(
        distinct=ch.anno_barplot(cs_dist, gp=ch.gpar(fill="black"), height=ch.unit(1.8, "cm")),
        intersect=ch.anno_barplot(cs_int, gp=ch.gpar(fill="black"), height=ch.unit(1.8, "cm")),
        union=ch.anno_barplot(cs_uni, gp=ch.gpar(fill="black"), height=ch.unit(1.8, "cm")),
        gap=ch.unit(2, "mm"),
        annotation_name_side="left",
        annotation_name_rot=0,
    )
    save_ht(ch.UpSet(m1, top_annotation=top_ha), "ch08_upset_three_bars_top.png", width=6.0, height=6.4)
    reset_state()

    right_ha = ch.rowAnnotation(
        distinct=ch.anno_barplot(cs_dist, gp=ch.gpar(fill="black"), width=ch.unit(1.8, "cm")),
        intersect=ch.anno_barplot(cs_int, gp=ch.gpar(fill="black"), width=ch.unit(1.8, "cm")),
        union=ch.anno_barplot(cs_uni, gp=ch.gpar(fill="black"), width=ch.unit(1.8, "cm")),
        gap=ch.unit(2, "mm"),
        annotation_name_side="bottom",
    )
    save_ht(ch.UpSet(ch.t(m1), right_annotation=right_ha), "ch08_upset_three_bars_right.png", width=6.8, height=5.2)
    reset_state()


def render_movies_examples() -> None:
    """Render Section 8.6 movies dataset figures."""
    movies = _build_movies_df()
    all_genres = ["Action", "Romance", "Horror", "Children", "SciFi", "Documentary", "Comedy", "Drama"]

    m = ch.make_comb_mat(movies[all_genres], top_n_sets=6)
    m = m[:, ch.comb_degree(m) > 0]
    save_ht(ch.UpSet(m), "ch08_movies_basic.png", width=9.5, height=4.2)
    reset_state()

    ss = ch.set_size(m)
    cs = ch.comb_size(m)
    ht_fmt = ch.UpSet(
        m,
        set_order=np.argsort(ss).tolist(),
        comb_order=sorted(range(m.ncol), key=lambda j: (int(ch.comb_degree(m)[j]), -float(cs[j]))),
        top_annotation=ch.HeatmapAnnotation(
            **{
                "Genre Intersections": ch.anno_barplot(
                    cs,
                    ylim=(0, float(cs.max()) * 1.15),
                    border=False,
                    gp=ch.gpar(fill="black"),
                    height=ch.unit(4.0, "cm"),
                )
            },
            annotation_name_side="left",
            annotation_name_rot=90,
        ),
        left_annotation=ch.rowAnnotation(
            **{
                "Movies Per Genre": ch.anno_barplot(
                    -ss,
                    baseline=0.0,
                    axis_param={"at": [0, -500, -1000, -1500], "labels": ["0", "500", "1000", "1500"], "labels_rot": 0},
                    border=False,
                    gp=ch.gpar(fill="black"),
                    width=ch.unit(4.0, "cm"),
                ),
                "set_name": ch.anno_text(
                    ch.set_name(m),
                    location=0.5,
                    just="center",
                    width=ch.max_text_width(ch.set_name(m)) + ch.unit(4, "mm"),
                ),
            },
            annotation_name_side="top",
        ),
        right_annotation=None,
        show_row_names=False,
    )

    def _dec_movies(res) -> None:
        od = ch.column_order(res)
        with ch.decorate_annotation("Genre Intersections"):
            ch.grid_text(
                [str(int(cs[j])) for j in od],
                x=np.arange(1, len(cs) + 1, dtype=float),
                y=ch.unit(cs[od], "native") + ch.unit(2, "pt"),
                default_units="native",
                just=["left", "bottom"],
                gp=ch.gpar(fontsize=7, col="#404040"),
                rot=45,
            )

    save_decorated(
        ht_fmt,
        "ch08_movies_formatted.png",
        _dec_movies,
        width=10.2,
        height=4.8,
        padding=ch.unit([6, 6, 8, 4], "mm"),
    )
    reset_state()

    # Split movies by rating
    genre = ["Action", "Romance", "Horror", "Children", "SciFi", "Documentary"]
    bins = [0, 1, 2, 3, 4, 5]
    rating_cat = pd.cut(movies["AvgRating"], bins=bins)
    m_list = []
    rating_labels = []
    for cat in rating_cat.cat.categories:
        sub = movies.loc[rating_cat == cat, genre]
        mc = ch.make_comb_mat(sub)
        mc = mc[:, ch.comb_degree(mc) > 0]
        m_list.append(mc)
        rating_labels.append(str(cat))

    m_list = ch.normalize_comb_mat(m_list)
    tot_cs = np.sum([ch.comb_size(mk) for mk in m_list], axis=0)
    tot_deg = ch.comb_degree(m_list[0])
    shared_comb_ord = sorted(range(m_list[0].ncol), key=lambda j: (int(tot_deg[j]), -float(tot_cs[j])))
    shared_set_ord = list(range(m_list[0].nrow))

    max_set_size = float(max(ch.set_size(mk).max() for mk in m_list))
    max_comb_size = float(max(ch.comb_size(mk).max() for mk in m_list))

    ht_list = None
    for lbl, mk in zip(rating_labels, m_list):
        u = ch.UpSet(
            mk,
            row_title=f"rating in {lbl}",
            set_order=shared_set_ord,
            comb_order=shared_comb_ord,
            top_annotation=ch.upset_top_annotation(mk, ylim=(0, max_comb_size), height=ch.unit(1.8, "cm")),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_size)),
        )
        ht_list = ch.add_heatmap(ht_list, u, direction="vertical")

    save_ht(ht_list, "ch08_movies_by_rating.png", width=7.5, height=11.0)
    reset_state()

    # Relative fraction (excluding (0, 1])
    m_list_rel = m_list[1:]
    labels_rel = rating_labels[1:]
    max_set_rel = float(max(ch.set_size(mk).max() for mk in m_list_rel))

    ht_rel = None
    for lbl, mk in zip(labels_rel, m_list_rel):
        s = ch.comb_size(mk)
        rel = s / max(1.0, float(s.sum()))
        u = ch.UpSet(
            mk,
            row_title=f"rating in {lbl}",
            set_order=shared_set_ord,
            comb_order=shared_comb_ord,
            top_annotation=ch.HeatmapAnnotation(
                **{
                    "Relative\nfraction": ch.anno_barplot(
                        rel, ylim=(0, 0.5), gp=ch.gpar(fill="black"), border=False, height=ch.unit(1.8, "cm")
                    )
                },
                annotation_name_side="left",
                annotation_name_rot=0,
            ),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_rel)),
        )
        ht_rel = ch.add_heatmap(ht_rel, u, direction="vertical")

    save_ht(ht_rel, "ch08_movies_by_rating_rel.png", width=7.5, height=9.2)
    reset_state()

    # Split movies by decade
    decades_col = (movies["ReleaseDate"] // 10) * 10
    dec_vals = sorted(decades_col.unique().tolist())
    m_dec = []
    dec_labels = []
    for dv in dec_vals:
        sub = movies.loc[decades_col == dv, genre]
        md = ch.make_comb_mat(sub)
        md = md[:, ch.comb_degree(md) > 0]
        m_dec.append(md)
        dec_labels.append(f"{int(dv)}s")

    m_dec = ch.normalize_comb_mat(m_dec)
    tot_cs_d = np.sum([ch.comb_size(mk) for mk in m_dec], axis=0)
    tot_deg_d = ch.comb_degree(m_dec[0])
    dec_comb_ord = sorted(range(m_dec[0].ncol), key=lambda j: (int(tot_deg_d[j]), -float(tot_cs_d[j])))
    dec_set_ord = list(range(m_dec[0].nrow))
    max_set_d = float(max(ch.set_size(mk).max() for mk in m_dec))
    max_comb_d = float(max(ch.comb_size(mk).max() for mk in m_dec))

    ht_dec1 = None
    for lbl, mk in zip(dec_labels[:5], m_dec[:5]):
        u = ch.UpSet(
            mk,
            row_title=lbl,
            set_order=dec_set_ord,
            comb_order=dec_comb_ord,
            top_annotation=ch.upset_top_annotation(mk, ylim=(0, max_comb_d), height=ch.unit(1.6, "cm")),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_d)),
        )
        ht_dec1 = ch.add_heatmap(ht_dec1, u, direction="vertical")

    ht_dec2 = None
    for lbl, mk in zip(dec_labels[5:], m_dec[5:]):
        u = ch.UpSet(
            mk,
            row_title=lbl,
            set_order=dec_set_ord,
            comb_order=dec_comb_ord,
            top_annotation=ch.upset_top_annotation(mk, ylim=(0, max_comb_d), height=ch.unit(1.6, "cm")),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_d)),
        )
        ht_dec2 = ch.add_heatmap(ht_dec2, u, direction="vertical")

    g_dec1 = ch.grid_grabExpr(
        lambda: ch.draw(ht_dec1, newpage=False, padding=ch.unit([14, 4, 4, 4], "mm")),
        width_in=7.8,
        height_in=11.0,
    )
    g_dec2 = ch.grid_grabExpr(
        lambda: ch.draw(ht_dec2, newpage=False, padding=ch.unit([14, 4, 4, 4], "mm")),
        width_in=7.8,
        height_in=11.0,
    )

    def _draw_decades() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.49, just="left"))
        ch.grid_draw(g_dec1)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, just="left"))
        ch.grid_draw(g_dec2)
        ch.popViewport()

    save_Standalone(_draw_decades, "ch08_movies_by_decade.png", width=16.0, height=11.2)
    reset_state()

    # Relative fraction by decade
    ht_dec_rel1 = None
    for lbl, mk in zip(dec_labels[:5], m_dec[:5]):
        s = ch.comb_size(mk)
        rel = s / max(1.0, float(s.sum()))
        u = ch.UpSet(
            mk,
            row_title=lbl,
            set_order=dec_set_ord,
            comb_order=dec_comb_ord,
            top_annotation=ch.HeatmapAnnotation(
                **{
                    "Relative\nfraction": ch.anno_barplot(
                        rel, ylim=(0, 0.5), gp=ch.gpar(fill="black"), border=False, height=ch.unit(1.6, "cm")
                    )
                },
                annotation_name_side="left",
                annotation_name_rot=0,
            ),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_d)),
        )
        ht_dec_rel1 = ch.add_heatmap(ht_dec_rel1, u, direction="vertical")

    ht_dec_rel2 = None
    for lbl, mk in zip(dec_labels[5:], m_dec[5:]):
        s = ch.comb_size(mk)
        rel = s / max(1.0, float(s.sum()))
        u = ch.UpSet(
            mk,
            row_title=lbl,
            set_order=dec_set_ord,
            comb_order=dec_comb_ord,
            top_annotation=ch.HeatmapAnnotation(
                **{
                    "Relative\nfraction": ch.anno_barplot(
                        rel, ylim=(0, 0.5), gp=ch.gpar(fill="black"), border=False, height=ch.unit(1.6, "cm")
                    )
                },
                annotation_name_side="left",
                annotation_name_rot=0,
            ),
            right_annotation=ch.upset_right_annotation(mk, ylim=(0, max_set_d)),
        )
        ht_dec_rel2 = ch.add_heatmap(ht_dec_rel2, u, direction="vertical")

    g_dec_rel1 = ch.grid_grabExpr(
        lambda: ch.draw(ht_dec_rel1, newpage=False, padding=ch.unit([14, 4, 4, 4], "mm")),
        width_in=7.8,
        height_in=11.0,
    )
    g_dec_rel2 = ch.grid_grabExpr(
        lambda: ch.draw(ht_dec_rel2, newpage=False, padding=ch.unit([14, 4, 4, 4], "mm")),
        width_in=7.8,
        height_in=11.0,
    )

    def _draw_decades_rel() -> None:
        ch.pushViewport(ch.viewport(x=0.0, width=0.49, just="left"))
        ch.grid_draw(g_dec_rel1)
        ch.popViewport()
        ch.pushViewport(ch.viewport(x=0.51, width=0.49, just="left"))
        ch.grid_draw(g_dec_rel2)
        ch.popViewport()

    save_Standalone(_draw_decades_rel, "ch08_movies_by_decade_rel.png", width=16.0, height=11.2)
    reset_state()

    # Boxplots of years, ratings, watches per combination
    m_sub = ch.make_comb_mat(movies[genre])
    m_sub = m_sub[:, ch.comb_degree(m_sub) > 0]
    comb_elements = [ch.extract_comb(m_sub, nm) for nm in ch.comb_name(m_sub)]
    years_list = [movies.loc[ind, "ReleaseDate"].to_numpy() for ind in comb_elements]
    rating_list = [movies.loc[ind, "AvgRating"].to_numpy() for ind in comb_elements]
    watches_list = [movies.loc[ind, "Watches"].to_numpy() for ind in comb_elements]

    save_ht(
        ch.UpSet(
            ch.t(m_sub),
            right_annotation=ch.upset_right_annotation(
                ch.t(m_sub),
                annotation_name_side="top",
                axis_param={"side": "top"},
            ),
        )
        + ch.rowAnnotation(
            years=ch.anno_boxplot(years_list, width=ch.unit(2.2, "cm")),
            rating=ch.anno_boxplot(rating_list, width=ch.unit(2.2, "cm")),
            watches=ch.anno_boxplot(watches_list, width=ch.unit(2.2, "cm")),
            gap=ch.unit(3, "mm"),
            annotation_name_offset=ch.unit(8, "mm"),
        ),
        "ch08_movies_boxplots.png",
        width=9.2,
        height=7.8,
        padding=ch.unit([6, 14, 6, 4], "mm"),
    )
    reset_state()


def render_genomic_examples() -> None:
    """Render Section 8.7 genomic intervals UpSet examples."""
    peak_list = _build_genomic_peaks()
    m = ch.make_comb_mat(peak_list)
    m = m[:, ch.comb_size(m) > 500_000]

    save_ht(ch.UpSet(m), "ch08_genomic_basic.png", width=8.5, height=4.2)
    reset_state()

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.upset_top_annotation(
                m,
                axis_param={"at": [0, 5_000_000, 10_000_000, 15_000_000], "labels": ["0Mb", "5Mb", "10Mb", "15Mb"]},
                height=ch.unit(3.8, "cm"),
            ),
            right_annotation=ch.upset_right_annotation(
                m,
                axis_param={"at": [0, 10_000_000, 20_000_000, 30_000_000], "labels": ["0Mb", "10Mb", "20Mb", "30Mb"], "labels_rot": 0},
                width=ch.unit(3.8, "cm"),
            ),
        ),
        "ch08_genomic_axes.png",
        width=8.5,
        height=4.4,
    )
    reset_state()

    rng = np.random.default_rng(101)
    subgroup = {
        "ESC": "group1",
        "ES-deriv1": "group1",
        "ES-deriv2": "group1",
        "Brain": "group2",
        "Muscle": "group2",
        "Heart": "group2",
    }
    comb_sets = [ch.extract_comb(m, nm) for nm in ch.comb_name(m)]
    dist_to_tss = [
        np.abs(rng.normal(loc=float(rng.uniform(500, 2000)), scale=900, size=max(10, len(gr))))
        for gr in comb_sets
    ]
    mean_meth = [float(np.clip(np.mean(np.abs(rng.normal(0.12, 0.08, size=max(10, len(gr))))), 0.0, 1.0)) for gr in comb_sets]

    save_ht(
        ch.UpSet(
            m,
            top_annotation=ch.upset_top_annotation(
                m,
                axis_param={"at": [0, 5_000_000, 10_000_000, 15_000_000], "labels": ["0Mb", "5Mb", "10Mb", "15Mb"]},
                height=ch.unit(3.8, "cm"),
            ),
            right_annotation=ch.upset_right_annotation(
                m,
                axis_param={"at": [0, 10_000_000, 20_000_000, 30_000_000], "labels": ["0Mb", "10Mb", "20Mb", "30Mb"], "labels_rot": 0},
                width=ch.unit(3.8, "cm"),
            ),
            left_annotation=ch.rowAnnotation(
                group=[subgroup[s] for s in ch.set_name(m)],
                col={"group": {"group1": "#1B9E77", "group2": "#D95F02"}},
                show_annotation_name=False,
            ),
            bottom_annotation=ch.HeatmapAnnotation(
                dist_to_tss=ch.anno_boxplot(dist_to_tss, outline=False, height=ch.unit(1.8, "cm")),
                mean_meth=mean_meth,
                annotation_name_side="left",
                annotation_name_offset=ch.unit(9, "mm"),
            ),
        ),
        "ch08_genomic_annotated.png",
        width=10.5,
        height=5.6,
        padding=ch.unit([4, 14, 4, 4], "mm"),
    )
    reset_state()


def main() -> None:
    reset_state()
    render_upset_mode_diagram()
    render_basic_and_styling_plots()
    render_movies_examples()
    render_genomic_examples()
    print(f"Chapter 8 figures written to {IMAGES_DIR}")


if __name__ == "__main__":
    main()
