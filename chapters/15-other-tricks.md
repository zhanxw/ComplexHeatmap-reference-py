# Other Tricks {#other-tricks}

## Set the same cell size for different heatmaps with different dimensions {#same-cell-size}

When exporting multiple heatmaps of different row or column counts to separate PNG or PDF files, you often want every cell to have the exact same physical dimensions (e.g., `5 mm` × `5 mm`). Because every non-body component in `ComplexHeatmap-py` (titles, dendrograms, annotations, axis labels, legends) has fixed physical size, the total figure height grows linearly with `nr`: `height_in = a * nr + b`.

`ht_size(ht, unit_name="inches", valueOnly=True)` computes the exact total `[width, height]` of a `Heatmap` or `HeatmapList` layout without even needing to render pixels to disk. Testing two row counts (`nr = 10` and `nr = 20`) recovers the exact slope (`5 mm / 25.4 mm/in = 0.1969 in/row`) and fixed overhead `b`:

```python
import numpy as np
import pandas as pd
import complexheatmap as ch

def make_fixed_cell_heatmap(nr: int, title: str | None = None) -> ch.Heatmap:
    m = pd.DataFrame(
        np.random.default_rng(123 + nr).normal(size=(nr, 10)),
        index=[f"G{i + 1}" for i in range(nr)],
        columns=list("abcdefghij"),
    )
    return ch.Heatmap(
        m,
        name="mat",
        height=ch.unit(5, "mm") * nr,
        width=ch.unit(5, "mm") * 10,
        column_title=title or f"{nr} rows, each row is 5 mm",
        top_annotation=ch.HeatmapAnnotation(bar=np.arange(1, 11)),
    )

xs = np.array([10.0, 20.0])
sizes = [ch.ht_size(make_fixed_cell_heatmap(int(nr), title="foo"), "inches", valueOnly=True) for nr in xs]
heights_in = [s[1] for s in sizes]
dev_w = max(s[0] for s in sizes) + 4.0 / 25.4  # add 4 mm horizontal page padding

slope, intercept = np.polyfit(xs, heights_in, 1)
print(f"height_in = {slope:.4f} * nr + {intercept:.4f}")
# Output: height_in = 0.1969 * nr + 1.0770
```

Adding `4.0 / 25.4` inches (`2 mm` top + `2 mm` bottom padding) to the figure height ensures both exported PNGs have identical `5 mm` cell heights:

```python
pad_in = 4.0 / 25.4
for nr, out_name in ((10, "ch15_cell_size_10rows.png"), (20, "ch15_cell_size_20rows.png")):
    ht = make_fixed_cell_heatmap(nr, title=f"{nr} rows, each row is 5mm")
    ch.draw(ht, save=out_name, width=float(dev_w), height=float(slope * nr + intercept + pad_in), padding=ch.unit([2, 2, 2, 2], "mm"))
```

![10-row heatmap exported at 5 mm per row](images/ch15_cell_size_10rows.png)

![20-row heatmap exported with the exact same 5 mm per row cell height](images/ch15_cell_size_20rows.png)

## Global options with ht_opt {#global-options-ht-opt}

`ht_opt` stores package-wide styling defaults for borders, fonts, paddings, and clustering behavior. In Python, you can:

- Inspect all options with `ht_opt()`.
- Set global options persistently via `ht_opt.heatmap_border = True` or `ht_opt(heatmap_border=True, ...)`.
- Scope temporary overrides safely using the `with ht_opt.local(...):` context manager (automatically restored on block exit).
- Reset all options back to defaults with `ht_opt(RESET=True)`.

```python
rng = np.random.default_rng(42)
m = pd.DataFrame(rng.normal(size=(10, 8)), index=[f"Gene{i + 1}" for i in range(10)], columns=[f"S{j + 1}" for j in range(8)])
ha = ch.HeatmapAnnotation(
    Group=["Ctrl"] * 4 + ["Treat"] * 4,
    Score=ch.anno_barplot(rng.uniform(0.2, 1.0, 8)),
    col={"Group": {"Ctrl": "#66C2A5", "Treat": "#FC8D62"}},
)

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
    ch.draw(ht, width=7.5, height=5.8)

ch.ht_opt(RESET=True)
```

![Heatmap styled inside an ht_opt.local() context manager](images/ch15_ht_opt_scoped.png)

## Capturing and composing heatmaps with grid_grabExpr {#grid-grabexpr-composition}

When you need to arrange completely independent heatmaps (with unrelated rows and columns) into multi-panel figures—for instance, pairing a 2D expression heatmap in panel A with a 3D pathway activity heatmap in panel B—`grid_grabExpr()` records any `draw(..., newpage=False)` call into a replayable `CapturedGrob` that can be drawn inside arbitrary `Viewport` coordinates via `grid_draw()`:

```python
rng = np.random.default_rng(7)
m1 = pd.DataFrame(rng.normal(size=(12, 8)), index=[f"R{i + 1}" for i in range(12)], columns=[f"A{j + 1}" for j in range(8)])
m2 = pd.DataFrame(rng.uniform(0, 10, size=(8, 6)), index=[f"P{i + 1}" for i in range(8)], columns=[f"B{j + 1}" for j in range(6)])

ht_a = ch.Heatmap(
    m1, name="Expr", col=ch.colorRamp2([-2, 0, 2], ["#2166AC", "white", "#B2182B"]),
    column_title="Panel A: Z-score Expression", row_km=2,
    row_names_gp=ch.gpar(fontsize=8), column_names_gp=ch.gpar(fontsize=8),
)
ht_b = ch.Heatmap3D(
    m2, name="Score", col=ch.colorRamp2([0, 5, 10], ["#FFFFCC", "#FD8D3C", "#BD0026"]),
    column_title="Panel B: 3D Pathway Activity",
    top_annotation=ch.HeatmapAnnotation(_spacer=ch.anno_empty(border=False, height=ch.unit(6, "mm")), show_annotation_name=False),
    bar_max_length=ch.unit(7, "mm"),
    row_names_gp=ch.gpar(fontsize=8), column_names_gp=ch.gpar(fontsize=8),
)

grob_a = ch.grid_grabExpr(lambda: ch.draw(ht_a, padding=ch.unit([3, 3, 3, 3], "mm"), newpage=False), width_in=5.5, height_in=5.0)
grob_b = ch.grid_grabExpr(lambda: ch.draw(ht_b, padding=ch.unit([3, 3, 5, 5], "mm"), newpage=False), width_in=5.5, height_in=5.0)

dev = ch.grid_newpage(width_in=11.0, height_in=5.0, dpi=140)
ch.grid_rect(gp=ch.gpar(fill="white", col=None))

ch.pushViewport(ch.Viewport(x=ch.unit(0.0, "npc"), y=ch.unit(0.0, "npc"), width=ch.unit(0.5, "npc"), height=ch.unit(1.0, "npc"), just=("left", "bottom")))
ch.grid_draw(grob_a)
ch.grid_text("A", x=ch.unit(4, "mm"), y=ch.unit(1, "npc") - ch.unit(4, "mm"), just=("left", "top"), gp=ch.gpar(fontsize=14, fontface="bold"))
ch.popViewport()

ch.pushViewport(ch.Viewport(x=ch.unit(0.5, "npc"), y=ch.unit(0.0, "npc"), width=ch.unit(0.5, "npc"), height=ch.unit(1.0, "npc"), just=("left", "bottom")))
ch.grid_draw(grob_b)
ch.grid_text("B", x=ch.unit(4, "mm"), y=ch.unit(1, "npc") - ch.unit(4, "mm"), just=("left", "top"), gp=ch.gpar(fontsize=14, fontface="bold"))
ch.popViewport()

dev.fig.savefig("ch15_grid_grab_compose.png", dpi=140, facecolor="white")
```

![Multi-panel figure combining a 2D heatmap and a 3D heatmap via grid_grabExpr](images/ch15_grid_grab_compose.png)
