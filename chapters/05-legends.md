# Legends {#legends}

Heatmaps and simple annotations automatically generate legends, which by default are placed on the right side of the plot. Complex annotations (such as `anno_points()`, `anno_lines()`, or `anno_barplot()`) do not generate legends automatically, though custom legends can be constructed and attached manually (see [Add customized legends](#add-customized-legends)). Internally, all legends are built with the `Legend()` constructor. In the sections below, we cover continuous legends, discrete legends, custom graphic legends, packing multiple legends, and customizing or positioning legends on heatmaps and annotations.

Every single legend or packed group of legends is an instance of the `Legends` class. It exposes a `.grob` attribute (`GTree`) containing the underlying grid primitives, along with `.width` and `.height` properties (also accessible via `ch.width(lgd)` and `ch.height(lgd)`). Because a `Legends` object encapsulates its layout and size, it can be drawn anywhere inside a grid viewport by specifying `x`, `y`, and `just`.

Legends created automatically for heatmaps and annotations are configured via `heatmap_legend_param` in `Heatmap()` and `annotation_legend_param` in `HeatmapAnnotation()`. **Most parameters accepted by `Legend()` can be passed directly in these two dictionaries**, as detailed in [Heatmap and annotation legends](#heatmap-and-annotation-legends).

## Continuous legends {#continuous-legends}

Because most heatmaps display continuous numeric matrices, we begin with continuous color-bar legends.

A continuous legend requires a color mapping function created with `colorRamp2()`. In heatmaps and simple annotations, this function is passed to `col` and forwarded to the legend builder automatically; when building a standalone legend, you pass `col_fun` directly to `Legend()`.

Notice that the breakpoints supplied to `colorRamp2()` (e.g., `[0, 0.5, 1]`) define the color interpolation scale, not necessarily the exact tick positions displayed in the legend. By default, `Legend()` chooses pretty break values aiming for roughly 5–6 labels.

First, here is the default vertical continuous legend:

```python
import complexheatmap as ch
from complexheatmap import Legend, colorRamp2, gpar, grid, gt_render, unit

col_fun = colorRamp2([0, 0.5, 1], ["blue", "white", "red"])
lgd = Legend(col_fun=col_fun, title="foo")
ch.draw(lgd)
```

![Default vertical continuous legend](images/ch05_01_cont_default.png)

`lgd` is a `Legends` object. Its physical dimensions can be queried with `ch.width()` and `ch.height()` (or `.width` and `.height`):

```python
ch.width(lgd)
# 10.1798mm
ch.height(lgd)
# 31.1887mm
```

Because a `Legends` object wraps a packed `GTree` of rectangles, tick segments, and text grobs, you can place it at arbitrary coordinates inside any active viewport using `draw(lgd, x=..., y=..., just=...)` (or `lgd.draw(...)`):

```python
ch.pushViewport(ch.viewport(width=0.9, height=0.9))
grid.rect()  # viewport border
ch.draw(lgd, x=unit(1, "cm"), y=unit(1, "cm"), just=["left", "bottom"])
ch.draw(lgd, x=unit(0.5, "npc"), y=unit(0.5, "npc"))
ch.draw(lgd, x=unit(1, "npc"), y=unit(1, "npc"), just=["right", "top"])
ch.popViewport()
```

![Drawing a standalone legend at multiple viewport positions](images/ch05_02_draw_positions.png)

If you only want to style legends produced by heatmaps or annotations, you do not need to call `Legend()` manually—pass the same parameters inside `heatmap_legend_param` or `annotation_legend_param`:

```python
import numpy as np

rng = np.random.default_rng(123)
mat = rng.standard_normal((10, 10))
ht = ch.Heatmap(
    mat,
    heatmap_legend_param={
        "title": "rnorm",
        "at": [-2, 0, 2],
        "labels": ["neg_two", "zero", "pos_two"],
    },
    top_annotation=ch.HeatmapAnnotation(
        foo=np.arange(1, 11),
        annotation_legend_param={"foo": {"title": "foo_top_anno"}},
    ),
)
ch.draw(ht)
```

![Customizing heatmap and annotation legends via legend_param](images/ch05_03_heatmap_anno_param.png)

In the remainder of this section, we focus on constructing `Legend()` objects; remember that calling `ch.draw(lgd)` (or `lgd.draw()`) renders any standalone `Legends` object.

Manually specify tick breakpoints via `at`:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.25, 0.5, 0.75, 1])
```

![Continuous legend with custom break positions](images/ch05_04_cont_at.png)

Set custom tick labels via `labels`:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.5, 1], labels=["low", "median", "high"])
```

![Continuous legend with custom break labels](images/ch05_05_cont_labels.png)

Control the height of a vertical continuous legend body (excluding the title) with `legend_height`:

```python
lgd = Legend(col_fun=col_fun, title="foo", legend_height=unit(6, "cm"))
```

![Vertical continuous legend with legend_height=6cm](images/ch05_06_cont_height.png)

Control the color-bar width with `grid_width`:

```python
lgd = Legend(col_fun=col_fun, title="foo", grid_width=unit(1, "cm"))
```

![Vertical continuous legend with grid_width=1cm](images/ch05_07_cont_width.png)

Customize label graphics with `labels_gp` (here `font=3` or `fontface="italic"` selects italics):

```python
lgd = Legend(col_fun=col_fun, title="foo", labels_gp=gpar(col="red", font=3))
```

![Continuous legend with custom labels_gp](images/ch05_08_cont_labels_gp.png)

Set the border and tick color with `border`:

```python
lgd = Legend(col_fun=col_fun, title="foo", border="red")
```

![Continuous legend with red border and ticks](images/ch05_09_cont_border.png)

`title_position` controls title placement. For vertical legends, valid values are `"topleft"`, `"topcenter"`, `"lefttop-rot"`, and `"leftcenter-rot"`:

```python
lgd = Legend(col_fun=col_fun, title="foooooooo", title_position="lefttop-rot", legend_height=unit(4, "cm"))
```

![Vertical continuous legend with title_position='lefttop-rot'](images/ch05_10_cont_title_lefttop_rot.png)

```python
lgd = Legend(col_fun=col_fun, title="foooooooo", title_position="leftcenter-rot", legend_height=unit(4, "cm"))
```

![Vertical continuous legend with title_position='leftcenter-rot'](images/ch05_11_cont_title_leftcenter_rot.png)

Mathematical symbols and equations can be rendered using Matplotlib's `$...$` mathtext syntax:

```python
lgd = Legend(
    col_fun=col_fun,
    title=r"$\hat{\beta} = (X^T X)^{-1} X^T y$",
    at=[0, 0.25, 0.5, 0.75, 1],
    labels=[r"$\alpha$", r"$\beta$", r"$\gamma$", r"$\delta$", r"$\epsilon$"],
)
```

![Continuous legend with mathematical formulas](images/ch05_12_cont_math.png)

Rich HTML/Markdown styling can be applied to the legend title using `gt_render()` together with `title_gp=gpar(box_fill=...)`:

```python
col_fun_gt = colorRamp2([-3, 0, 3], ["blue", "white", "red"])
lgd = Legend(
    col_fun=col_fun_gt,
    title=gt_render("<span style='color:orange'>**Legend title**</span>"),
    title_gp=gpar(box_fill="grey"),
    at=[-3, 0, 3],
    labels=["negative three", "zero", "positive three"],
)
```

![Continuous legend with gt_render title](images/ch05_13_cont_gridtext.png)

> **Note:** In the current `ComplexHeatmap-py` implementation, `gt_render()` renders rich HTML/Markdown formatting on `Legend()` titles, while `labels` inside `Legend()` are converted to plain strings. For full `gt_render()` label styling on heatmap titles, row/column names, and `anno_text()`, see [Integrate with other packages](integrate-with-other-packages.html).

Horizontal continuous legends share nearly the same parameters, except `legend_width` sets the bar width and `title_position` must be one of `"topleft"`, `"topcenter"`, `"lefttop"`, or `"leftcenter"`:

```python
lgd = Legend(col_fun=col_fun, title="foo", direction="horizontal")
```

![Default horizontal continuous legend](images/ch05_14_horiz_default.png)

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.25, 0.5, 0.75, 1], direction="horizontal")
```

![Horizontal continuous legend with custom breaks](images/ch05_15_horiz_at.png)

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.5, 1], labels=["low", "median", "high"], direction="horizontal")
```

![Horizontal continuous legend with custom labels](images/ch05_16_horiz_labels.png)

```python
lgd = Legend(col_fun=col_fun, title="foo", legend_width=unit(6, "cm"), direction="horizontal")
```

![Horizontal continuous legend with legend_width=6cm](images/ch05_17_horiz_width.png)

```python
lgd = Legend(col_fun=col_fun, title="foo", labels_gp=gpar(col="red", font=3), direction="horizontal")
```

![Horizontal continuous legend with custom labels_gp](images/ch05_18_horiz_labels_gp.png)

Rotate horizontal legend labels using `labels_rot`:

```python
lgd = Legend(col_fun=col_fun, title="foo", labels_rot=45, direction="horizontal")
```

![Horizontal continuous legend with labels_rot=45](images/ch05_19_horiz_rot45.png)

Horizontal title positions `"topcenter"` and `"lefttop"`:

```python
lgd = Legend(col_fun=col_fun, title="foooooooo", direction="horizontal", title_position="topcenter")
```

![Horizontal continuous legend with title_position='topcenter'](images/ch05_20_horiz_topcenter.png)

```python
lgd = Legend(col_fun=col_fun, title="foooooooo", direction="horizontal", title_position="lefttop")
```

![Horizontal continuous legend with title_position='lefttop'](images/ch05_21_horiz_lefttop.png)

`at` can contain breakpoints with unequal spacing. When labels would otherwise overlap, `Legend()` keeps the tick marks at their true positions along the color bar, shifts the text labels cleanly via `smartAlign`, and draws connector segments from the ticks to the labels:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1])
```

![Vertical continuous legend with unequal break intervals](images/ch05_22_unequal_vert.png)

If the legend is tall enough that labels do not collide, no displacement is performed:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.3, 1], legend_height=unit(4, "cm"))
```

![Tall vertical legend where unequal labels do not collide](images/ch05_23_unequal_no_shift.png)

The same automatic label dodging works for horizontal legends and rotated labels:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1], direction="horizontal")
```

![Horizontal continuous legend with unequal intervals](images/ch05_24_unequal_horiz.png)

```python
lgd = Legend(
    col_fun=col_fun,
    title="foo",
    at=[0, 0.1, 0.15, 0.5, 0.9, 0.95, 1],
    direction="horizontal",
    title_position="lefttop",
    labels_rot=90,
)
```

![Horizontal legend with rotated labels and unequal intervals](images/ch05_25_unequal_horiz_rot90.png)

When `title_position="lefttop"`, the space beneath the title is also available for shifted labels:

```python
lgd = Legend(
    col_fun=col_fun,
    title="foo",
    at=[0, 0.1, 0.5, 0.75, 1],
    labels=["minimal", "q10", "median", "q75", "maximal"],
    direction="horizontal",
    title_position="lefttop",
)
```

![Horizontal legend extending first label under lefttop title](images/ch05_26_unequal_horiz_lefttop.png)

Providing `at` in decreasing order reverses the continuous legend so the smallest value appears at the top:

```python
lgd = Legend(col_fun=col_fun, title="foo", at=[1, 0.8, 0.6, 0.4, 0.2, 0])
```

![Reversed vertical continuous legend](images/ch05_27_reversed_at.png)

Normally, visual distances along the color bar are proportional to numeric differences between breaks. Consider a proportion color scale with breakpoints `[0, 0.05, 0.1, 0.5, 1]`, where we want extra visual detail in `[0, 0.1]`. By default, `Legend()` chooses uniform breaks across `[0, 1]`:

```python
col_fun_prop = colorRamp2([0, 0.05, 0.1, 0.5, 1], ["green", "white", "red", "black", "blue"])
lgd = Legend(col_fun=col_fun_prop, title="Prop")
```

![Proportion color scale with default uniform breaks](images/ch05_28_break_dist_default.png)

Even if we specify `at=[0, 0.05, 0.1, 0.5, 1]`, the `[0.5, 1]` interval still occupies five times as much bar height as `[0, 0.1]`:

```python
lgd = Legend(col_fun=col_fun_prop, title="Prop", at=[0, 0.05, 0.1, 0.5, 1])
```

![Proportion color scale with manual at breaks](images/ch05_29_break_dist_manual_at.png)

The `break_dist` parameter controls the **visual distance between neighboring breaks in the legend**. Pass a scalar (`break_dist=1`) to give every break interval equal visual span, or pass a sequence of length `len(at) - 1`:

```python
lgd = Legend(col_fun=col_fun_prop, title="Prop", break_dist=1)
```

![Equal visual break spacing via break_dist=1](images/ch05_30_break_dist_equal.png)

Here the upper two intervals are given three times the visual distance of the lower two intervals (`break_dist=[1, 1, 3, 3]`):

```python
lgd = Legend(col_fun=col_fun_prop, title="Prop", break_dist=[1, 1, 3, 3])
```

![Non-uniform visual spacing via break_dist=[1, 1, 3, 3]](images/ch05_31_break_dist_1133.png)

Increasing `legend_height` eliminates the need for label displacement:

```python
lgd = Legend(col_fun=col_fun_prop, title="Prop", break_dist=[1, 1, 3, 3], legend_height=unit(4, "cm"))
```

![Taller legend with break_dist=[1, 1, 3, 3]](images/ch05_32_break_dist_1133_h4.png)

Combined with a tiny epsilon step (`0.1 + 1e-6`), `break_dist` also lets you visualize piecewise color schemes cleanly across `[0, 0.1]` and `(0.1, 1]`:

```python
col_fun2 = colorRamp2([0, 0.1, 0.1 + 1e-6, 1], ["white", "red", "yellow", "blue"])
lgd = Legend(
    col_fun=col_fun2,
    title="Prop",
    at=[0, 0.05, 0.1, 0.5, 1],
    break_dist=[1, 1, 3, 3],
    legend_height=unit(4, "cm"),
)
```

![Piecewise dual color scheme with break_dist](images/ch05_33_break_dist_twocolor.png)

## Discrete legends {#discrete-legends}

Discrete legends represent categorical levels—or discretized slices of a continuous scale. Supply `at` and/or `labels`, along with `legend_gp=gpar(fill=...)`:

```python
lgd = Legend(at=list(range(1, 7)), title="foo", legend_gp=gpar(fill=list(range(1, 7))))
```

![Basic discrete legend with numeric breaks](images/ch05_34_disc_at.png)

```python
MONTH_NAMES = ["January", "February", "March", "April", "May", "June",
               "July", "August", "September", "October", "November", "December"]
lgd = Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7))))
```

![Discrete legend with month name labels](images/ch05_35_disc_months.png)

Discretizing a continuous color mapping:

```python
at = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
lgd = Legend(at=at, title="foo", legend_gp=gpar(fill=col_fun(at).tolist()))
```

![Discrete legend sampling a continuous colorRamp2](images/ch05_36_disc_from_cont.png)

Title placement options on discrete legends (`"lefttop"` and `"leftcenter-rot"`):

```python
lgd = Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7))), title_position="lefttop")
```

![Discrete legend with title_position='lefttop'](images/ch05_37_disc_lefttop.png)

```python
lgd = Legend(labels=MONTH_NAMES[:6], title="foo", legend_gp=gpar(fill=list(range(1, 7))), title_position="leftcenter-rot")
```

![Discrete legend with title_position='leftcenter-rot'](images/ch05_38_disc_leftcenter_rot.png)

Control swatch sizes with `grid_height` and `grid_width`:

```python
lgd = Legend(
    at=list(range(1, 7)),
    legend_gp=gpar(fill=list(range(1, 7))),
    title="foo",
    grid_height=unit(1, "cm"),
    grid_width=unit(5, "mm"),
)
```

![Discrete legend with custom swatch dimensions](images/ch05_39_disc_grid_size.png)

Style label and title typography via `labels_gp` and `title_gp`:

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    legend_gp=gpar(fill=list(range(1, 7))),
    title="foo",
    labels_gp=gpar(col="red", fontsize=14),
)
```

![Discrete legend with customized labels_gp](images/ch05_40_disc_labels_gp.png)

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    legend_gp=gpar(fill=list(range(1, 7))),
    title="foo",
    title_gp=gpar(col="red", fontsize=14),
)
```

![Discrete legend with customized title_gp](images/ch05_41_disc_title_gp.png)

Rich HTML/Markdown titles with `gt_render()`:

```python
lgd = Legend(
    title=gt_render("<span style='color:orange'>**Legend title**</span>"),
    title_gp=gpar(box_fill="grey"),
    at=[-3, 0, 3],
    labels=["negative three", "zero", "positive three"],
    legend_gp=gpar(fill=[1, 2, 3]),
)
```

![Discrete legend with gt_render title](images/ch05_42_disc_gridtext.png)

Control swatch borders with `border`:

```python
lgd = Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foo", border="red")
```

![Discrete legend with red swatch borders](images/ch05_43_disc_border.png)

Discrete legends can be arranged across multiple columns (`ncol`) or rows (`nrow`). Notice you cannot set `ncol` and `nrow` simultaneously:

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3)
```

![Multi-column discrete legend with ncol=3](images/ch05_44_disc_ncol3.png)

Title centering adapts automatically to multi-column widths:

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, title_position="topcenter")
```

![Multi-column legend with centered title](images/ch05_45_disc_ncol3_topcenter.png)

Fill entries row-by-row with `by_row=True`:

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, by_row=True)
```

![Multi-column legend filled by row](images/ch05_46_disc_ncol3_byrow.png)

Adjust horizontal spacing between columns with `gap` (or `column_gap`), and vertical spacing between rows with `row_gap`:

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, gap=unit(1, "cm"))
```

![Multi-column legend with gap=1cm](images/ch05_47_disc_ncol3_gap.png)

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", ncol=3, row_gap=unit(5, "mm"))
```

![Multi-column legend with row_gap=5mm](images/ch05_48_disc_ncol3_rowgap.png)

Specify row count with `nrow`:

```python
lgd = Legend(labels=MONTH_NAMES[:10], legend_gp=gpar(fill=list(range(1, 11))), title="foo", nrow=3)
```

![Multi-row discrete legend with nrow=3](images/ch05_49_disc_nrow3.png)

Single-row horizontal discrete legends paired with `"lefttop-rot"` or `"leftcenter"` titles:

```python
lgd = Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foooooo", nrow=1, title_position="lefttop-rot")
```

![Single-row legend with lefttop-rot title](images/ch05_50_disc_nrow1_lefttoprot.png)

```python
lgd = Legend(labels=MONTH_NAMES[:6], legend_gp=gpar(fill=list(range(1, 7))), title="foooooo", nrow=1, title_position="leftcenter")
```

![Single-row legend with leftcenter title](images/ch05_51_disc_nrow1_leftcenter.png)

`Legend()` also supports simple geometric primitives via `type`:

- `type="points"` (or `"p"`) with numeric or character `pch`
- `type="lines"` (or `"l"`)
- `type="boxplot"` (or `"box"`)

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    title="foo",
    type="points",
    pch=list(range(1, 7)),
    legend_gp=gpar(col=list(range(1, 7))),
    background="#FF8080",
)
```

![Points legend with numeric pch](images/ch05_52_disc_points_num.png)

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    title="foo",
    type="points",
    pch=list("abcdef"),
    legend_gp=gpar(col=list(range(1, 7))),
    background="white",
)
```

![Points legend with character pch](images/ch05_53_disc_points_char.png)

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    title="foo",
    type="lines",
    legend_gp=gpar(col=list(range(1, 7)), lty=list(range(1, 7))),
    grid_width=unit(1, "cm"),
)
```

![Lines legend](images/ch05_54_disc_lines.png)

```python
lgd = Legend(
    labels=MONTH_NAMES[:6],
    title="foo",
    type="boxplot",
    legend_gp=gpar(fill=list(range(1, 7))),
)
```

![Boxplot legend](images/ch05_55_disc_boxplot.png)

Special integer `pch` codes `26`, `27`, and `28` draw diagonal `/` (`26`), `\` (`27`), and crossed `X` (`28`) segments:

```python
lgd = Legend(labels=[f"pch = {p}" for p in (26, 27, 28)], type="points", pch=[26, 27, 28])
```

![Points legend with special diagonal pch codes 26-28](images/ch05_56_disc_pch26_28.png)

Multi-line labels automatically enlarge the corresponding legend grid height (and equalize row heights across multi-column rows):

```python
lgd = Legend(labels=["aaaaa\naaaaa", "bbbbb\nbbbbb", "c", "d"], legend_gp=gpar(fill=[1, 2, 3, 4]))
```

![Discrete legend with multi-line labels](images/ch05_57_disc_multiline.png)

```python
lgd = Legend(labels=["aaaaa\naaaaa", "c", "d", "bbbbb\nbbbbb"], legend_gp=gpar(fill=[1, 2, 3, 4]), nrow=2)
```

![Multi-column legend with multi-line labels](images/ch05_58_disc_multiline_nrow2.png)

Finally, `graphics` accepts a list (or a dict keyed by `at`) of custom drawing callbacks with signature `(x, y, w, h)`, giving you full control over each swatch cell:

```python
lgd = Legend(
    labels=list("abcd"),
    graphics=[
        lambda x, y, w, h: grid.rect(x, y, w * 0.33, h, gp=gpar(fill="red")),
        lambda x, y, w, h: grid.rect(x, y, w, h * 0.33, gp=gpar(fill="blue")),
        lambda x, y, w, h: grid.text("A", x, y, gp=gpar(col="darkgreen")),
        lambda x, y, w, h: grid.points(x, y, gp=gpar(col="orange"), pch=16),
    ],
)
```

![Discrete legend with custom graphics callbacks](images/ch05_59_disc_custom_graphics.png)

## A list of legends {#a-list-of-legends}

Multiple `Legends` objects can be arranged together using `packLegend()`, either as positional arguments or via `list=[...]`. Internally, `ComplexHeatmap-py` uses `packLegend()` to lay out all heatmap and annotation legends.

```python
from complexheatmap import packLegend

lgd1 = Legend(at=list(range(1, 7)), legend_gp=gpar(fill=list(range(1, 7))), title="legend1")
lgd2 = Legend(col_fun=col_fun, title="legend2", at=[0, 0.25, 0.5, 0.75, 1])
lgd3 = Legend(labels=MONTH_NAMES[:3], legend_gp=gpar(fill=[7, 8, 9]), title="legend3")

pd_vert = packLegend(lgd1, lgd2, lgd3)
# Equivalent to:
# pd_vert = packLegend(list=[lgd1, lgd2, lgd3])
```

![Vertically packed legends](images/ch05_60_pack_vertical.png)

Check packed dimensions just like a single legend:

```python
ch.width(pd_vert)
# 20.181mm
ch.height(pd_vert)
# 83.8748mm
```

Pack legends horizontally with `direction="horizontal"`:

```python
pd_horiz = packLegend(lgd1, lgd2, lgd3, direction="horizontal")
```

![Horizontally packed legends](images/ch05_61_pack_horizontal.png)

When vertical legends exceed `max_height`, `packLegend()` wraps them into multiple columns separated by `column_gap`:

```python
pd_mh = packLegend(
    lgd1, lgd3, lgd2, lgd3, lgd2, lgd1,
    max_height=unit(10, "cm"),
    column_gap=unit(1, "cm"),
)
```

![Vertically packed legends wrapping into multiple columns](images/ch05_62_pack_max_height.png)

Similarly, horizontal packing wraps into multiple rows when exceeding `max_width`, spaced by `column_gap` and `row_gap`:

```python
lgd1_row = Legend(at=list(range(1, 7)), legend_gp=gpar(fill=list(range(1, 7))), title="legend1", nrow=1)
lgd2_row = Legend(col_fun=col_fun, title="legend2", at=[0, 0.25, 0.5, 0.75, 1], direction="horizontal")

pd_mw = packLegend(
    lgd1_row, lgd2_row, lgd3, lgd1_row, lgd2_row, lgd3,
    max_width=unit(10, "cm"),
    direction="horizontal",
    column_gap=unit(5, "mm"),
    row_gap=unit(1, "cm"),
)
```

![Horizontally packed legends wrapping into multiple rows](images/ch05_63_pack_max_width.png)

Because `packLegend()` returns a `Legends` object, it can also be drawn at any viewport position:

```python
pd_pos = packLegend(lgd1_row, lgd2_row, lgd3, direction="horizontal")
ch.pushViewport(ch.viewport(width=0.8, height=0.8))
grid.rect()
ch.draw(pd_pos, x=unit(1, "cm"), y=unit(1, "cm"), just=["left", "bottom"])
ch.draw(pd_pos, x=unit(1, "npc"), y=unit(1, "npc"), just=["right", "top"])
ch.popViewport()
```

![Drawing packed legends inside a viewport](images/ch05_64_pack_draw_positions.png)

## Heatmap and annotation legends {#heatmap-and-annotation-legends}

Pass legend options for a heatmap matrix via `heatmap_legend_param` in `Heatmap()`:

```python
rng = np.random.default_rng(123)
m = rng.standard_normal((10, 10))
ht = ch.Heatmap(
    m,
    name="mat",
    heatmap_legend_param={
        "at": [-2, 0, 2],
        "labels": ["low", "zero", "high"],
        "title": "Some values",
        "legend_height": unit(4, "cm"),
        "title_position": "lefttop-rot",
    },
)
ch.draw(ht)
```

![Heatmap with customized heatmap_legend_param](images/ch05_65_hm_legend_param.png)

For annotations, `annotation_legend_param` maps each annotation name to its parameter dictionary:

```python
ha = ch.HeatmapAnnotation(
    foo=rng.uniform(size=10),
    bar=rng.choice(["f", "m"], 10, replace=True),
    annotation_legend_param={
        "foo": {"title": "Fooooooh", "at": [0, 0.5, 1], "labels": ["zero", "median", "one"]},
        "bar": {"title": "Baaaaaaar", "at": ["f", "m"], "labels": ["Female", "Male"]},
    },
)
ch.draw(ch.Heatmap(m, name="mat", top_annotation=ha))
```

![Heatmap with customized annotation_legend_param](images/ch05_66_anno_legend_param.png)

When heatmaps are concatenated horizontally (`+`), heatmap legends and row-annotation legends are grouped together while column-annotation legends form a separate group:

```python
ha1 = ch.HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True))
ha2 = ch.HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
ht_list = (
    ch.Heatmap(m, name="mat1", top_annotation=ha1)
    + ch.rowAnnotation(sth=rng.uniform(size=10))
    + ch.Heatmap(m, name="mat2", top_annotation=ha2)
)
ch.draw(ht_list, padding=unit([10, 2, 2, 2], "mm"))
```

![Horizontal heatmap concatenation legend grouping](images/ch05_67_horiz_grouping.png)

Conversely, in vertical concatenation (`% v %`), heatmap and column-annotation legends are grouped together while row-annotation legends form the second group:

```python
from complexheatmap import v

ha1 = ch.HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True), annotation_name_side="left")
ha2 = ch.HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
ht_list_v = (
    ch.Heatmap(m, name="mat1", top_annotation=ha1)
    % v
    % ch.Heatmap(m, name="mat2", top_annotation=ha2, right_annotation=ch.rowAnnotation(sth=np.arange(1, 11)))
)
ch.draw(ht_list_v, padding=unit([10, 2, 2, 2], "mm"))
```

![Vertical heatmap concatenation legend grouping](images/ch05_68_vert_grouping.png)

Control legend visibility with `show_legend` in `HeatmapAnnotation()` (accepts a boolean, sequence, or dictionary keyed by annotation name) and `show_heatmap_legend=False` in `Heatmap()`:

```python
ha = ch.HeatmapAnnotation(
    foo=rng.uniform(size=10),
    bar=rng.choice(["f", "m"], 10, replace=True),
    show_legend=[True, False],  # or show_legend={"bar": False}
    annotation_name_side="left",
)
ch.draw(
    ch.Heatmap(m, name="mat1", top_annotation=ha)
    + ch.Heatmap(m, name="mat2", show_heatmap_legend=False)
)
```

![Hiding selected annotation and heatmap legends](images/ch05_69_hide_legends.png)

Set `merge_legends=True` in `draw()` to combine heatmap and annotation legends into a single packed group:

```python
ha1 = ch.HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True))
ha2 = ch.rowAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(list("abc"), 10, replace=True))
ha3 = ch.rowAnnotation(foo3=rng.uniform(size=10), bar3=rng.choice(MONTH_NAMES[:3], 10, replace=True))
ht_list = ch.Heatmap(m, name="mat1", top_annotation=ha1) + ch.Heatmap(m, name="mat2", left_annotation=ha2) + ha3
ch.draw(ht_list, merge_legends=True)
```

![Merging all heatmap and annotation legends with merge_legends=True](images/ch05_70_merge_legends.png)

Set `legend_grouping="original"` to separate all annotation legends (both row and column) from pure heatmap legends:

```python
ch.draw(ht_list, legend_grouping="original")
```

![Separating pure heatmap and annotation legends with legend_grouping='original'](images/ch05_71_legend_grouping_original.png)

A continuous color mapping can be rendered as discrete legend swatches by setting `color_bar="discrete"`:

```python
ch.draw(
    ch.Heatmap(
        m,
        name="mat",
        heatmap_legend_param={"color_bar": "discrete"},
        top_annotation=ch.HeatmapAnnotation(
            foo=np.arange(1, 11),
            annotation_legend_param={"foo": {"color_bar": "discrete"}},
        ),
    )
)
```

![Rendering continuous color scales with discrete legend swatches](images/ch05_72_color_bar_discrete.png)

For categorical annotation vectors, default colors preserve first-seen level order; you can enforce an explicit level order in the legend via a `col` dictionary or `at` inside `annotation_legend_param` (and `heatmap_legend_param`):

```python
chr_vec = rng.choice(list("abc"), 10, replace=True)
ch.draw(
    ch.Heatmap(
        m,
        name="mat",
        top_annotation=ch.HeatmapAnnotation(
            chr=chr_vec,
            fa_ordered=chr_vec,
            col={"fa_ordered": {"c": "#1B9E77", "a": "#D95F02", "b": "#7570B3"}},
            annotation_legend_param={"fa_ordered": {"at": ["b", "c", "a"]}},
        ),
    )
)
```

![Controlling discrete legend level ordering with col and at](images/ch05_73_level_order.png)

## Add customized legends {#add-customized-legends}

Custom `Legend()` objects can be appended to the heatmap legend group via `heatmap_legend_list` or to the annotation legend group via `annotation_legend_list` in `draw()` (see also [Other High-level Plots](other-high-level-plots.html)).

Because `anno_points()`, `anno_lines()`, and `anno_summary()` do not automatically create legends, drawing the combined plot alone omits legends for those tracks:

```python
ha1 = ch.HeatmapAnnotation(
    pt=ch.anno_points(np.arange(1, 11), gp=gpar(col=np.repeat([2, 3], 5)), height=unit(2, "cm")),
    show_annotation_name=False,
)
ha2 = ch.HeatmapAnnotation(
    ln=ch.anno_lines(np.column_stack([np.arange(1, 11), np.arange(10, 0, -1)]), gp=gpar(col=[4, 5], lty=[1, 2]), height=unit(2, "cm")),
    show_annotation_name=False,
)
ht_list = (
    ch.Heatmap(m, name="mat1", top_annotation=ha1)
    + ch.Heatmap(m, name="mat2", top_annotation=ha2)
    + ch.Heatmap(
        m[:, [0]],
        name="mat3",
        top_annotation=ch.HeatmapAnnotation(summary=ch.anno_summary(gp=gpar(fill=[2, 3]))),
        width=unit(1, "cm"),
    )
)
ch.draw(ht_list, ht_gap=unit(7, "mm"), row_km=2)
```

![Complex annotations without manually added legends](images/ch05_74_custom_annos_no_legend.png)

Construct matching `Legend()` objects and pass them to `annotation_legend_list`:

```python
lgd_list = [
    Legend(labels=["red", "green"], title="pt", type="points", pch=16, legend_gp=gpar(col=[2, 3])),
    Legend(labels=["darkblue", "lightblue"], title="ln", type="lines", legend_gp=gpar(col=[4, 5], lty=[1, 2])),
    Legend(labels=["group1", "group2"], title="km", type="boxplot", legend_gp=gpar(fill=[2, 3])),
]
ch.draw(ht_list, ht_gap=unit(7, "mm"), row_km=2, annotation_legend_list=lgd_list)
```

![Attaching custom Legends via annotation_legend_list](images/ch05_75_custom_annos_with_legend.png)

## The side of legends {#the-side-of-legends}

By default, both heatmap and annotation legends sit on the right side of the plot. Use `heatmap_legend_side` and `annotation_legend_side` in `draw()` to place them on `"left"`, `"right"`, `"bottom"`, or `"top"`:

```python
ha1 = ch.HeatmapAnnotation(foo1=rng.uniform(size=10), bar1=rng.choice(["f", "m"], 10, replace=True), annotation_name_side="left")
ha2 = ch.HeatmapAnnotation(foo2=rng.uniform(size=10), bar2=rng.choice(["f", "m"], 10, replace=True))
ht_list = ch.Heatmap(m, name="mat1", top_annotation=ha1) + ch.rowAnnotation(sth=rng.uniform(size=10), annotation_name_side="top") + ch.Heatmap(m, name="mat2", top_annotation=ha2)
ch.draw(
    ht_list,
    heatmap_legend_side="left",
    annotation_legend_side="bottom",
    ht_gap=unit(4, "mm"),
    padding=unit([2, 2, 2, 12], "mm"),
)
```

![Placing heatmap legends on the left and annotation legends at the bottom](images/ch05_76_legend_sides.png)

When placing legends at the top or bottom, you can also make each individual continuous legend horizontal (`"legend_direction": "horizontal"`) and each discrete legend single-row (`"nrow": 1`) inside `heatmap_legend_param` and `annotation_legend_param`:

> **Python note:** Standalone `Legend()` uses `direction="horizontal"`, whereas `heatmap_legend_param` and `annotation_legend_param` forward arguments to `ColorMapping.color_mapping_legend()`, where the parameter name is `legend_direction="horizontal"`.

```python
ha1 = ch.HeatmapAnnotation(
    foo1=rng.uniform(size=10),
    bar1=rng.choice(["f", "m"], 10, replace=True),
    annotation_name_side="left",
    annotation_legend_param={"foo1": {"legend_direction": "horizontal"}, "bar1": {"nrow": 1}},
)
ha2 = ch.HeatmapAnnotation(
    foo2=rng.uniform(size=10),
    bar2=rng.choice(["f", "m"], 10, replace=True),
    annotation_legend_param={"foo2": {"legend_direction": "horizontal"}, "bar2": {"nrow": 1}},
)
ht_list = (
    ch.Heatmap(m, name="mat1", top_annotation=ha1, heatmap_legend_param={"legend_direction": "horizontal"})
    + ch.rowAnnotation(sth=rng.uniform(size=10), annotation_name_side="top", annotation_legend_param={"sth": {"legend_direction": "horizontal"}})
    + ch.Heatmap(m, name="mat2", top_annotation=ha2, heatmap_legend_param={"legend_direction": "horizontal"})
)
ch.draw(
    ht_list,
    merge_legends=True,
    heatmap_legend_side="bottom",
    annotation_legend_side="bottom",
    ht_gap=unit(4, "mm"),
    padding=unit([2, 10, 2, 12], "mm"),
)
```

![All legends horizontal and merged at the bottom](images/ch05_77_all_horizontal_bottom.png)
