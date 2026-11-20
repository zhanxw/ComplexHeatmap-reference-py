# About {-}

This is the complete reference documentation for [**ComplexHeatmap-py**](https://github.com/zhanxw/ComplexHeatmap-py), a faithful Python port of the R/Bioconductor [**ComplexHeatmap**](https://github.com/jokergoo/ComplexHeatmap) package (upstream R version `2.25.3`). Examples throughout this book are generated under **complexheatmap** `0.2.0` backed by NumPy, SciPy, pandas, and Matplotlib.

## Installation

Install the latest release from PyPI or directly from GitHub:

```bash
pip install ComplexHeatmap-py
```

Or install the development version from GitHub:

```bash
pip install "git+https://github.com/zhanxw/ComplexHeatmap-py.git"
```

In Python scripts and notebooks, import the package using:

```python
import complexheatmap as ch
from complexheatmap import Heatmap, HeatmapAnnotation, colorRamp2, draw, gpar, unit
```

**ComplexHeatmap-py** mirrors the R **ComplexHeatmap** API closely so existing workflows translate nearly line-for-line. Like its R predecessor—which was originally inspired by [**pheatmap**](https://CRAN.R-project.org/package=pheatmap)—**ComplexHeatmap-py** provides familiar argument names and even includes a drop-in `ch.pheatmap()` translator (see [Chapter 10](integrate-with-other-packages.html)) that maps `pheatmap` arguments directly onto a composable `Heatmap` object.

## R-to-Python Migration Highlights

| Feature | R (`ComplexHeatmap`) | Python (`complexheatmap`) |
|---|---|---|
| Imports | `library(ComplexHeatmap)`, `library(circlize)`, `library(grid)` | `import complexheatmap as ch` |
| Indexing & slicing | 1-based (`1:n`, `ht[1:5, ]`, `pindex(m, 1:2, ...)`) | **0-based** (`0..n-1`, `ht[:5, :]`, `row_order`, `column_order`, `cell_fun(j, i, ...)`, `layer_fun(j, i, ...)`, `anno_mark(at=...)`, `pindex`) |
| Dots in identifiers | `grid.rect(...)`, `default.units = "native"`, `unit.c(...)` | `grid.rect(...)` or `grid_rect(...)`, `default_units="native"`, `unit_c(...)` |
| Horizontal & vertical concatenation | `ht1 + ht2 + rowAnnotation(...)`, `ht1 %v% ht2` | `ht1 + ht2 + rowAnnotation(...)`, `ht1 %v% ht2` (with `from complexheatmap import v`) or `ht1 / ht2` |
| Booleans, missing & named maps | `TRUE`, `FALSE`, `NULL`, `NA`, `c(A = "red")` | `True`, `False`, `None`, `np.nan`, `{"A": "red"}`, `heatmap_legend_param={...}` |
| Post-draw decorations | `decorate_heatmap_body("mat", { ... })` | `with decorate_heatmap_body("mat"): ...` (1-based `slice=1`, `row_slice=1`) or `lambda: ...` |
| Rich text & math | `gridtext::richtext_grob`, `gt_render()`, `expression(...)` | `gt_render(...)` (HTML/Markdown spans & `<br>`) and Matplotlib mathtext strings (`r"$\hat{\beta}$"`) |
| Saving figures | `png("out.png"); draw(ht); dev.off()` | `draw(ht, save="out.png", width=8, height=6, dpi=140)` or `ht.draw(save="out.png")` |

> **Indexing convention note**: All matrix row/column positions (`row_order()`, `column_order()`, `pindex()`, `cell_fun(j, i, ...)`, `layer_fun(j, i, ...)`, `anno_mark(at=...)`, and `ht[0:5, 0:5]`) are **0-based** in Python. Only visual slice identifiers in `decorate_*` calls (`slice=1`, `row_slice=1`, `column_slice=1`) and `kr`/`kc` in `layer_fun` remain 1-based slice numbers.

## Citation

If you use **ComplexHeatmap-py** or **ComplexHeatmap** in your publications, please cite:

Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics*, DOI: [10.1093/bioinformatics/btw313](https://doi.org/10.1093/bioinformatics/btw313).

![ComplexHeatmap-py Logo](images/ch00_logo.png)

## Python Environment

```python
import sys
import matplotlib
import numpy as np
import pandas as pd
import scipy
import complexheatmap as ch

print(f"Python: {sys.version.split()[0]}")
print(f"complexheatmap: {ch.__version__} (upstream R {ch.__upstream_r_version__})")
print(
    f"numpy: {np.__version__}, pandas: {pd.__version__}, "
    f"scipy: {scipy.__version__}, matplotlib: {matplotlib.__version__}"
)
```

```
Python: 3.14.10
complexheatmap: 0.2.0 (upstream R 2.25.3)
numpy: 2.5.4, pandas: 3.1.1, scipy: 1.18.1, matplotlib: 3.11.2
```
