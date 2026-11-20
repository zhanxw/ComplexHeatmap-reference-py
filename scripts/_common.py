"""Shared utilities for generating ComplexHeatmap-reference-py book figures."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Callable

# Ensure MPLCONFIGDIR points to a writable directory before importing matplotlib
_mpl_dir = Path("/tmp/mplconfig")
_mpl_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_dir))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import complexheatmap as ch
from complexheatmap.grid import GridDevice as _GridDevice

_orig_draw_raster = _GridDevice.draw_raster


def _patched_draw_raster(self: Any, *args: Any, **kwargs: Any) -> Any:
    n_before = len(self.ax.images)
    res = _orig_draw_raster(self, *args, **kwargs)
    for im in self.ax.images[n_before:]:
        if im.get_zorder() < 1:
            im.set_zorder(1)
    return res


_GridDevice.draw_raster = _patched_draw_raster  # type: ignore[method-assign]

REPO_ROOT = Path(__file__).resolve().parent.parent
CHAPTERS_DIR = REPO_ROOT / "chapters"
DOCS_DIR = REPO_ROOT / "docs"
BOOK_DIR = DOCS_DIR / "book"
IMAGES_DIR = BOOK_DIR / "images"
ASSETS_DIR = BOOK_DIR / "assets"

IMAGES_DIR.mkdir(parents=True, exist_ok=True)
ASSETS_DIR.mkdir(parents=True, exist_ok=True)


def reset_state() -> None:
    """Reset ComplexHeatmap global options and close open Matplotlib figures."""
    ch.ht_opt(RESET=True)
    plt.close("all")


def save_ht(
    obj: Any,
    filename: str,
    width: float = 7.5,
    height: float = 5.5,
    dpi: int = 140,
    **draw_kwargs: Any,
) -> Path:
    """Draw a Heatmap, HeatmapList, or high-level plot and save to docs/book/images/."""
    out = IMAGES_DIR / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        ch.draw(obj, save=str(out), width=width, height=height, dpi=dpi, **draw_kwargs)
    finally:
        plt.close("all")
    return out


def save_decorated(
    obj: Any,
    filename: str,
    decorate_fn: Callable[[Any], Any] | Callable[[], Any],
    width: float = 7.5,
    height: float = 5.5,
    dpi: int = 140,
    **draw_kwargs: Any,
) -> Path:
    """Draw obj, execute decorate_fn (taking 0 or 1 arg), then save to docs/book/images/."""
    import inspect

    out = IMAGES_DIR / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    draw_kwargs.pop("save", None)
    try:
        res = ch.draw(obj, width=width, height=height, dpi=dpi, **draw_kwargs)
        sig = inspect.signature(decorate_fn)
        if len(sig.parameters) == 0:
            decorate_fn()
        else:
            decorate_fn(res)
        fig = None
        if hasattr(res, "_device") and res._device is not None:
            fig = res._device.fig
        elif ch._ENV.last is not None and hasattr(ch._ENV.last, "_device"):
            fig = ch._ENV.last._device.fig
        else:
            fig = plt.gcf()
        fig.savefig(out, dpi=dpi, facecolor="white")
    finally:
        plt.close("all")
    return out


def save_Standalone(
    draw_fn: Callable[[], Any],
    filename: str,
    width: float = 6.0,
    height: float = 4.5,
    dpi: int = 140,
    tight: bool = False,
) -> Path:
    """Open a fresh grid device, run draw_fn(), and save to docs/book/images/."""
    out = IMAGES_DIR / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        dev = ch.grid_newpage(width_in=width, height_in=height, dpi=dpi)
        draw_fn()
        save_kwargs: dict[str, Any] = {"dpi": dpi, "facecolor": "white"}
        if tight:
            save_kwargs["bbox_inches"] = "tight"
        dev.fig.savefig(out, **save_kwargs)
    finally:
        plt.close("all")
    return out


save_standalone = save_Standalone
