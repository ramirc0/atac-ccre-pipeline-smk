"""Shared matplotlib style for pipeline plots."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.style  # noqa: F401  (registers mpl.style)

FONT_STACK = [
    "Anthropic Sans Text",
    "Google Sans Flex",
    "Arimo",
    "Arial",
    "DejaVu Sans",
]

FORMATS = ("svg", "png")


def apply_style():
    """Set the project plotting style (call before importing pyplot)."""
    mpl.use("Agg")
    mpl.style.use("default")
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = FONT_STACK
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["figure.constrained_layout.use"] = True
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.spines.right"] = False


def save_figure(fig, path, **kwargs):
    """Write `fig` to `path` as both SVG and PNG.

    `path` MUST NOT carry a suffix; one is added per format. Returns the
    list of paths written.
    """
    stem = Path(path).with_suffix("")
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in FORMATS:
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(out, format=fmt, **kwargs)
        written.append(out)
    return written
