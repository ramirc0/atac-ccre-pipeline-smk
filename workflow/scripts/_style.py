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
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"


def despine(ax, categorical_x=False):
    """Offset left and bottom spines by 10 pt and trim them to the end ticks.

    Each continuous axis is fitted to its data without margins, then widened
    to the nearest ticks enclosing the data, so the trimmed spine never ends
    short of the data. A categorical x axis has no spine or tick marks; its
    labels carry the categories.
    """
    import seaborn as sns

    ax.margins(0)
    ax.autoscale_view()
    axes = [(ax.yaxis, ax.get_ylim, ax.set_ylim)]
    if not categorical_x:
        axes.append((ax.xaxis, ax.get_xlim, ax.set_xlim))
    for axis, get_lim, set_lim in axes:
        lo, hi = sorted(get_lim())
        ticks = axis.get_major_locator().tick_values(lo, hi)
        lo = max((t for t in ticks if t <= lo), default=lo)
        hi = min((t for t in ticks if t >= hi), default=hi)
        axis.set_ticks([t for t in ticks if lo <= t <= hi])
        set_lim(lo, hi)
    sns.despine(ax=ax, bottom=categorical_x, offset=10, trim=True)
    if categorical_x:
        ax.tick_params(axis="x", length=0)


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
