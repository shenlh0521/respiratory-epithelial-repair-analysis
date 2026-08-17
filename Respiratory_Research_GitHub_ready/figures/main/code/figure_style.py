"""Shared visual system for revised BMC Figures 1–8."""
from __future__ import annotations

import matplotlib as mpl

COLORS = {
    "repair": "#D55E00",
    "repair_light": "#E69F73",
    "imbalance": "#A44200",
    "transition": "#CC79A7",
    "transition_dark": "#8E4B87",
    "cilia": "#0072B2",
    "ciliating": "#56B4E9",
    "maturation": "#005F73",
    "secretory": "#5F8F7B",
    "persistent": "#B54A2B",
    "regressive": "#4C78A8",
    "neutral": "#666666",
    "neutral_light": "#D9DDE0",
    "null": "#C9CDD0",
    "stress": "#7A5195",
    "oxidative": "#E69F00",
    "paper": "#FFFFFF",
    "ink": "#1F2529",
    "grid": "#E6E8EA",
}

FONT = {"panel": 13, "heading": 10, "axis": 9, "tick": 8, "legend": 8, "annotation": 8, "footnote": 7}
SIZE = {1: (7.09, 6.0), 2: (7.09, 5.5), 3: (7.09, 5.7), 4: (7.09, 5.2), 5: (7.09, 5.5), 6: (7.09, 5.5), 7: (7.09, 5.4), 8: (7.09, 3.8)}


def apply_style() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": FONT["tick"],
        "axes.titlesize": FONT["heading"], "axes.titleweight": "bold",
        "axes.labelsize": FONT["axis"], "axes.linewidth": 0.75,
        "axes.edgecolor": COLORS["ink"], "axes.labelcolor": COLORS["ink"],
        "xtick.labelsize": FONT["tick"], "ytick.labelsize": FONT["tick"],
        "xtick.color": COLORS["ink"], "ytick.color": COLORS["ink"],
        "legend.fontsize": FONT["legend"], "legend.frameon": False,
        "text.color": COLORS["ink"], "figure.facecolor": COLORS["paper"],
        "axes.facecolor": COLORS["paper"], "savefig.facecolor": COLORS["paper"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def panel(ax, label: str, heading: str | None = None) -> None:
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=FONT["panel"], clip_on=False)
    if heading:
        ax.set_title(heading, loc="left", pad=6)


def clean(ax, grid: bool = True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color=COLORS["grid"], linewidth=0.55, zorder=0)

