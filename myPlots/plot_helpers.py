"""
plot_helpers.py

Higher-level plotting building blocks layered on top of plotTOOL (pto).

The original analysis script repeated the same matplotlib scaffolding many
times: single-axis overlays, 4x2 / 4x1 E_miss panels, and 2x2 Q^2 panels,
each with the same loop of pto.plotTGE* calls. These helpers capture that
scaffolding once.

A `Series` describes one thing to draw (data histogram, a simulation line, or
a scaled+offset data marker set). The same Series objects can be reused across
every figure, which is what makes the driver script short.

This module deliberately avoids dataclasses / f-strings so it stays compatible
with the Python 2-style `from __future__ import division` header the original
used; it works unchanged on Python 3.
"""

from __future__ import division
import matplotlib.pyplot as plt
import plotTOOL as pto


class Series(object):
    """One drawable element of a figure.

    kind:
        'data'    -> pto.plotTGEStep   (unscaled reference histogram)
        'sim'     -> pto.plotTGELine   (scaled simulation curve)
        'data_ex' -> pto.plotTGEStepEx (scaled data with x-offset + marker)

    scale_ep / scale_epp let the same Series carry the (e,e'p) and (e,e'pp)
    normalisation factors; `draw` picks the right one based on the channel.
    `offset` is a base x-offset; figures can shrink it via offset_scale.
    """

    def __init__(self, label, root_file, color, kind,
                 scale_ep=1.0, scale_epp=1.0, offset=0.0, marker=None):
        self.label = label
        self.file = root_file
        self.color = color
        self.kind = kind
        self.scale_ep = scale_ep
        self.scale_epp = scale_epp
        self.offset = offset
        self.marker = marker


def _scale_for(series, channel, unit_scale):
    if unit_scale:
        return 1
    return series.scale_ep if channel == 'ep' else series.scale_epp


def draw(ax, series, hist, channel='epp', offset_scale=1.0, unit_scale=False):
    """Draw a single Series onto an axis for the given histogram/channel.

    unit_scale=True forces the normalisation to 1, used for ratio histograms
    (e.g. pMiss_epp) where no scaling should be applied.
    """
    if series.kind == 'data':
        pto.plotTGEStep(ax, series.file, hist, series.color)
    elif series.kind == 'sim':
        scale = _scale_for(series, channel, unit_scale)
        pto.plotTGELine(ax, series.file, hist, series.color, scale)
    elif series.kind == 'data_ex':
        scale = _scale_for(series, channel, unit_scale)
        pto.plotTGEStepEx(ax, series.file, hist, series.color,
                          scale, series.offset * offset_scale, series.marker)
    else:
        raise ValueError("Unknown Series kind: %r" % series.kind)


def annotate(ax, items):
    """Place text labels. Each item is a dict with keys:
    x, y, text, and optional color, fontsize, bbox.
    """
    for it in (items or []):
        ax.text(it['x'], it['y'], it['text'],
                color=it.get('color', 'black'),
                fontsize=it.get('fontsize', 15),
                bbox=it.get('bbox'))


def plot_overlay(pdf, hist, series, xlabel, ylabel, xlim, ylim,
                 annotations=None, channel='epp',
                 offset_scale=1.0, unit_scale=False,
                 xlabel_size=15, ylabel_size=15):
    """Single-axis figure with any number of overlaid Series."""
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(xlabel, fontsize=xlabel_size)
    ax.set_ylabel(ylabel, fontsize=ylabel_size)
    fig.tight_layout()
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    for s in series:
        draw(ax, s, hist, channel,
             offset_scale=offset_scale, unit_scale=unit_scale)
    annotate(ax, annotations)
    pdf.savefig(fig)
    return fig


def plot_emiss_4x2(pdf, series, var_label, hist_ep_prefix, hist_epp_prefix,
                   ylim_ep, ylim_epp, labely_ep, labely_epp, pmiss_labels,
                   xlim=(-0.15, 0.4), offset_scale=1.0):
    """4-row x 2-col E_miss panel: (e,e'p) on the left, (e,e'pp) on the right,
    one pmiss bin per row. Each Series in `series` is drawn in every panel.
    Histogram names are built as '<prefix>_<bin>' for bin = 1..4.
    """
    fig, ax = plt.subplots(4, 2, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=False)
    for r in range(4):
        ax[r, 0].set_ylabel('Counts', fontsize=15)
        ax[r, 1].set_ylabel('Counts', fontsize=15)
        ax[r, 0].set_ylim(0, ylim_ep[r])
        ax[r, 1].set_ylim(0, ylim_epp[r])
        ax[r, 0].text(-0.14, labely_ep[r], pmiss_labels[r], fontsize=9)
        ax[r, 1].text(-0.14, labely_epp[r], pmiss_labels[r], fontsize=9)
    ax[3, 0].set_xlabel(var_label, fontsize=15)
    ax[3, 1].set_xlabel(var_label, fontsize=15)
    ax[3, 0].set_xlim(*xlim)
    ax[3, 1].set_xlim(*xlim)
    ax[0, 0].set_title(r'$(e,e^{\prime}p)$', fontsize=15)
    ax[0, 1].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)
    for r in range(4):
        hep = '{0}_{1}'.format(hist_ep_prefix, r + 1)
        hepp = '{0}_{1}'.format(hist_epp_prefix, r + 1)
        for s in series:
            draw(ax[r, 0], s, hep, 'ep', offset_scale=offset_scale)
            draw(ax[r, 1], s, hepp, 'epp', offset_scale=offset_scale)
    pdf.savefig(fig)
    return fig


def plot_emiss_4x1(pdf, series, var_label, hist_prefix,
                   ylim, labely, pmiss_labels,
                   xlim=(-0.2, 0.4), offset_scale=1.0):
    """4-row x 1-col (e,e'pp) E_miss panel, one pmiss bin per row."""
    fig, ax = plt.subplots(4, 1, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=False)
    for r in range(4):
        ax[r].set_ylabel('Counts', fontsize=15)
        ax[r].set_ylim(0, ylim[r])
        ax[r].text(-0.14, labely[r], pmiss_labels[r], fontsize=9)
    ax[3].set_xlabel(var_label, fontsize=15)
    ax[3].set_xlim(*xlim)
    ax[0].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)
    for r in range(4):
        hist = '{0}_{1}'.format(hist_prefix, r + 1)
        for s in series:
            draw(ax[r], s, hist, 'epp', offset_scale=offset_scale)
    pdf.savefig(fig)
    return fig


def plot_q2_single(pdf, ref, sim, hist, ylabel,
                   xlim=(1.5, 5.0), ylim=(0.0, 0.34)):
    """Single quantity vs Q^2: data step + (optional) simulation line."""
    fig, ax = plt.subplots(1, 1)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel(r'$Q^{2} GeV$', fontsize=15)
    fig.tight_layout()
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    pto.plotTGEStepQ2(ax, ref.file, hist, ref.color)
    if sim is not None:
        pto.plotTGELineQ2(ax, sim.file, hist, sim.color)
    pdf.savefig(fig)
    return fig


def plot_q2_2x2(pdf, ref, sim, hist_prefix, ylabel, pmiss_labels,
                xlim=(1.5, 5.0), ylim=(0.0, 0.25),
                label_xy=(2.3, 0.2), ylabel_size=15,
                big_label=None, big_label_xy=(1.7, 0.12), draw_sim=True):
    """2x2 panel vs Q^2, one pmiss bin per panel in row-major order
    (bins 1..4 -> [0,0], [0,1], [1,0], [1,1]).

    NOTE: the original script placed the pmiss text labels for the [0,1] and
    [1,0] panels swapped in the g_sigma_* figures (they matched the data only
    in the h_Q2_epp figure). Labels here always follow the histogram drawn in
    each panel, which corrects that inconsistency.
    """
    fig, ax = plt.subplots(2, 2, gridspec_kw={'wspace': 0, 'hspace': 0},
                           sharex=True, sharey=True)
    ax[0, 0].set_ylabel(ylabel, fontsize=ylabel_size)
    ax[1, 0].set_ylabel(ylabel, fontsize=ylabel_size)
    ax[1, 0].set_xlabel(r'$Q^{2}$', fontsize=15)
    ax[1, 1].set_xlabel(r'$Q^{2}$', fontsize=15)
    ax[0, 0].set_xlim(*xlim)
    ax[0, 0].set_ylim(*ylim)
    if big_label is not None:
        ax[0, 0].text(big_label_xy[0], big_label_xy[1], big_label, fontsize=25)
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for (r, c), label, b in zip(positions, pmiss_labels, range(1, 5)):
        hist = '{0}_{1}'.format(hist_prefix, b)
        pto.plotTGEStepQ2(ax[r, c], ref.file, hist, ref.color)
        if draw_sim and sim is not None:
            pto.plotTGELineQ2(ax[r, c], sim.file, hist, sim.color)
        ax[r, c].text(label_xy[0], label_xy[1], label, fontsize=13)
    pdf.savefig(fig)
    return fig