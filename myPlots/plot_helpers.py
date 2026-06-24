"""
plot_helpers.py

Higher-level plotting building blocks layered on top of graph_io (ROOT I/O)
and plot_primitives (matplotlib drawing). Replaces the previous version of
this file, which called into plotTOOL.py's TGraph-name-based functions
directly; now every Series knows the (task_name, selection, axis_bin) triple
that names its graph (see graph_names.py), and plot_helpers reads, scales,
and draws plain (x, y, yerr) lists.

A `Series` describes one thing to draw (a data histogram, a simulation
line, or a scaled+offset data marker set). The same Series objects can be
reused across every figure, which is what keeps driver scripts short.
"""

from __future__ import division
import matplotlib.pyplot as plt

import graph_io
import graph_names
import plot_primitives as pp
import ratio


class Series(object):
    """One drawable element of a figure.

    kind:
        'data'    -> step + error bars   (unscaled reference histogram)
        'sim'     -> line + shaded band  (scaled simulation curve)
        'data_ex' -> scaled marker set with an x-offset (other nuclei)

    scale_ep / scale_epp let the same Series carry the (e,e'p) and (e,e'pp)
    normalization factors; draw() picks the right one based on `selection`.
    Graphs whose selection is 'ratio' (already a ratio quantity, e.g.
    pMiss_epp_over_pMiss_ep) are never scaled, regardless of these values.
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


def _scale_for(series, selection):
    if selection == 'ratio':
        return 1.0
    return series.scale_ep if selection == 'ep' else series.scale_epp


def get_xy_err(series, task_name, selection, axis_bin=(), offset_scale=1.0,
               integrated=False, pattern=None):
    """Reads, scales, and x-shifts one Series' graph. Returns (x, y, yerr).

    Set integrated=True and pass `pattern` to read from the integrated
    table's graphs instead of the differential ones.
    """
    if integrated:
        name = graph_names.integrated_graph_name(task_name, selection, pattern)
    else:
        name = graph_names.diff_graph_name(task_name, selection, list(axis_bin))

    x, y, yerr = graph_io.read_graph(series.file, name)
    scale = _scale_for(series, selection)
    shift = series.offset * offset_scale

    x = [xi + shift for xi in x]
    y = [yi * scale for yi in y]
    yerr = [ei * scale for ei in yerr]
    return x, y, yerr


def draw(ax, series, task_name, selection, axis_bin=(), offset_scale=1.0,
         integrated=False, pattern=None, q2_panel=False):
    """Draws one Series onto an axis for the given (task_name, selection,
    axis_bin) graph. q2_panel=True uses the wide-step/band padding used by
    the vs-Q^2 panels (matches the old plotTGEStepQ2/plotTGELineQ2 look)."""
    x, y, yerr = get_xy_err(series, task_name, selection, axis_bin,
                            offset_scale, integrated, pattern)

    if series.kind == 'data':
        if q2_panel:
            pp.step_with_error_q2(ax, x, y, yerr, series.color)
        else:
            pp.step_with_error(ax, x, y, yerr, series.color)
    elif series.kind == 'sim':
        if q2_panel:
            pp.line_with_band_q2(ax, x, y, yerr, series.color)
        else:
            pp.line_with_band(ax, x, y, yerr, series.color)
    elif series.kind == 'data_ex':
        pp.errorbar_marker(ax, x, y, yerr, series.color, series.marker)
    else:
        raise ValueError("Unknown Series kind: %r" % series.kind)

    return x, y, yerr


def annotate(ax, items):
    """Place text labels. Each item is a dict with keys:
    x, y, text, and optional color, fontsize, bbox."""
    for it in (items or []):
        ax.text(it['x'], it['y'], it['text'],
                color=it.get('color', 'black'),
                fontsize=it.get('fontsize', 15),
                bbox=it.get('bbox'))


def _draw_ratio_panel(ax_ratio, ref_xy, series_list, task_name, selection,
                      axis_bin, offset_scale, integrated, pattern):
    """Draws data/sim for each 'sim'-kind Series in series_list against the
    reference (data) series' (x, y, yerr), already read by the caller."""
    ref_x, ref_y, ref_err = ref_xy
    for s in series_list:
        if s.kind != 'sim':
            continue
        sim_x, sim_y, sim_err = get_xy_err(s, task_name, selection, axis_bin,
                                           offset_scale, integrated, pattern)
        if sim_x != ref_x:
            raise ValueError(
                "Cannot build a ratio panel: data and sim graphs for "
                "task '%s' have different binning (%d vs %d points). They "
                "should come from the same FillTask definition run over "
                "different input files." % (task_name, len(ref_x), len(sim_x)))
        r, e = ratio.ratio_series_with_error(ref_y, ref_err, sim_y, sim_err)
        ax_ratio.errorbar(ref_x, r, e, color=s.color, linestyle='', marker='o',
                          markersize=3)
    ax_ratio.axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
    ax_ratio.set_ylabel('Data/Sim', fontsize=10)


def plot_overlay(pdf, task_name, series, xlabel, ylabel, xlim, ylim,
                 selection='epp', axis_bin=(), annotations=None,
                 offset_scale=1.0, with_ratio=True, ratio_ylim=(0.5, 1.5),
                 xlabel_size=15, ylabel_size=15):
    """Single-axis figure with any number of overlaid Series, plus an
    optional ratio subpanel (data/sim) underneath.

    with_ratio is on by default for any figure that has both a 'data' and
    at least one 'sim' Series; set it False for ratio-quantity plots (where
    a second ratio-of-a-ratio panel wouldn't mean anything) or A-dependence
    plots that have no single simulation reference.
    """
    has_data = any(s.kind == 'data' for s in series)
    has_sim = any(s.kind == 'sim' for s in series)
    show_ratio = with_ratio and has_data and has_sim

    if show_ratio:
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    else:
        fig, ax = plt.subplots(1, 1)
        ax_ratio = None

    ax.set_ylabel(ylabel, fontsize=ylabel_size)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    ref_xy = None
    for s in series:
        x, y, yerr = draw(ax, s, task_name, selection, axis_bin, offset_scale)
        if s.kind == 'data':
            ref_xy = (x, y, yerr)
    annotate(ax, annotations)

    if show_ratio:
        _draw_ratio_panel(ax_ratio, ref_xy, series, task_name, selection,
                          axis_bin, offset_scale, integrated=False, pattern=None)
        ax_ratio.set_xlabel(xlabel, fontsize=xlabel_size)
        ax_ratio.set_xlim(*xlim)
        ax_ratio.set_ylim(*ratio_ylim)
    else:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)

    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def plot_emiss_4x2(pdf, series, var_label, task_ep, task_epp,
                   ylim_ep, ylim_epp, labely_ep, labely_epp, pmiss_labels,
                   xlim=(-0.15, 0.4), offset_scale=1.0):
    """4-row x 2-col E_miss panel: (e,e'p) on the left, (e,e'pp) on the
    right, one pMiss bin per row. task_ep/task_epp are the FillTask names
    (e.g. 'E1miss_ep_SRC_pmiss' / 'E1miss_epp_SRC_pmiss'); pMiss bin index
    `r` (0..3) is passed as axis_bin=[r] when reading each graph."""
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
        for s in series:
            draw(ax[r, 0], s, task_ep, 'ep', axis_bin=[r], offset_scale=offset_scale)
            draw(ax[r, 1], s, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
    pdf.savefig(fig)
    return fig


def plot_emiss_4x2_ratio_only(pdf, ref, sim, var_label, task_ep, task_epp,
                              pmiss_labels, xlim=(-0.15, 0.4),
                              ylim=(0.5, 1.5), labely=1.35,
                              offset_scale=1.0):
    """4-row x 2-col data/sim-only panel for E_miss plots.

    Errors are propagated bin-by-bin assuming independent data and
    simulation uncertainties.
    """
    fig, ax = plt.subplots(4, 2, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=True)
    for r in range(4):
        ax[r, 0].set_ylabel('Data/Sim', fontsize=13)
        ax[r, 1].set_ylabel('Data/Sim', fontsize=13)
        ax[r, 0].set_ylim(*ylim)
        ax[r, 1].set_ylim(*ylim)
        ax[r, 0].text(-0.14, labely, pmiss_labels[r], fontsize=9)
        ax[r, 1].text(-0.14, labely, pmiss_labels[r], fontsize=9)
        ax[r, 0].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
        ax[r, 1].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
    ax[3, 0].set_xlabel(var_label, fontsize=15)
    ax[3, 1].set_xlabel(var_label, fontsize=15)
    ax[3, 0].set_xlim(*xlim)
    ax[3, 1].set_xlim(*xlim)
    ax[0, 0].set_title(r'$(e,e^{\prime}p)$', fontsize=15)
    ax[0, 1].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)

    for r in range(4):
        ref_x_ep, ref_y_ep, ref_err_ep = get_xy_err(
            ref, task_ep, 'ep', axis_bin=[r], offset_scale=offset_scale)
        sim_x_ep, sim_y_ep, sim_err_ep = get_xy_err(
            sim, task_ep, 'ep', axis_bin=[r], offset_scale=offset_scale)
        if sim_x_ep != ref_x_ep:
            raise ValueError(
                "Cannot build data/sim panel for task '%s' bin %d: "
                "data and sim have different x binning (%d vs %d points)."
                % (task_ep, r, len(ref_x_ep), len(sim_x_ep)))
        ratio_ep, ratio_err_ep = ratio.ratio_series_with_error(
            ref_y_ep, ref_err_ep, sim_y_ep, sim_err_ep)
        ax[r, 0].errorbar(ref_x_ep, ratio_ep, ratio_err_ep,
                          color=sim.color, linestyle='', marker='o',
                          markersize=3)

        ref_x_epp, ref_y_epp, ref_err_epp = get_xy_err(
            ref, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
        sim_x_epp, sim_y_epp, sim_err_epp = get_xy_err(
            sim, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
        if sim_x_epp != ref_x_epp:
            raise ValueError(
                "Cannot build data/sim panel for task '%s' bin %d: "
                "data and sim have different x binning (%d vs %d points)."
                % (task_epp, r, len(ref_x_epp), len(sim_x_epp)))
        ratio_epp, ratio_err_epp = ratio.ratio_series_with_error(
            ref_y_epp, ref_err_epp, sim_y_epp, sim_err_epp)
        ax[r, 1].errorbar(ref_x_epp, ratio_epp, ratio_err_epp,
                          color=sim.color, linestyle='', marker='o',
                          markersize=3)

    pdf.savefig(fig)
    return fig


def plot_emiss_4x1(pdf, series, var_label, task_epp,
                   ylim, labely, pmiss_labels,
                   xlim=(-0.2, 0.4), offset_scale=1.0):
    """4-row x 1-col (e,e'pp) E_miss panel, one pMiss bin per row."""
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
        for s in series:
            draw(ax[r], s, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
    pdf.savefig(fig)
    return fig


def plot_emiss_4x1_ratio_only(pdf, ref, sim, var_label, task_epp,
                              pmiss_labels, xlim=(-0.2, 0.4),
                              ylim=(0.5, 1.5), labely=1.35,
                              offset_scale=1.0):
    """4-row x 1-col data/sim-only panel for E_miss(e,e'pp) plots.

    Errors are propagated bin-by-bin assuming independent data and
    simulation uncertainties.
    """
    fig, ax = plt.subplots(4, 1, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=True)
    for r in range(4):
        ax[r].set_ylabel('Data/Sim', fontsize=13)
        ax[r].set_ylim(*ylim)
        ax[r].text(-0.14, labely, pmiss_labels[r], fontsize=9)
        ax[r].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
    ax[3].set_xlabel(var_label, fontsize=15)
    ax[3].set_xlim(*xlim)
    ax[0].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)

    for r in range(4):
        ref_x, ref_y, ref_err = get_xy_err(
            ref, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
        sim_x, sim_y, sim_err = get_xy_err(
            sim, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
        if sim_x != ref_x:
            raise ValueError(
                "Cannot build data/sim panel for task '%s' bin %d: "
                "data and sim have different x binning (%d vs %d points)."
                % (task_epp, r, len(ref_x), len(sim_x)))
        ratio_y, ratio_err = ratio.ratio_series_with_error(
            ref_y, ref_err, sim_y, sim_err)
        ax[r].errorbar(ref_x, ratio_y, ratio_err,
                       color=sim.color, linestyle='', marker='o',
                       markersize=3)

    pdf.savefig(fig)
    return fig


def plot_q2_2x2_ratio(pdf, ref, sim, numerator_task, denominator_task, ylabel,
                     pmiss_labels, xlim=(1.5, 5.0), ylim=(0.0, 0.25),
                     label_xy=(2.3, 0.2), ylabel_size=15,
                     big_label=None, big_label_xy=(1.7, 0.12), draw_sim=True):
    """2x2 panel of a ratio quantity (e.g. epp/ep) vs Q^2, one pMiss bin per
    panel in row-major order (bins 0..3 -> [0,0], [0,1], [1,0], [1,1]).

    numerator_task/denominator_task must be a pair registered in
    Main_Figs_Binned.cpp's ratioSpecs (e.g. 'Q2_epp_SRC_pmiss' /
    'Q2_ep_SRC_pmiss'), sharing a single selector axis (pMiss) -- axis_bin
    picks the panel, and Q2 (the value axis of that task) is what ends up
    on the x-axis, read from the differential ratio graph BuildGraphs.cpp
    already computed (buildRatioDiffRows in BinnedHistStore.h).
    """
    ratio_task = numerator_task + '_over_' + denominator_task
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
    for bin_idx, ((r, c), label) in enumerate(zip(positions, pmiss_labels)):
        draw(ax[r, c], ref, ratio_task, 'ratio', axis_bin=[bin_idx], q2_panel=True)
        if draw_sim and sim is not None:
            draw(ax[r, c], sim, ratio_task, 'ratio', axis_bin=[bin_idx], q2_panel=True)
        ax[r, c].text(label_xy[0], label_xy[1], label, fontsize=13)
    pdf.savefig(fig)
    return fig


def plot_q2_2x2_data_over_sim_ratio_only(pdf, ref, sim,
                                         numerator_task, denominator_task,
                                         pmiss_labels,
                                         xlim=(1.5, 5.0), ylim=(0.5, 1.5),
                                         label_xy=(2.3, 1.35), ylabel_size=13,
                                         big_label=None,
                                         big_label_xy=(1.7, 1.12)):
    """2x2 panel of (data/sim) for a ratio quantity vs Q^2.

    Each panel uses one pMiss bin. Errors are propagated with the
    independent-uncertainty formula using the precomputed ratio graph errors
    from data and simulation.
    """
    ratio_task = numerator_task + '_over_' + denominator_task
    fig, ax = plt.subplots(2, 2, gridspec_kw={'wspace': 0, 'hspace': 0},
                           sharex=True, sharey=True)
    ax[0, 0].set_ylabel('Data/Sim', fontsize=ylabel_size)
    ax[1, 0].set_ylabel('Data/Sim', fontsize=ylabel_size)
    ax[1, 0].set_xlabel(r'$Q^{2}$', fontsize=15)
    ax[1, 1].set_xlabel(r'$Q^{2}$', fontsize=15)
    ax[0, 0].set_xlim(*xlim)
    ax[0, 0].set_ylim(*ylim)
    if big_label is not None:
        ax[0, 0].text(big_label_xy[0], big_label_xy[1], big_label, fontsize=25)

    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for bin_idx, ((r, c), label) in enumerate(zip(positions, pmiss_labels)):
        ref_x, ref_y, ref_err = get_xy_err(
            ref, ratio_task, 'ratio', axis_bin=[bin_idx])
        sim_x, sim_y, sim_err = get_xy_err(
            sim, ratio_task, 'ratio', axis_bin=[bin_idx])
        if sim_x != ref_x:
            raise ValueError(
                "Cannot build data/sim panel for task '%s' bin %d: "
                "data and sim have different x binning (%d vs %d points)."
                % (ratio_task, bin_idx, len(ref_x), len(sim_x)))
        ratio_y, ratio_err = ratio.ratio_series_with_error(
            ref_y, ref_err, sim_y, sim_err)
        ax[r, c].errorbar(ref_x, ratio_y, ratio_err,
                          color=sim.color, linestyle='', marker='o',
                          markersize=3)
        ax[r, c].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
        ax[r, c].text(label_xy[0], label_xy[1], label, fontsize=13)

    pdf.savefig(fig)
    return fig


# NOTE: the old g_sigma_pcmx / g_sigma_E1miss_*_pmiss / g_mean_E1miss_*_pmiss
# plots (C.M. width vs Q^2, E_miss mean/width vs Q^2) came from Gaussian
# fits across the 100 toy histograms, not from simple bin sums. That data
# was deliberately NOT put into diffTable/integratedTable -- only the
# nominal+toy histograms themselves were persisted (hists/nominal/...,
# hists/toy_NNN/...) for exactly this purpose. Reproducing those plots needs
# a separate fit-extraction step (read each toy's histogram, fit, collect
# mean/sigma across toys) that hasn't been built yet -- intentionally left
# out of plot_data_vs_sim.py / plot_A_dependence.py rather than faked here.
