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
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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


def _drop_extreme_ytick(ax, which):
    """Drops the min or max y tick to avoid overlap across stacked panels."""
    ticks = list(ax.get_yticks())
    if len(ticks) <= 2:
        return
    if which == 'min':
        ticks = ticks[1:]
    elif which == 'max':
        ticks = ticks[:-1]
    ax.set_yticks(ticks)


def _auto_panel_label(ax, text, fontsize=9, manual_pos=None):
    """Places a panel label in a low-density corner in axes coordinates.

    manual_pos can be a tuple (x, y) or (x, y, ha, va) in axes coordinates.
    """
    if manual_pos is not None:
        if len(manual_pos) == 2:
            x, y = manual_pos
            ha, va = 'left', 'top'
        else:
            x, y, ha, va = manual_pos
        ax.text(x, y, text, transform=ax.transAxes, fontsize=fontsize,
                ha=ha, va=va)
        return

    # Candidate corners: (x, y, ha, va)
    candidates = [
        (0.04, 0.88, 'left', 'top'),
        (0.96, 0.88, 'right', 'top'),
        (0.04, 0.12, 'left', 'bottom'),
        (0.96, 0.12, 'right', 'bottom'),
    ]

    # Build point cloud from all Line2D artists already drawn on this axis.
    points_axes = []
    for line in ax.lines:
        xdata = line.get_xdata()
        ydata = line.get_ydata()
        for x, y in zip(xdata, ydata):
            try:
                xa, ya = ax.transAxes.inverted().transform(ax.transData.transform((x, y)))
                points_axes.append((xa, ya))
            except Exception:
                continue

    # Approximate text box size in axes coordinates; pick least-overlapping corner.
    box_w = 0.34
    box_h = 0.14
    best = candidates[0]
    best_score = None
    for x, y, ha, va in candidates:
        left = x if ha == 'left' else (x - box_w)
        right = left + box_w
        bottom = (y - box_h) if va == 'top' else y
        top = bottom + box_h
        score = 0
        for px, py in points_axes:
            if left <= px <= right and bottom <= py <= top:
                score += 1
        if best_score is None or score < best_score:
            best = (x, y, ha, va)
            best_score = score

    ax.text(best[0], best[1], text, transform=ax.transAxes, fontsize=fontsize,
            ha=best[2], va=best[3])


def _scale_for(series, selection):
    if selection == 'ratio':
        return 1.0
    return series.scale_ep if selection == 'ep' else series.scale_epp


def get_xy_err(series, task_name, selection, axis_bin=(), offset_scale=1.0,
               integrated=False, pattern=None, unit_scale=False):
    """Reads, scales, and x-shifts one Series' graph. Returns (x, y, yerr).

    Set integrated=True and pass `pattern` to read from the integrated
    table's graphs instead of the differential ones.
    """
    if integrated:
        name = graph_names.integrated_graph_name(task_name, selection, pattern)
    else:
        name = graph_names.diff_graph_name(task_name, selection, list(axis_bin))

    x, y, yerr = graph_io.read_graph(series.file, name)
    scale = 1.0 if unit_scale else _scale_for(series, selection)
    shift = series.offset * offset_scale

    x = [xi + shift for xi in x]
    y = [yi * scale for yi in y]
    yerr = [ei * scale for ei in yerr]
    return x, y, yerr


def draw(ax, series, task_name, selection, axis_bin=(), offset_scale=1.0,
         integrated=False, pattern=None, q2_panel=False, unit_scale=False):
    """Draws one Series onto an axis for the given (task_name, selection,
    axis_bin) graph. q2_panel=True uses the wide-step/band padding used by
    the vs-Q^2 panels (matches the old plotTGEStepQ2/plotTGELineQ2 look)."""
    x, y, yerr = get_xy_err(series, task_name, selection, axis_bin,
                            offset_scale, integrated, pattern, unit_scale)

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


def _save(pdf, fig, save_as=None):
    """Always saves fig as the next page of the combined book `pdf`. If
    `save_as` is given (a path), also writes fig out as its own single-figure
    PDF -- the per-call opt-in for "save this one plot individually" without
    restructuring the driver script."""
    pdf.savefig(fig)
    if save_as:
        fig.savefig(save_as)
    return fig


def _draw_ratio_panel(ax_ratio, ref_xy, series_list, task_name, selection,
                      axis_bin, offset_scale, integrated, pattern,
                      unit_scale=False):
    """Draws data/sim for each 'sim'-kind Series in series_list against the
    reference (data) series' (x, y, yerr), already read by the caller."""
    ref_x, ref_y, ref_err = ref_xy
    for s in series_list:
        if s.kind != 'sim':
            continue
        sim_x, sim_y, sim_err = get_xy_err(s, task_name, selection, axis_bin,
                                           offset_scale, integrated, pattern,
                                           unit_scale)
        if len(sim_x) != len(ref_x):
            raise ValueError(
                "Cannot build a ratio panel: data and sim graphs for "
                "task '%s' have different point counts (%d vs %d). They "
                "should come from the same FillTask definition run over "
                "different input files." % (task_name, len(ref_x), len(sim_x)))

        # For some derived quantities (e.g. fit-derived vs-Q^2 graphs), x can
        # differ slightly between files even when points are index-aligned.
        # In that case, compute ratios point-by-point and plot at the data x.
        r, e = ratio.ratio_series_with_error(ref_y, ref_err, sim_y, sim_err)
        ax_ratio.errorbar(ref_x, r, e, color=s.color, linestyle='', marker='o',
                          markersize=3)
    ax_ratio.axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
    ax_ratio.set_ylabel('Data/Sim', fontsize=10)


def plot_overlay(pdf, task_name, series, xlabel, ylabel, xlim, ylim,
                 selection='epp', axis_bin=(), annotations=None,
                 offset_scale=1.0, with_ratio=True, ratio_ylim=(0.5, 1.5),
                 xlabel_size=15, ylabel_size=15, unit_scale=False,
                 corner_label=None, save_as=None):
    """Single-axis figure with any number of overlaid Series, plus an
    optional ratio subpanel (data/sim) underneath.

    with_ratio is on by default for any figure that has both a 'data' and
    at least one 'sim' Series; set it False for ratio-quantity plots (where
    a second ratio-of-a-ratio panel wouldn't mean anything) or A-dependence
    plots that have no single simulation reference.

    corner_label, if given, is drawn inside the top axis (e.g. "(e,e'p)")
    in the top-left corner, in axes coordinates.
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
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    ref_xy = None
    for s in series:
        x, y, yerr = draw(ax, s, task_name, selection, axis_bin, offset_scale,
                          unit_scale=unit_scale)
        if s.kind == 'data':
            ref_xy = (x, y, yerr)
    annotate(ax, annotations)
    if corner_label is not None:
        ax.text(0.04, 0.92, corner_label, transform=ax.transAxes,
                fontsize=15, ha='left', va='top')

    if show_ratio:
        _draw_ratio_panel(ax_ratio, ref_xy, series, task_name, selection,
                          axis_bin, offset_scale, integrated=False,
                          pattern=None, unit_scale=unit_scale)
        ax_ratio.set_xlabel(xlabel, fontsize=xlabel_size)
        ax_ratio.set_xlim(*xlim)
        ax_ratio.set_ylim(*ratio_ylim)
        ax_ratio.yaxis.set_major_locator(MaxNLocator(nbins=4))
        _drop_extreme_ytick(ax, 'min')
        _drop_extreme_ytick(ax_ratio, 'max')
    else:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)

    fig.tight_layout()
    return _save(pdf, fig, save_as)


def plot_emiss_4x2(pdf, series, var_label, task_ep, task_epp,
                   ylim_ep, ylim_epp, labely_ep, labely_epp, pmiss_labels,
                   xlim=(-0.15, 0.4), offset_scale=1.0,
                   label_positions=None, save_as=None):
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
        ax[r, 0].yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax[r, 1].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual_ep = None if label_positions is None else label_positions.get(('ep', r))
        manual_epp = None if label_positions is None else label_positions.get(('epp', r))
        _auto_panel_label(ax[r, 0], pmiss_labels[r], fontsize=9, manual_pos=manual_ep)
        _auto_panel_label(ax[r, 1], pmiss_labels[r], fontsize=9, manual_pos=manual_epp)
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
    return _save(pdf, fig, save_as)


def plot_emiss_4x2_ratio_only(pdf, ref, sims, var_label, task_ep, task_epp,
                              pmiss_labels, xlim=(-0.15, 0.4),
                              ylim=(0.5, 1.5), labely=1.35,
                              offset_scale=1.0, label_positions=None,
                              save_as=None):
    """4-row x 2-col data/sim-only panel for E_miss plots.

    `sims` is a list of 'sim'-kind Series (0 or more); each gets its own
    overlaid data/sim trace per panel, colored by that Series' color.
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
        ax[r, 0].yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax[r, 1].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual_ep = None if label_positions is None else label_positions.get(('ep', r))
        manual_epp = None if label_positions is None else label_positions.get(('epp', r))
        _auto_panel_label(ax[r, 0], pmiss_labels[r], fontsize=9, manual_pos=manual_ep)
        _auto_panel_label(ax[r, 1], pmiss_labels[r], fontsize=9, manual_pos=manual_epp)
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
        for sim in sims:
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
        for sim in sims:
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

    return _save(pdf, fig, save_as)


def plot_emiss_4x1(pdf, series, var_label, task_epp,
                   ylim, labely, pmiss_labels,
                   xlim=(-0.2, 0.4), offset_scale=1.0,
                   label_positions=None, save_as=None):
    """4-row x 1-col (e,e'pp) E_miss panel, one pMiss bin per row."""
    fig, ax = plt.subplots(4, 1, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=False)
    for r in range(4):
        ax[r].set_ylabel('Counts', fontsize=15)
        ax[r].set_ylim(0, ylim[r])
        ax[r].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(r)
        _auto_panel_label(ax[r], pmiss_labels[r], fontsize=9, manual_pos=manual)
    ax[3].set_xlabel(var_label, fontsize=15)
    ax[3].set_xlim(*xlim)
    ax[0].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)
    for r in range(4):
        for s in series:
            draw(ax[r], s, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
    return _save(pdf, fig, save_as)


def plot_emiss_4x1_ratio_only(pdf, ref, sims, var_label, task_epp,
                              pmiss_labels, xlim=(-0.2, 0.4),
                              ylim=(0.5, 1.5), labely=1.35,
                              offset_scale=1.0, label_positions=None,
                              save_as=None):
    """4-row x 1-col data/sim-only panel for E_miss(e,e'pp) plots.

    `sims` is a list of 'sim'-kind Series (0 or more), each overlaid as its
    own data/sim trace per panel. Errors are propagated bin-by-bin assuming
    independent data and simulation uncertainties.
    """
    fig, ax = plt.subplots(4, 1, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=True)
    for r in range(4):
        ax[r].set_ylabel('Data/Sim', fontsize=13)
        ax[r].set_ylim(*ylim)
        ax[r].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(r)
        _auto_panel_label(ax[r], pmiss_labels[r], fontsize=9, manual_pos=manual)
        ax[r].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
    ax[3].set_xlabel(var_label, fontsize=15)
    ax[3].set_xlim(*xlim)
    ax[0].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)

    for r in range(4):
        ref_x, ref_y, ref_err = get_xy_err(
            ref, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)
        for sim in sims:
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

    return _save(pdf, fig, save_as)


def _row_box_label(fig, ax_left, ax_right, text, fontsize=11):
    """Places one boxed label centered above a row, spanning ax_left and
    ax_right (or just ax_left if ax_right is None), sitting just inside the
    top of the panel -- used by the *_note Emiss panels below in place of
    the corner label plot_emiss_4x2/4x1 use, to match an analysis-note
    figure where the pMiss-bin range is a single centered box per row."""
    pos_l = ax_left.get_position()
    pos_r = ax_right.get_position() if ax_right is not None else pos_l
    x_center = 0.5 * (pos_l.x0 + pos_r.x1)
    y = pos_l.y1 - 0.08 * (pos_l.y1 - pos_l.y0)
    fig.text(x_center, y, text, ha='center', va='center', fontsize=fontsize,
             bbox=dict(boxstyle='square', facecolor='white', edgecolor='black'))


def plot_emiss_4x2_note(pdf, series, var_label, task_ep, task_epp,
                        pmiss_centers_ep, pmiss_centers_epp, box_labels, mN,
                        xlim=(-0.15, 0.4), offset_scale=1.0, save_as=None):
    """Like plot_emiss_4x2, but styled for an analysis-note figure instead
    of the existing makePlots.py/plot_data_vs_sim.py/plot_A_dependence.py
    look (kept as a separate function so those callers' figures don't
    change): a single boxed pMiss-bin label centered above each row instead
    of a corner label, a thin vertical line at the quasi-free
    single-nucleon-knockout threshold E_miss = sqrt(pMiss^2 + mN^2) - mN
    (evaluated at each row's mean pMiss -- pmiss_centers_ep/epp are
    per-row averages measured from data, not bin-edge midpoints, since the
    e'p and e'pp populations' pMiss distributions within a nominal bin can
    differ slightly), per-row y-limits computed from the data+sim peak
    instead of caller-tuned ylim arrays, and fewer/dropped y-ticks to
    reduce clutter between stacked rows.
    """
    fig, ax = plt.subplots(4, 2, gridspec_kw={'hspace': 0, 'wspace': 0},
                           sharex=True, sharey=False)
    ax[3, 0].set_xlabel(var_label, fontsize=15)
    ax[3, 1].set_xlabel(var_label, fontsize=15)
    ax[3, 0].set_xlim(*xlim)
    ax[3, 1].set_xlim(*xlim)
    ax[0, 0].set_title(r'$(e,e^{\prime}p)$', fontsize=15)
    ax[0, 1].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)

    for r in range(4):
        ax[r, 0].set_ylabel('Counts', fontsize=13)
        ax[r, 1].set_ylabel('Counts', fontsize=13)
        # wspace=0 above butts the two columns together at a shared
        # boundary line; move the right column's ticks/label to its own
        # right-hand side so they don't collide with that boundary.
        ax[r, 1].yaxis.tick_right()
        ax[r, 1].yaxis.set_label_position('right')
        ax[r, 0].yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax[r, 1].yaxis.set_major_locator(MaxNLocator(nbins=3))

        peak_ep, peak_epp = 0.0, 0.0
        for s in series:
            _, y_ep, e_ep = get_xy_err(s, task_ep, 'ep', axis_bin=[r],
                                       offset_scale=offset_scale)
            _, y_epp, e_epp = get_xy_err(s, task_epp, 'epp', axis_bin=[r],
                                         offset_scale=offset_scale)
            peak_ep = max([peak_ep] + [v + e for v, e in zip(y_ep, e_ep)])
            peak_epp = max([peak_epp] + [v + e for v, e in zip(y_epp, e_epp)])
            draw(ax[r, 0], s, task_ep, 'ep', axis_bin=[r], offset_scale=offset_scale)
            draw(ax[r, 1], s, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)

        ax[r, 0].set_ylim(0, peak_ep * 1.2 if peak_ep > 0.0 else 1.0)
        ax[r, 1].set_ylim(0, peak_epp * 1.2 if peak_epp > 0.0 else 1.0)
        if r != 3:
            _drop_extreme_ytick(ax[r, 0], 'max')
            _drop_extreme_ytick(ax[r, 1], 'max')

        thresh_ep = math.sqrt(pmiss_centers_ep[r] ** 2 + mN ** 2) - mN
        thresh_epp = math.sqrt(pmiss_centers_epp[r] ** 2 + mN ** 2) - mN
        ax[r, 0].axvline(thresh_ep, color='black', linewidth=0.8)
        ax[r, 1].axvline(thresh_epp, color='black', linewidth=0.8)

        _row_box_label(fig, ax[r, 0], ax[r, 1], box_labels[r])

    return _save(pdf, fig, save_as)


def plot_emiss_4x1_note(pdf, series, var_label, task_epp,
                        pmiss_centers, box_labels, mN,
                        xlim=(-0.2, 0.4), offset_scale=1.0, save_as=None):
    """Like plot_emiss_4x1, but styled to match plot_emiss_4x2_note above
    (see there for what changed and why) -- used for E2miss, which has no
    e'p counterpart since it requires a detected recoil."""
    fig, ax = plt.subplots(4, 1, gridspec_kw={'hspace': 0, 'wspace': 0.3},
                           sharex=True, sharey=False)
    ax[3].set_xlabel(var_label, fontsize=15)
    ax[3].set_xlim(*xlim)
    ax[0].set_title(r'$(e,e^{\prime}pp)$', fontsize=15)

    for r in range(4):
        ax[r].set_ylabel('Counts', fontsize=13)
        ax[r].yaxis.set_major_locator(MaxNLocator(nbins=3))

        peak = 0.0
        for s in series:
            _, y, e = get_xy_err(s, task_epp, 'epp', axis_bin=[r],
                                 offset_scale=offset_scale)
            peak = max([peak] + [v + ei for v, ei in zip(y, e)])
            draw(ax[r], s, task_epp, 'epp', axis_bin=[r], offset_scale=offset_scale)

        ax[r].set_ylim(0, peak * 1.2 if peak > 0.0 else 1.0)
        if r != 3:
            _drop_extreme_ytick(ax[r], 'max')

        thresh = math.sqrt(pmiss_centers[r] ** 2 + mN ** 2) - mN
        ax[r].axvline(thresh, color='black', linewidth=0.8)

        _row_box_label(fig, ax[r], None, box_labels[r])

    return _save(pdf, fig, save_as)


def plot_q2_2x2_ratio(pdf, ref, sims, numerator_task, denominator_task, ylabel,
                     pmiss_labels, xlim=(1.5, 5.0), ylim=(0.0, 0.25),
                     label_xy=(2.3, 0.2), ylabel_size=15,
                     big_label=None, big_label_xy=(1.7, 0.12),
                     label_positions=None, save_as=None):
    """2x2 panel of a ratio quantity (e.g. epp/ep) vs Q^2, one pMiss bin per
    panel in row-major order (bins 0..3 -> [0,0], [0,1], [1,0], [1,1]).

    `sims` is a list of 'sim'-kind Series (0 or more) overlaid on top of
    `ref`. numerator_task/denominator_task must be a pair registered in
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
        for sim in sims:
            draw(ax[r, c], sim, ratio_task, 'ratio', axis_bin=[bin_idx], q2_panel=True)
        ax[r, c].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(bin_idx)
        _auto_panel_label(ax[r, c], label, fontsize=13, manual_pos=manual)
    return _save(pdf, fig, save_as)


def plot_q2_2x2_data_over_sim_ratio_only(pdf, ref, sims,
                                         numerator_task, denominator_task,
                                         pmiss_labels,
                                         xlim=(1.5, 5.0), ylim=(0.5, 1.5),
                                         label_xy=(2.3, 1.35), ylabel_size=13,
                                         big_label=None,
                                         big_label_xy=(1.7, 1.12),
                                         label_positions=None, save_as=None):
    """2x2 panel of (data/sim) for a ratio quantity vs Q^2.

    `sims` is a list of 'sim'-kind Series (0 or more); each gets its own
    overlaid trace per panel. Each panel uses one pMiss bin. Errors are
    propagated with the independent-uncertainty formula using the
    precomputed ratio graph errors from data and simulation.
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
        for sim in sims:
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
        ax[r, c].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(bin_idx)
        _auto_panel_label(ax[r, c], label, fontsize=13, manual_pos=manual)

    return _save(pdf, fig, save_as)


def plot_q2_single(pdf, ref, sims, task_name, ylabel,
                   selection='epp', axis_bin=(),
                   xlim=(1.5, 5.0), ylim=None,
                   with_ratio=False, ratio_ylim=(0.5, 1.5),
                   xlabel=r'$Q^{2}$', ylabel_size=15, unit_scale=False,
                   suptitle=None, save_as=None):
    """Single-panel vs-Q^2 overlay for one task.

    `sims` is a list of 'sim'-kind Series (0 or more) overlaid on `ref`.
    Intended for fit-derived quantities such as sigma_pcmx/sigma_pcmy.
    """
    if with_ratio:
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    else:
        fig, ax = plt.subplots(1, 1)
        ax_ratio = None

    ax.set_ylabel(ylabel, fontsize=ylabel_size)
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    ref_xy = draw(ax, ref, task_name, selection, axis_bin=axis_bin,
                  q2_panel=True, unit_scale=unit_scale)
    for sim in sims:
        draw(ax, sim, task_name, selection, axis_bin=axis_bin,
             q2_panel=True, unit_scale=unit_scale)

    if with_ratio:
        _draw_ratio_panel(ax_ratio, ref_xy, sims, task_name, selection,
                          axis_bin, offset_scale=1.0,
                          integrated=False, pattern=None,
                          unit_scale=unit_scale)
        ax_ratio.set_xlabel(xlabel, fontsize=15)
        ax_ratio.set_xlim(*xlim)
        ax_ratio.set_ylim(*ratio_ylim)
        ax_ratio.yaxis.set_major_locator(MaxNLocator(nbins=4))
        _drop_extreme_ytick(ax, 'min')
        _drop_extreme_ytick(ax_ratio, 'max')
    else:
        ax.set_xlabel(xlabel, fontsize=15)

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout()
    return _save(pdf, fig, save_as)


def plot_q2_2x2(pdf, ref, sims, task_name, ylabel,
                pmiss_labels, selection='epp',
                xlim=(1.5, 5.0), ylim=(0.0, 0.25),
                label_xy=(2.3, 0.2), ylabel_size=15,
                big_label=None, big_label_xy=(1.7, 0.12),
                unit_scale=False, suptitle=None,
                label_positions=None, save_as=None):
    """2x2 panel of one vs-Q^2 quantity in pMiss bins 0..3.

    `sims` is a list of 'sim'-kind Series (0 or more) overlaid on `ref`.
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
    for bin_idx, ((r, c), label) in enumerate(zip(positions, pmiss_labels)):
        draw(ax[r, c], ref, task_name, selection, axis_bin=[bin_idx],
             q2_panel=True, unit_scale=unit_scale)
        for sim in sims:
            draw(ax[r, c], sim, task_name, selection, axis_bin=[bin_idx],
                   q2_panel=True, unit_scale=unit_scale)
        ax[r, c].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(bin_idx)
        _auto_panel_label(ax[r, c], label, fontsize=13, manual_pos=manual)

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

    return _save(pdf, fig, save_as)


def plot_q2_2x2_data_over_sim(pdf, ref, sims, task_name,
                              pmiss_labels, selection='epp',
                              xlim=(1.5, 5.0), ylim=(0.5, 1.5),
                              label_xy=(2.3, 1.35), ylabel_size=13,
                              big_label=None, big_label_xy=(1.7, 1.12),
                              unit_scale=False, suptitle=None,
                              label_positions=None, save_as=None):
    """2x2 panel of Data/Sim for one vs-Q^2 quantity in pMiss bins 0..3.

    `sims` is a list of 'sim'-kind Series (0 or more); each gets its own
    overlaid trace per panel.
    """
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
            ref, task_name, selection, axis_bin=[bin_idx], unit_scale=unit_scale)
        for sim in sims:
            sim_x, sim_y, sim_err = get_xy_err(
                sim, task_name, selection, axis_bin=[bin_idx], unit_scale=unit_scale)
            if len(sim_x) != len(ref_x):
                raise ValueError(
                    "Cannot build data/sim panel for task '%s' bin %d: "
                    "data and sim have different point counts (%d vs %d)."
                    % (task_name, bin_idx, len(ref_x), len(sim_x)))

            ratio_y, ratio_err = ratio.ratio_series_with_error(ref_y, ref_err, sim_y, sim_err)
            ax[r, c].errorbar(ref_x, ratio_y, ratio_err,
                              color=sim.color, linestyle='', marker='o',
                              markersize=3)
        ax[r, c].axhline(1.0, color='gray', linewidth=0.5, linestyle='--')
        ax[r, c].yaxis.set_major_locator(MaxNLocator(nbins=4))
        manual = None if label_positions is None else label_positions.get(bin_idx)
        _auto_panel_label(ax[r, c], label, fontsize=13, manual_pos=manual)

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

    return _save(pdf, fig, save_as)


def plot_sig_waterfall(pdf, q2_bins, ref, sims, axis='x',
                       xlim=(-0.6, 0.6), ylim=(0.0, 120.0), save_as=None):
    """Legacy-style pseudo-3D waterfall panel used in old plotTOOL.plotSig.

    `sims` is a list of 'sim'-kind Series (0 or more); each is normalized
    and drawn independently. Draws one C.M. histogram-vs-Q2 panel per Q2 bin
    on slanted coordinates. Uses per-Q2 data/sim normalization (sum over y
    bins in each panel) to match the original visual comparison behavior.
    """
    if len(q2_bins) < 2:
        raise ValueError('q2_bins must contain at least two edges')

    q2_start = [0.4, 0.5]
    q2_end = [0.07, 0.02]

    fig = plt.figure()
    ax_back = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax_back.set_ylim(0.0, 1.0)
    ax_back.set_xlim(0.0, 1.0)
    ax_back.plot([q2_start[0], q2_end[0]], [q2_start[1], q2_end[1]], 'k', linewidth=1)
    ax_back.text(0.1, 0.35, r'$Q^{2} [GeV]$', fontsize=15)
    ax_back.text(0.05, 0.9, r'$(e,e^{\prime}pp)$', fontsize=25)

    ticks = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    for t in ticks:
        tx = q2_start[0] + ((t - 1.5) / (5.0 - 1.5)) * (q2_end[0] - q2_start[0])
        ty = q2_start[1] + ((t - 1.5) / (5.0 - 1.5)) * (q2_end[1] - q2_start[1])
        ax_back.plot([tx, tx - 0.01], [ty, ty], 'k', linewidth=1.0)
        ax_back.text(tx - 0.02, ty, str(t), horizontalalignment='right',
                     verticalalignment='center')

    task_name = 'pcm%s_epp_SRC_Q2' % axis
    xlabel = r'$\vec{p}_{C.M.}\cdot\hat{v}_{%s}[GeV]$' % axis

    for i in range(0, len(q2_bins) - 1):
        q2 = 0.5 * (q2_bins[i] + q2_bins[i + 1])
        panel_x = q2_start[0] + ((q2 - 1.5) / (5.0 - 1.5)) * (q2_end[0] - q2_start[0])
        panel_y = q2_start[1] + ((q2 - 1.5) / (5.0 - 1.5)) * (q2_end[1] - q2_start[1])

        ax = fig.add_axes([panel_x, panel_y, 0.7, 0.5])

        dx, dy, de = get_xy_err(ref, task_name, 'epp', axis_bin=[i], unit_scale=True)
        draw(ax, ref, task_name, 'epp', axis_bin=[i], q2_panel=True, unit_scale=True)

        # Use panel-by-panel scaling exactly like the old plotSig's factor[i],
        # computed independently for each sim.
        panel_top = max([y + e for y, e in zip(dy, de)] or [0.0])
        for sim in sims:
            sx, sy, se = get_xy_err(sim, task_name, 'epp', axis_bin=[i], unit_scale=True)
            sim_scale = 1.0
            sum_sim = sum(sy)
            if sum_sim != 0:
                sim_scale = float(sum(dy)) / float(sum_sim)
            sy = [v * sim_scale for v in sy]
            se = [v * sim_scale for v in se]
            pp.line_with_band_q2(ax, sx, sy, se, sim.color)
            sim_top = max([y + e for y, e in zip(sy, se)] or [0.0])
            panel_top = max(panel_top, sim_top)

        # Avoid clipping when a panel's peak exceeds the nominal y-limit.
        y_low = ylim[0]
        y_high = max(ylim[1], panel_top * 1.08 if panel_top > 0.0 else ylim[1])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i != 0:
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])
        if i != (len(q2_bins) - 2):
            ax.set_xticks([])
        ax.patch.set_alpha(0.0)
        ax.set_ylim(y_low, y_high)
        ax.set_xlim(*xlim)
        if i == 0:
            ax.set_ylabel(r'Counts', fontsize=15)
        if i == (len(q2_bins) - 2):
            ax.set_xlabel(xlabel, fontsize=15)

    return _save(pdf, fig, save_as)


def plot_emiss_waterfall(pdf, pmiss_bins, ref, sims, task_name,
                         selection='epp', xlim=(-0.2, 0.4),
                         ylim=(0.0, 120.0), xlabel=r'$E_{miss} [GeV]$',
                         panel_label=None, save_as=None):
    """Legacy-style pseudo-3D waterfall for E_miss vs pMiss.

    `sims` is a list of 'sim'-kind Series (0 or more), each normalized and
    drawn independently. Draws one E_miss panel per pMiss bin on slanted
    coordinates, with pMiss shown on the far-left axis and Counts as the
    in-panel y axis. Each simulation is normalized to data separately in
    each pMiss panel for visual shape comparison, matching the old
    waterfall behavior used for Q2 slices.
    """
    if len(pmiss_bins) < 2:
        raise ValueError('pmiss_bins must contain at least two edges')

    pmiss_start = [0.4, 0.5]
    pmiss_end = [0.07, 0.12]

    fig = plt.figure()
    ax_back = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax_back.set_ylim(0.0, 1.0)
    ax_back.set_xlim(0.0, 1.0)
    ax_back.plot([pmiss_start[0], pmiss_end[0]],
                 [pmiss_start[1], pmiss_end[1]],
                 'k', linewidth=1)
    ax_back.text(0.1, 0.35, r'$p_{miss} [GeV]$', fontsize=15)
    if panel_label is not None:
        ax_back.text(0.05, 0.9, panel_label, fontsize=22)

    # Tick marks on the slanted pMiss guide axis.
    for t in pmiss_bins:
        tx = pmiss_start[0] + ((t - pmiss_bins[0]) / (pmiss_bins[-1] - pmiss_bins[0])) * (pmiss_end[0] - pmiss_start[0])
        ty = pmiss_start[1] + ((t - pmiss_bins[0]) / (pmiss_bins[-1] - pmiss_bins[0])) * (pmiss_end[1] - pmiss_start[1])
        ax_back.plot([tx, tx - 0.01], [ty, ty], 'k', linewidth=1.0)
        ax_back.text(tx - 0.02, ty, ('%.2g' % t), horizontalalignment='right',
                     verticalalignment='center')

    for i in range(0, len(pmiss_bins) - 1):
        pmiss = 0.5 * (pmiss_bins[i] + pmiss_bins[i + 1])
        panel_x = pmiss_start[0] + ((pmiss - pmiss_bins[0]) / (pmiss_bins[-1] - pmiss_bins[0])) * (pmiss_end[0] - pmiss_start[0])
        panel_y = pmiss_start[1] + ((pmiss - pmiss_bins[0]) / (pmiss_bins[-1] - pmiss_bins[0])) * (pmiss_end[1] - pmiss_start[1])

        ax = fig.add_axes([panel_x, panel_y, 0.7, 0.5])

        _, dy, de = get_xy_err(ref, task_name, selection, axis_bin=[i], unit_scale=True)
        draw(ax, ref, task_name, selection, axis_bin=[i],
             q2_panel=True, unit_scale=True)

        # Per-panel data/sim normalization for visual shape comparison,
        # computed independently for each sim.
        panel_top = max([y + e for y, e in zip(dy, de)] or [0.0])
        for sim in sims:
            sx, sy, se = get_xy_err(sim, task_name, selection, axis_bin=[i],
                        unit_scale=True)
            sim_scale = 1.0
            sum_sim = sum(sy)
            if sum_sim != 0:
                sim_scale = float(sum(dy)) / float(sum_sim)
            sy = [v * sim_scale for v in sy]
            se = [v * sim_scale for v in se]
            pp.line_with_band_q2(ax, sx, sy, se, sim.color)
            sim_top = max([y + e for y, e in zip(sy, se)] or [0.0])
            panel_top = max(panel_top, sim_top)

        # Avoid clipping when a panel's peak exceeds the nominal y-limit.
        y_low = ylim[0]
        y_high = max(ylim[1], panel_top * 1.08 if panel_top > 0.0 else ylim[1])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i != 0:
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])
        if i != (len(pmiss_bins) - 2):
            ax.set_xticks([])
        ax.patch.set_alpha(0.0)
        ax.set_ylim(y_low, y_high)
        ax.set_xlim(*xlim)
        if i == 0:
            ax.set_ylabel(r'Counts', fontsize=15)
        if i == (len(pmiss_bins) - 2):
            ax.set_xlabel(xlabel, fontsize=15)

    return _save(pdf, fig, save_as)


# NOTE: the old g_sigma_pcmx / g_sigma_E1miss_*_pmiss / g_mean_E1miss_*_pmiss
# plots (C.M. width vs Q^2, E_miss mean/width vs Q^2) came from Gaussian
# fits across the 100 toy histograms, not from simple bin sums. That data
# was deliberately NOT put into diffTable/integratedTable -- only the
# nominal+toy histograms themselves were persisted (hists/nominal/...,
# hists/toy_NNN/...) for exactly this purpose. Reproducing those plots needs
# a separate fit-extraction step (read each toy's histogram, fit, collect
# mean/sigma across toys) that hasn't been built yet -- intentionally left
# out of plot_data_vs_sim.py / plot_A_dependence.py rather than faked here.
