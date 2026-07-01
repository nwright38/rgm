"""
Scratch p_cm,z vs Q2-bin overlays.

Reads Main_Figs_Binned output after BuildGraphs has been run and overlays the
pcmz_epp_SRC_Q2 distributions: one panel per established Q2 bin. Simulation
curves are scaled independently in each Q2 bin to the data yield in that same
panel, so the figure emphasizes shape differences.

Usage:
    python myPlots/plot_pcmz_q2_scratch.py Data.root Sim.root \\
        --sim-label GCF --out myPlots/pdf/pcmz_q2_scratch.pdf
"""

from __future__ import division

import argparse
import math
import os
import tempfile

# Keep matplotlib/fontconfig from trying to write cache files under $HOME on
# systems where that directory is read-only or unavailable in batch jobs.
os.environ.setdefault('MPLCONFIGDIR', os.path.join(tempfile.gettempdir(), 'matplotlib-cache'))
os.environ.setdefault('XDG_CACHE_HOME', os.path.join(tempfile.gettempdir(), 'xdg-cache'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import graph_io
import graph_names
import plot_primitives as pp


Q2_BINS = [1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0]
TASK_NAME = 'pcmx_epp_SRC_Q2'
SELECTION = 'epp'
DEFAULT_SIM_COLORS = ['red', 'blue', 'green', 'darkorange', 'purple', 'brown']
Q2_COLORS = ['black', 'red', 'blue', 'green', 'darkorange', 'purple', 'brown']
Q2_MARKERS = ['o', 's', '^', 'v', 'D', 'P', 'X']
DATA_OVERLAY_X_OFFSET = 0.006


class Series(object):
    def __init__(self, label, root_file, color, kind):
        self.label = label
        self.file = root_file
        self.color = color
        self.kind = kind


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('data_file', help='Data ROOT file after Main_Figs_Binned + BuildGraphs')
    p.add_argument('sim_files', nargs='*',
                   help='Simulation ROOT files after Main_Figs_Binned + BuildGraphs')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-label', action='append', default=[])
    p.add_argument('--sim-color', action='append', default=[])
    p.add_argument('--out', default='myPlots/pdf/pcmz_q2_scratch.pdf')
    p.add_argument('--xlim', nargs=2, type=float, default=(-0.5, .5))
    p.add_argument('--unit-area', action='store_true',
                   help='Normalize data and sim to unit area in each Q2 bin instead of scaling sim to data yield.')
    return p.parse_args()


def read_panel(series, q2_bin):
    name = graph_names.diff_graph_name(TASK_NAME, SELECTION, [q2_bin])
    return graph_io.read_graph_asymm(series.file, name)


def scale_graph(y, yerr_low, yerr_high, scale):
    return (
        [v * scale for v in y],
        [v * scale for v in yerr_low],
        [v * scale for v in yerr_high],
    )


def graph_sum(y):
    return float(sum(y))


def panel_scale(reference_y, target_y, unit_area=False):
    target_sum = graph_sum(target_y)
    if target_sum == 0.0:
        return 1.0
    if unit_area:
        return 1.0 / target_sum
    return graph_sum(reference_y) / target_sum


def q2_label(i):
    return r'$%.2g<Q^{2}<%.2g$' % (Q2_BINS[i], Q2_BINS[i + 1])


def shifted_x(x, bin_index, n_bins):
    center = 0.5 * (n_bins - 1)
    offset = (bin_index - center) * DATA_OVERLAY_X_OFFSET
    return [xi + offset for xi in x]


def build_sims(args):
    sims = []
    for i, path in enumerate(args.sim_files):
        label = args.sim_label[i] if i < len(args.sim_label) else 'Sim %d' % (i + 1)
        color = args.sim_color[i] if i < len(args.sim_color) else \
            DEFAULT_SIM_COLORS[i % len(DEFAULT_SIM_COLORS)]
        sims.append(Series(label, graph_io.open_file(path), color, 'sim'))
    return sims


def draw_panel(ax, data, sims, q2_bin, args):
    dx, dy, de_low, de_high = read_panel(data, q2_bin)

    data_scale = panel_scale(dy, dy, unit_area=args.unit_area)
    dy, de_low, de_high = scale_graph(dy, de_low, de_high, data_scale)
    pp.step_with_error(ax, dx, dy, de_low, de_high, data.color)

    ymax = max([y + e for y, e in zip(dy, de_high)] or [0.0])
    for sim in sims:
        sx, sy, se_low, se_high = read_panel(sim, q2_bin)
        scale = panel_scale(dy, sy, unit_area=args.unit_area)
        sy, se_low, se_high = scale_graph(sy, se_low, se_high, scale)
        pp.line_with_band(ax, sx, sy, se_low, se_high, sim.color, alpha=0.22)
        ymax = max(ymax, max([y + e for y, e in zip(sy, se_high)] or [0.0]))

    ax.set_xlim(*args.xlim)
    ax.set_ylim(0.0, ymax * 1.18 if ymax > 0.0 else 1.0)
    ax.text(0.04, 0.92, q2_label(q2_bin), transform=ax.transAxes,
            ha='left', va='top', fontsize=10)
    ax.tick_params(labelsize=9)


def draw_all_data_overlay(data, args):
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    ymax = 0.0
    handles = []

    n_q2_bins = len(Q2_BINS) - 1
    for i in range(n_q2_bins):
        x, y, yerr_low, yerr_high = read_panel(data, i)
        scale = panel_scale(y, y, unit_area=True)
        y, yerr_low, yerr_high = scale_graph(y, yerr_low, yerr_high, scale)

        color = Q2_COLORS[i % len(Q2_COLORS)]
        marker = Q2_MARKERS[i % len(Q2_MARKERS)]
        x_plot = shifted_x(x, i, n_q2_bins)
        ax.errorbar(x_plot, y, yerr=[yerr_low, yerr_high], color=color,
                    linestyle='', marker=marker, markersize=4.5,
                    capsize=2.0, linewidth=1.0, alpha=0.9)
        handles.append(plt.Line2D([0], [0], color=color, marker=marker,
                                  linestyle='', markersize=5,
                                  label=q2_label(i)))
        ymax = max(ymax, max([v + e for v, e in zip(y, yerr_high)] or [0.0]))

    ax.set_xlim(*args.xlim)
    ax.set_ylim(0.0, ymax * 1.18 if ymax > 0.0 else 1.0)
    ax.set_xlabel(r'$\vec{p}_{C.M.}\cdot\hat{v}_{y}$ [GeV]', fontsize=13)
    ax.set_ylabel('Unit area', fontsize=13)
    ax.set_title(r'%s: $p_{C.M.,y}$ data slices by $Q^{2}$' % args.data_label,
                 fontsize=14)
    ax.legend(handles=handles, loc='best', frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


def draw_data_sim_panels(data, sims, args):
    n_panels = len(Q2_BINS) - 1
    n_cols = 4
    n_rows = int(math.ceil(n_panels / float(n_cols)))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(13.0, 6.8), sharex=True)
    axes = list(axes.flat)

    for i in range(n_panels):
        draw_panel(axes[i], data, sims, i, args)

    for ax in axes[n_panels:]:
        ax.axis('off')

    for i, ax in enumerate(axes[:n_panels]):
        if i % n_cols == 0:
            ax.set_ylabel('Unit area' if args.unit_area else 'Counts', fontsize=11)
        if i >= n_panels - n_cols:
            ax.set_xlabel(r'$\vec{p}_{C.M.}\cdot\hat{v}_{y}$ [GeV]', fontsize=11)

    handles = [
        plt.Line2D([0], [0], color=data.color, linewidth=2, label=data.label)
    ]
    handles += [
        plt.Line2D([0], [0], color=sim.color, linewidth=2, label=sim.label)
        for sim in sims
    ]
    fig.legend(handles=handles, loc='upper center', ncol=max(1, len(handles)),
               frameon=False, fontsize=11)
    fig.suptitle(r'$p_{C.M.,y}$ in established $Q^{2}$ bins, $(e,e^\prime pp)$',
                 y=0.98, fontsize=14)
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.92])
    return fig


def main():
    args = parse_args()
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    data = Series(args.data_label, graph_io.open_file(args.data_file),
                  args.data_color, 'data')
    sims = build_sims(args)

    with PdfPages(args.out) as pdf:
        fig = draw_all_data_overlay(data, args)
        pdf.savefig(fig)
        plt.close(fig)

        fig = draw_data_sim_panels(data, sims, args)
        pdf.savefig(fig)
        plt.close(fig)

    print('Wrote %s' % args.out)


if __name__ == '__main__':
    main()
