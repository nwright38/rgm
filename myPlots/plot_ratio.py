"""
plot_ratio.py

(e,e'pp)/(e,e'p) ratio vs pMiss and vs kMiss, written out as individual
single-figure PDFs under pdf/ratio/.

Uses the pMiss_epp_over_pMiss_ep / kMiss_epp_over_kMiss_ep ratio tasks
Main_Figs_Binned.cpp already computes from the long-binned pMiss_ep/epp and
kMiss_ep/epp tasks (the "old bins" -- not the uniform-binned *_note variants
added for the 1D plots, which are deliberately kept separate so changing
their binning doesn't also rebin this ratio).

Usage:
    python plot_ratio.py Data_He.root
    python plot_ratio.py Data_He.root Sim_He_AV18.root --sim-label AV18
    python plot_ratio.py Data_He.root Sim_AV18.root Sim_N2LO.root \\
        --sim-label AV18 --sim-label N2LO --sim-color red --sim-color blue
"""

import argparse
import os

import graph_io
from plot_helpers import Series
import plot_helpers as ph


# (output key, task_name, xlabel, xlim)
VARIABLES = [
    ('pMiss_ratio', 'pMiss_epp_over_pMiss_ep', r'$p_{miss} [GeV]$', (0.4, 1.)),
    ('kMiss_ratio', 'kMiss_epp_over_kMiss_ep', r'$k_{miss} [GeV]$', (0.3, .9)),
]

DEFAULT_SIM_COLORS = ['red', 'blue', 'green', 'darkorange']


class _NoOpPdf(object):
    """Stands in for a PdfPages book when we only want save_as's individual
    file, not a combined multi-page PDF."""
    def savefig(self, fig):
        pass


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('data_file', help='Main_Figs_Binned+BuildGraphs output for data')
    p.add_argument('sim_files', nargs='*',
                   help='0 to 4 Main_Figs_Binned+BuildGraphs sim outputs (theory overlay)')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-label', action='append', default=[])
    p.add_argument('--sim-color', action='append', default=[])
    p.add_argument('--out-dir', default='pdf/ratio')
    return p.parse_args()


def _build_sims(args):
    if len(args.sim_files) > 4:
        raise ValueError('plot_ratio.py supports at most 4 simulation files, got %d'
                         % len(args.sim_files))
    sims = []
    for i, sim_file in enumerate(args.sim_files):
        label = args.sim_label[i] if i < len(args.sim_label) else 'Sim %d' % (i + 1)
        color = args.sim_color[i] if i < len(args.sim_color) else \
            DEFAULT_SIM_COLORS[i % len(DEFAULT_SIM_COLORS)]
        # No normalization scale factors: a ratio task's selection is
        # 'ratio', and plot_helpers._scale_for always returns 1.0 for that
        # selection regardless of what's passed here (the epp/ep ratio is
        # already self-normalized, independent of luminosity scaling).
        sims.append(Series(label, graph_io.open_file(sim_file), color, 'sim'))
    return sims


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    f_data = graph_io.open_file(args.data_file)
    data = Series(args.data_label, f_data, args.data_color, 'data')
    sims = _build_sims(args)

    pdf = _NoOpPdf()
    missing = []
    for key, task_name, xlabel, xlim in VARIABLES:
        try:
            _, y, yerr = ph.get_xy_err(data, task_name, 'ratio')
        except Exception as exc:
            missing.append('%s: %s' % (key, exc))
            continue
        peak = max([v + e for v, e in zip(y, yerr)] or [0.0])
        ylim = (0.0, peak * 1.2 if (peak > 0.0 and peak * 1.2 < 1.) else 1.0)
        ylim = (0.0,.3)

        out_path = os.path.join(args.out_dir, key + '.pdf')
        ph.plot_overlay(
            pdf, task_name, [data] + sims, selection='ratio',
            xlabel=xlabel, ylabel=r"$(e,e^{\prime}pp)/(e,e^{\prime}p)$",
            xlim=xlim, ylim=ylim, with_ratio=False, save_as=out_path)
        print('Wrote %s' % out_path)

    if missing:
        print('Skipped variables with missing graphs:')
        for m in missing:
            print(' - %s' % m)


if __name__ == '__main__':
    main()
