"""
plot_1D.py

Single-variable 1D distributions, data vs 0-4 simulation models, each
written out as its own single-figure PDF under pdf/1D/ (no combined book --
these are meant to be dropped individually into different note sections).

Variables: pMiss, kMiss (e'p and e'pp), pRel (e'pp), pcm_x/y/z (e'pp),
theta_pmiss (e'p and e'pp), q (e'p and e'pp), theta_pLead,q (e'p and e'pp).

Usage:
    python plot_1D.py Data_He.root
    python plot_1D.py Data_He.root Sim_He_AV18.root --sim-label AV18
    python plot_1D.py Data_He.root Sim_AV18.root Sim_N2LO.root \\
        --sim-label AV18 --sim-label N2LO --sim-color red --sim-color blue
"""

import argparse
import os

import graph_io
import normalization
from plot_helpers import Series
import plot_helpers as ph


# (output key, task_name, selection, xlabel, xlim)
VARIABLES = [
    ('pMiss_ep', 'pMiss_ep_note', 'ep', r'$p_{miss} [GeV]$', (0.4, 1.2)),
    ('pMiss_epp', 'pMiss_epp_note', 'epp', r'$p_{miss} [GeV]$', (0.4, 1.2)),
    ('kMiss_ep', 'kMiss_ep_note', 'ep', r'$k_{miss} [GeV]$', (0.3, 1.1)),
    ('kMiss_epp', 'kMiss_epp_note', 'epp', r'$k_{miss} [GeV]$', (0.3, 1.1)),
    ('pRel_epp', 'pRel_epp', 'epp', r'$p_{rel} [GeV]$', (0.15, 1.0)),
    ('pcmz_epp', 'pcmz_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{z} [GeV]$', (-0.75, 0.75)),
    ('pcmx_epp', 'pcmx_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{x} [GeV]$', (-0.75, 0.75)),
    ('pcmy_epp', 'pcmy_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{y} [GeV]$', (-0.75, 0.75)),
    ('theta_pmiss_ep', 'theta_pmiss_ep', 'ep', r'$\theta_{p_{miss},q} [deg]$', (100, 180)),
    ('theta_pmiss_epp', 'theta_pmiss_epp', 'epp', r'$\theta_{p_{miss},q} [deg]$', (100, 180)),
    ('q_ep', 'q_ep', 'ep', r'$|\vec{q}| [GeV]$', (0, 4)),
    ('q_epp', 'q_epp', 'epp', r'$|\vec{q}| [GeV]$', (0, 4)),
    ('theta_pLeadq_ep', 'theta_pLeadq_ep', 'ep', r'$\theta_{p_{lead},q} [deg]$', (0, 40)),
    ('theta_pLeadq_epp', 'theta_pLeadq_epp', 'epp', r'$\theta_{p_{lead},q} [deg]$', (0, 40)),
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
                   help='0 to 4 Main_Figs_Binned+BuildGraphs sim outputs')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-label', action='append', default=[])
    p.add_argument('--sim-color', action='append', default=[])
    p.add_argument('--out-dir', default='pdf/1D')
    p.add_argument('--normalization-mode',
                   choices=['legacy-last-q2', 'integrated'],
                   default='legacy-last-q2')
    return p.parse_args()


def _build_sims(args, f_data):
    if len(args.sim_files) > 4:
        raise ValueError('plot_1D.py supports at most 4 simulation files, got %d'
                         % len(args.sim_files))
    sims = []
    for i, sim_file in enumerate(args.sim_files):
        label = args.sim_label[i] if i < len(args.sim_label) else 'Sim %d' % (i + 1)
        color = args.sim_color[i] if i < len(args.sim_color) else \
            DEFAULT_SIM_COLORS[i % len(DEFAULT_SIM_COLORS)]
        f_sim = graph_io.open_file(sim_file)

        scale_ep, _ = normalization.scale_factor(
            f_sim, f_data, 'ep', mode=args.normalization_mode)

        # epp events are normalized to the same e'p-event-count scale factor
        # as ep events, rather than their own independent epp-event count.
        # To go back to normalizing epp plots by their own epp yield, comment
        # the line below and uncomment the one above it.
        # scale_epp, _ = normalization.scale_factor(
        #     f_sim, f_data, 'epp', mode=args.normalization_mode)
        scale_epp = scale_ep

        sims.append(Series(label, f_sim, color, 'sim', scale_ep, scale_epp))
    return sims


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    f_data = graph_io.open_file(args.data_file)
    data = Series(args.data_label, f_data, args.data_color, 'data')
    sims = _build_sims(args, f_data)

    pdf = _NoOpPdf()
    missing = []
    for key, task_name, selection, xlabel, xlim in VARIABLES:
        try:
            _, y, yerr = ph.get_xy_err(data, task_name, selection)
        except Exception as exc:
            missing.append('%s: %s' % (key, exc))
            continue
        peak = max([v + e for v, e in zip(y, yerr)] or [0.0])
        ylim = (0.0, peak * 1.2 if peak > 0.0 else 1.0)

        out_path = os.path.join(args.out_dir, key + '.pdf')
        ph.plot_overlay(
            pdf, task_name, [data] + sims, selection=selection,
            xlabel=xlabel, ylabel='Counts', xlim=xlim, ylim=ylim,
            save_as=out_path)
        print('Wrote %s' % out_path)

    if missing:
        print('Skipped variables with missing graphs:')
        for m in missing:
            print(' - %s' % m)


if __name__ == '__main__':
    main()
