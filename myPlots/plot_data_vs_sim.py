"""
plot_data_vs_sim.py

Data vs a single simulation model, every variable, each overlay with a
data/sim ratio subpanel underneath (plot_overlay's with_ratio=True default).

Takes exactly two input files (both already run through Main_Figs_Binned
and BuildGraphs) -- meant to be copied off the cluster and run locally
without also needing nuclei_config.py's full multi-nucleus/multi-model file
list. For the 5-model He comparison or the cross-nucleus comparison, see
plot_A_dependence.py / re-add a loop here if you do want several models at
once -- this script intentionally only handles one data file + one sim file.

Usage:
    python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root
    python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root \\
        --data-label '${}^{4}He$' --sim-label 'AV18' --out plots/he_av18.pdf

NOTE: the old C.M.-width-vs-Q^2 and E_miss-mean/width-vs-Q^2 plots
(g_sigma_pcmx, g_sigma_E1miss_ep_pmiss, g_mean_E1miss_ep_pmiss, ...) are not
reproduced here -- they came from Gaussian fits across the 100 toy
histograms, which is a separate fit-extraction step not yet built on top of
the persisted hists/nominal + hists/toy_NNN histograms. See the note in
plot_helpers.py.
"""

import argparse
import matplotlib.backends.backend_pdf

import graph_io
import normalization
from labels import PMISS_LABELS
from plot_helpers import Series
import plot_helpers as ph


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('data_file', help='Main_Figs_Binned+BuildGraphs output for data')
    p.add_argument('sim_file', help='Main_Figs_Binned+BuildGraphs output for one sim model')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--sim-label', default='Sim')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-color', default='red')
    p.add_argument('--out', default='plots/data_vs_sim.pdf')
    return p.parse_args()


def main():
    args = parse_args()
    pdf = matplotlib.backends.backend_pdf.PdfPages(args.out)

    f_data = graph_io.open_file(args.data_file)
    f_sim = graph_io.open_file(args.sim_file)

    data = Series(args.data_label, f_data, args.data_color, 'data')

    scale_ep, _ = normalization.scale_factor(f_sim, f_data, 'ep')
    scale_epp, _ = normalization.scale_factor(f_sim, f_data, 'epp')
    sim = Series(args.sim_label, f_sim, args.sim_color, 'sim', scale_ep, scale_epp)

    # ---------------------------------------------------------------------
    # C.M. projections
    # ---------------------------------------------------------------------
    for axis in ('x', 'y'):
        ph.plot_overlay(
            pdf, 'pcm{0}_epp'.format(axis), [data, sim], selection='epp',
            xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
            ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 350))

    # ---------------------------------------------------------------------
    # epp/ep vs pMiss (already a ratio quantity -- no rescaling)
    # ---------------------------------------------------------------------
    ph.plot_overlay(
        pdf, 'pMiss_epp_over_pMiss_ep', [data, sim], selection='ratio',
        xlabel=r'$p_{miss} [GeV]$', ylabel=r'$epp/ep$',
        xlim=(0.4, 1.0), ylim=(0, 0.2), with_ratio=True,
        xlabel_size=25, ylabel_size=25)

    # ---------------------------------------------------------------------
    # E1miss 4x2 (ep | epp), E2miss 4x1 (epp)
    # ---------------------------------------------------------------------
    ph.plot_emiss_4x2(
        pdf, [data, sim], r'$E_{1,miss} [GeV]$',
        'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
        ylim_ep=[1000, 1000, 850, 600], ylim_epp=[60, 110, 110, 60],
        labely_ep=[850, 850, 700, 500], labely_epp=[50, 95, 95, 50],
        pmiss_labels=PMISS_LABELS)

    ph.plot_emiss_4x1(
        pdf, [data, sim], r'$E_{2,miss} [GeV]$',
        'E2miss_epp_SRC_pmiss',
        ylim=[60, 110, 110, 80], labely=[50, 95, 95, 65],
        pmiss_labels=PMISS_LABELS)

    # ---------------------------------------------------------------------
    # Ratio-only pages for paneled plots
    # ---------------------------------------------------------------------
    ph.plot_emiss_4x2_ratio_only(
        pdf, data, sim, r'$E_{1,miss} [GeV]$',
        'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
        pmiss_labels=PMISS_LABELS,
        xlim=(-0.15, 0.4), ylim=(0.5, 1.5), labely=1.35)

    ph.plot_emiss_4x1_ratio_only(
        pdf, data, sim, r'$E_{2,miss} [GeV]$',
        'E2miss_epp_SRC_pmiss', pmiss_labels=PMISS_LABELS,
        xlim=(-0.2, 0.4), ylim=(0.5, 1.5), labely=1.35)

    # ---------------------------------------------------------------------
    # epp/ep vs Q^2 in pMiss bins
    # ---------------------------------------------------------------------
    ph.plot_q2_2x2_ratio(
        pdf, data, sim, 'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', r'$epp/ep$',
        PMISS_LABELS, ylim=(0, 0.25), label_xy=(2.3, 0.2), ylabel_size=15)

    ph.plot_q2_2x2_data_over_sim_ratio_only(
        pdf, data, sim, 'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', PMISS_LABELS,
        ylim=(0.5, 1.5), label_xy=(2.3, 1.35), ylabel_size=13)

    pdf.close()
    print('Wrote %s' % args.out)


if __name__ == '__main__':
    main()
