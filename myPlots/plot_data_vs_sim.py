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
(sigma_pcmx, stddev_E1miss_ep_pmiss, mean_E1miss_ep_pmiss, ...) are included
when those fit-derived graphs are present in the file's "graphs" directory
(written by ExtractFitQuantities.C after Main_Figs_Binned + BuildGraphs).
"""

import argparse
import matplotlib.backends.backend_pdf

import graph_io
import normalization
from labels import PMISS_LABELS
from plot_helpers import Series
import plot_helpers as ph


Q2_BINS = [1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0]
PMISS_BINS = [0.4, 0.55, 0.7, 0.85, 1.0]


def _require_graphs(series_list, needed):
    """Checks that every (task, selection, axis_bin) graph exists for every
    Series in series_list. Returns (ok, missing_messages)."""
    missing = []
    for s in series_list:
        for task_name, selection, axis_bin in needed:
            try:
                ph.get_xy_err(s, task_name, selection, axis_bin=axis_bin)
            except Exception as exc:
                missing.append(
                    "%s missing graph (%s, %s, %s): %s"
                    % (s.label, task_name, selection, list(axis_bin), exc)
                )
    return (len(missing) == 0), missing


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

    missing_notes = []

    # ---------------------------------------------------------------------
    # C.M. projections
    # ---------------------------------------------------------------------
    for axis in ('x', 'y'):
        ph.plot_overlay(
            pdf, 'pcm{0}_epp'.format(axis), [data, sim], selection='epp',
            xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
            ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 350))

    # Legacy-style pseudo-3D Q2 waterfall pages for C.M. distributions.
    ph.plot_sig_waterfall(pdf, Q2_BINS, data, sim, axis='x', xlim=(-0.6, 0.6), ylim=(0.0, 120.0))
    ph.plot_sig_waterfall(pdf, Q2_BINS, data, sim, axis='y', xlim=(-0.6, 0.6), ylim=(0.0, 120.0))

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

    # Legacy-style pseudo-3D pMiss waterfall pages for E_miss distributions.
    ph.plot_emiss_waterfall(
        pdf, PMISS_BINS, data, sim,
        task_name='E1miss_epp_SRC_pmiss', selection='epp',
        xlim=(-0.15, 0.4), ylim=(0.0, 120.0),
        xlabel=r'$E_{1,miss} [GeV]$',
        panel_label=r'$(e,e^{\prime}pp)$')

    ep_waterfall_needed = [
        ('E1miss_ep_SRC_pmiss', 'ep', (0,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (1,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (2,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (3,)),
    ]
    ok, miss = _require_graphs([data, sim], ep_waterfall_needed)
    if ok:
        ph.plot_emiss_waterfall(
            pdf, PMISS_BINS, data, sim,
            task_name='E1miss_ep_SRC_pmiss', selection='ep',
            xlim=(-0.15, 0.4), ylim=(0.0, 1000.0),
            xlabel=r'$E_{1,miss} [GeV]$',
            panel_label=r'$(e,e^{\prime}p)$')
    else:
        missing_notes.append(
            "Skipped pseudo-3D E1miss vs pMiss page for (e,e'p) because required graphs were missing.")
        missing_notes.extend(miss)

    ph.plot_emiss_waterfall(
        pdf, PMISS_BINS, data, sim,
        task_name='E2miss_epp_SRC_pmiss', selection='epp',
        xlim=(-0.2, 0.4), ylim=(0.0, 120.0),
        xlabel=r'$E_{2,miss} [GeV]$',
        panel_label=r'$(e,e^{\prime}pp)$')

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

    # ---------------------------------------------------------------------
    # Fit-derived quantities (from ExtractFitQuantities.C)
    # ---------------------------------------------------------------------
    fit_q2_single = [
        ('sigma_pcmx', 'epp', ()),
        ('sigma_pcmy', 'epp', ()),
    ]
    ok, miss = _require_graphs([data, sim], fit_q2_single)
    if ok:
        ph.plot_q2_single(
            pdf, data, sim, 'sigma_pcmx', r'$\sigma_{x,C.M.} [GeV]$',
            selection='epp', xlim=(1.5, 5.0), with_ratio=True,
            unit_scale=True,
            suptitle='Extracted quantity: sigma_pcmx (overlay + Data/Sim)')
        ph.plot_q2_single(
            pdf, data, sim, 'sigma_pcmy', r'$\sigma_{y,C.M.} [GeV]$',
            selection='epp', xlim=(1.5, 5.0), with_ratio=True,
            unit_scale=True,
            suptitle='Extracted quantity: sigma_pcmy (overlay + Data/Sim)')
    else:
        missing_notes.append(
            'Skipped sigma_{x/y,C.M.} vs Q^2 pages because required extracted graphs were missing.')
        missing_notes.extend(miss)

    fit_q2_panels = [
        ('stddev_E1miss_ep_pmiss', 'ep',
         r'Std Dev $E_{1,miss} [GeV]$', (0, 0.15), (2.3, 0.01), 10, r'$(e,e^{\prime}p)$'),
        ('stddev_E1miss_epp_pmiss', 'epp',
         r'Std Dev $E_{1,miss} [GeV]$', (0, 0.15), (2.3, 0.01), 10, r'$(e,e^{\prime}pp)$'),
        ('stddev_E2miss_epp_pmiss', 'epp',
         r'Std Dev $E_{2,miss} [GeV]$', (0, 0.15), (2.3, 0.01), 10, r'$(e,e^{\prime}pp)$'),
        ('mean_E1miss_ep_pmiss', 'ep',
         r'Mean $E_{1,miss} [GeV]$', (-0.02, 0.3), (2.3, 0.04), 10, r'$(e,e^{\prime}p)$'),
        ('mean_E1miss_epp_pmiss', 'epp',
         r'Mean $E_{1,miss} [GeV]$', (-0.02, 0.3), (2.3, 0.04), 10, r'$(e,e^{\prime}pp)$'),
        ('mean_E2miss_epp_pmiss', 'epp',
         r'Mean $E_{2,miss} [GeV]$', (-0.02, 0.3), (2.3, 0.04), 10, r'$(e,e^{\prime}pp)$'),
    ]

    for task_name, selection, ylabel, ylim, label_xy, ylabel_size, big_label in fit_q2_panels:
        needed = [(task_name, selection, (0,)),
                  (task_name, selection, (1,)),
                  (task_name, selection, (2,)),
                  (task_name, selection, (3,))]
        ok, miss = _require_graphs([data, sim], needed)
        if not ok:
            missing_notes.append(
                'Skipped panel for %s because required extracted graphs were missing.' % task_name)
            missing_notes.extend(miss)
            continue
        ph.plot_q2_2x2(
            pdf, data, sim, task_name, ylabel,
            PMISS_LABELS, selection=selection,
            ylim=ylim, label_xy=label_xy, ylabel_size=ylabel_size,
            big_label=big_label, unit_scale=True,
            suptitle='Extracted quantity: %s (overlay)' % task_name)

        ph.plot_q2_2x2_data_over_sim(
            pdf, data, sim, task_name,
            PMISS_LABELS, selection=selection,
            ylim=(0.5, 1.5), label_xy=(2.3, 1.35), ylabel_size=13,
            unit_scale=True,
            suptitle='Extracted quantity: %s (Data/Sim)' % task_name)

    pdf.close()
    if missing_notes:
        print('Some requested helium plots were skipped due to missing graphs:')
        for m in missing_notes:
            print(' - %s' % m)
    print('Wrote %s' % args.out)


if __name__ == '__main__':
    main()
