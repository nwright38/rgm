"""
plot_data_vs_sim.py

Data vs zero, one, or several simulation models, every variable, each
overlay with a data/sim ratio subpanel underneath (plot_overlay's
with_ratio=True default).

Takes one data file and any number of sim files (all already run through
Main_Figs_Binned and BuildGraphs) -- meant to be copied off the cluster and
run locally without also needing nuclei_config.py's full multi-nucleus/
multi-model file list. For the cross-nucleus comparison, see
plot_A_dependence.py.

Usage:
    # Data only, no simulation overlay.
    python plot_data_vs_sim.py Data_He.root

    # Data vs one simulation.
    python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root \\
        --sim-label 'AV18' --sim-color red

    # Data vs several simulations, overlaid together.
    python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root Sim_He_N2LO.root \\
        --sim-label 'AV18' --sim-label 'N2LO' \\
        --sim-color red --sim-color blue \\
        --out plots/he_multi.pdf

NOTE: the old C.M.-width-vs-Q^2 and E_miss-mean/width-vs-Q^2 plots
(sigma_pcmx, stddev_E1miss_ep_pmiss, mean_E1miss_ep_pmiss, ...) are included
when those fit-derived graphs are present in the file's "graphs" directory
(written by ExtractFitQuantities.C after Main_Figs_Binned + BuildGraphs).

SAVE_AS lets you flag individual plots to also be written out as their own
single-figure PDFs (in addition to the combined --out book), without having
to restructure any of the calls below. Add a 'plot_key': 'path/to/file.pdf'
entry and the matching call will pick it up automatically.
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

DEFAULT_SIM_COLORS = ['red', 'blue', 'green', 'darkorange', 'purple', 'brown']

# Flag a plot for individual saving by adding its key here, e.g.:
#   SAVE_AS = {'pcmx_epp': 'plots/individual/pcmx_epp.pdf'}
# Keys are named at each ph.plot_* call site below via SAVE_AS.get('key').
SAVE_AS = {}


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
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('data_file', help='Main_Figs_Binned+BuildGraphs output for data')
    p.add_argument('sim_files', nargs='*',
                   help='Zero or more Main_Figs_Binned+BuildGraphs sim outputs '
                        '(omit for data-only plots, pass several for a '
                        'multi-model overlay)')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-label', action='append', default=[],
                   help='One per --sim-label, matched in order to sim_files. '
                        'Defaults to "Sim 1", "Sim 2", ... if not given.')
    p.add_argument('--sim-color', action='append', default=[],
                   help='One per --sim-color, matched in order to sim_files. '
                        'Defaults to a fixed color cycle if not given.')
    p.add_argument('--out', default='plots/data_vs_sim.pdf')
    return p.parse_args()


def _build_sims(args, f_data):
    """Opens each sim file and builds its Series, normalized to f_data."""
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
    pdf = matplotlib.backends.backend_pdf.PdfPages(args.out)

    f_data = graph_io.open_file(args.data_file)
    data = Series(args.data_label, f_data, args.data_color, 'data')
    sims = _build_sims(args, f_data)
    has_sim = len(sims) > 0

    missing_notes = []

    # ---------------------------------------------------------------------
    # C.M. projections
    # ---------------------------------------------------------------------
    for axis in ('x', 'y'):
        ph.plot_overlay(
            pdf, 'pcm{0}_epp'.format(axis), [data] + sims, selection='epp',
            xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
            ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 350),
            save_as=SAVE_AS.get('pcm%s_epp' % axis))

    # Legacy-style pseudo-3D Q2 waterfall pages for C.M. distributions.
    ph.plot_sig_waterfall(pdf, Q2_BINS, data, sims, axis='x', xlim=(-0.6, 0.6), ylim=(0.0, 120.0),
                          save_as=SAVE_AS.get('sig_waterfall_x'))
    ph.plot_sig_waterfall(pdf, Q2_BINS, data, sims, axis='y', xlim=(-0.6, 0.6), ylim=(0.0, 120.0),
                          save_as=SAVE_AS.get('sig_waterfall_y'))

    # ---------------------------------------------------------------------
    # epp/ep vs pMiss (already a ratio quantity -- no rescaling)
    # ---------------------------------------------------------------------
    ph.plot_overlay(
        pdf, 'pMiss_epp_over_pMiss_ep', [data] + sims, selection='ratio',
        xlabel=r'$p_{miss} [GeV]$', ylabel=r'$epp/ep$',
        xlim=(0.4, 1.0), ylim=(0, 0.25), with_ratio=True,
        xlabel_size=25, ylabel_size=25,
        save_as=SAVE_AS.get('pMiss_epp_over_pMiss_ep'))

    # ---------------------------------------------------------------------
    # E1miss 4x2 (ep | epp), E2miss 4x1 (epp)
    # ---------------------------------------------------------------------
    ph.plot_emiss_4x2(
        pdf, [data] + sims, r'$E_{1,miss} [GeV]$',
        'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
        ylim_ep=[1000, 1000, 850, 600], ylim_epp=[60, 110, 110, 60],
        labely_ep=[850, 850, 700, 500], labely_epp=[50, 95, 95, 50],
        pmiss_labels=PMISS_LABELS, save_as=SAVE_AS.get('E1miss_4x2'))

    ph.plot_emiss_4x1(
        pdf, [data] + sims, r'$E_{2,miss} [GeV]$',
        'E2miss_epp_SRC_pmiss',
        ylim=[60, 110, 110, 80], labely=[50, 95, 95, 65],
        pmiss_labels=PMISS_LABELS, save_as=SAVE_AS.get('E2miss_4x1'))

    # Legacy-style pseudo-3D pMiss waterfall pages for E_miss distributions.
    ph.plot_emiss_waterfall(
        pdf, PMISS_BINS, data, sims,
        task_name='E1miss_epp_SRC_pmiss', selection='epp',
        xlim=(-0.15, 0.4), ylim=(0.0, 120.0),
        xlabel=r'$E_{1,miss} [GeV]$',
        panel_label=r'$(e,e^{\prime}pp)$',
        save_as=SAVE_AS.get('E1miss_waterfall_epp'))

    ep_waterfall_needed = [
        ('E1miss_ep_SRC_pmiss', 'ep', (0,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (1,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (2,)),
        ('E1miss_ep_SRC_pmiss', 'ep', (3,)),
    ]
    ok, miss = _require_graphs([data] + sims, ep_waterfall_needed)
    if ok:
        ph.plot_emiss_waterfall(
            pdf, PMISS_BINS, data, sims,
            task_name='E1miss_ep_SRC_pmiss', selection='ep',
            xlim=(-0.15, 0.4), ylim=(0.0, 1000.0),
            xlabel=r'$E_{1,miss} [GeV]$',
            panel_label=r'$(e,e^{\prime}p)$',
            save_as=SAVE_AS.get('E1miss_waterfall_ep'))
    else:
        missing_notes.append(
            "Skipped pseudo-3D E1miss vs pMiss page for (e,e'p) because required graphs were missing.")
        missing_notes.extend(miss)

    ph.plot_emiss_waterfall(
        pdf, PMISS_BINS, data, sims,
        task_name='E2miss_epp_SRC_pmiss', selection='epp',
        xlim=(-0.2, 0.4), ylim=(0.0, 120.0),
        xlabel=r'$E_{2,miss} [GeV]$',
        panel_label=r'$(e,e^{\prime}pp)$',
        save_as=SAVE_AS.get('E2miss_waterfall_epp'))

    # ---------------------------------------------------------------------
    # Ratio-only pages for paneled plots (need at least one sim to compare)
    # ---------------------------------------------------------------------
    if has_sim:
        ph.plot_emiss_4x2_ratio_only(
            pdf, data, sims, r'$E_{1,miss} [GeV]$',
            'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
            pmiss_labels=PMISS_LABELS,
            xlim=(-0.15, 0.4), ylim=(0.5, 1.5), labely=1.35,
            save_as=SAVE_AS.get('E1miss_4x2_ratio_only'))

        ph.plot_emiss_4x1_ratio_only(
            pdf, data, sims, r'$E_{2,miss} [GeV]$',
            'E2miss_epp_SRC_pmiss', pmiss_labels=PMISS_LABELS,
            xlim=(-0.2, 0.4), ylim=(0.5, 1.5), labely=1.35,
            save_as=SAVE_AS.get('E2miss_4x1_ratio_only'))

    # ---------------------------------------------------------------------
    # epp/ep vs Q^2 in pMiss bins
    # ---------------------------------------------------------------------
    ph.plot_q2_2x2_ratio(
        pdf, data, sims, 'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', r'$epp/ep$',
        PMISS_LABELS, ylim=(0, 0.25), label_xy=(2.3, 0.2), ylabel_size=15,
        save_as=SAVE_AS.get('q2_2x2_ratio'))

    if has_sim:
        ph.plot_q2_2x2_data_over_sim_ratio_only(
            pdf, data, sims, 'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', PMISS_LABELS,
            ylim=(0.5, 1.5), label_xy=(2.3, 1.35), ylabel_size=13,
            save_as=SAVE_AS.get('q2_2x2_data_over_sim_ratio_only'))

    # ---------------------------------------------------------------------
    # Fit-derived quantities (from ExtractFitQuantities.C)
    # ---------------------------------------------------------------------
    fit_q2_single = [
        ('sigma_pcmx', 'epp', ()),
        ('sigma_pcmy', 'epp', ()),
    ]
    ok, miss = _require_graphs([data] + sims, fit_q2_single)
    if ok:
        ph.plot_q2_single(
            pdf, data, sims, 'sigma_pcmx', r'$\sigma_{x,C.M.} [GeV]$',
            selection='epp', xlim=(1.5, 5.0),ylim=(0,.35), with_ratio=has_sim,
            unit_scale=True,
            suptitle=r'$\sigma_{x,C.M.} [GeV]$',
            save_as=SAVE_AS.get('sigma_pcmx'))
        ph.plot_q2_single(
            pdf, data, sims, 'sigma_pcmy', r'$\sigma_{y,C.M.} [GeV]$',
            selection='epp', xlim=(1.5, 5.0),ylim=(0,.35), with_ratio=has_sim,
            unit_scale=True,
            suptitle=r'$\sigma_{y,C.M.} [GeV]$',
            save_as=SAVE_AS.get('sigma_pcmy'))
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
        ok, miss = _require_graphs([data] + sims, needed)
        if not ok:
            missing_notes.append(
                'Skipped panel for %s because required extracted graphs were missing.' % task_name)
            missing_notes.extend(miss)
            continue
        ph.plot_q2_2x2(
            pdf, data, sims, task_name, ylabel,
            PMISS_LABELS, selection=selection,
            ylim=ylim, label_xy=label_xy, ylabel_size=ylabel_size,
            big_label=big_label, unit_scale=True,
            suptitle='Extracted quantity: %s (overlay)' % task_name,
            save_as=SAVE_AS.get('q2_2x2_%s' % task_name))

        if has_sim:
            ph.plot_q2_2x2_data_over_sim(
                pdf, data, sims, task_name,
                PMISS_LABELS, selection=selection,
                ylim=(0.5, 1.5), label_xy=(2.3, 1.35), ylabel_size=13,
                unit_scale=True,
                suptitle=ylabel,
                save_as=SAVE_AS.get('q2_2x2_data_over_sim_%s' % task_name))

    pdf.close()
    if missing_notes:
        print('Some requested helium plots were skipped due to missing graphs:')
        for m in missing_notes:
            print(' - %s' % m)
    print('Wrote %s' % args.out)


if __name__ == '__main__':
    main()
