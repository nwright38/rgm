"""
plot_Q2_pmiss_extracted.py

Extracted-quantity plotting driver focused on Q2 and pMiss-binned Q2 panels.
Outputs individual PDFs into pdf/Q2_pmiss/ (or --out-dir).

This mirrors the fit-derived and ratio pages from plot_data_vs_sim.py, but
writes each requested figure to its own file for note-ready usage.

Usage:
    python plot_Q2_pmiss_extracted.py Data_He.root
    python plot_Q2_pmiss_extracted.py Data_He.root Sim_He_AV18.root --sim-label AV18
"""

import argparse
import os
import subprocess

import graph_io
import normalization
from labels import PMISS_LABELS
from plot_helpers import Series
import plot_helpers as ph

DEFAULT_SIM_COLORS = ['red', 'blue', 'green', 'darkorange', 'purple', 'brown']


class _NoOpPdf(object):
    """Placeholder PdfPages-like object when using save_as-only outputs."""

    def savefig(self, fig):
        pass


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('data_file', help='Main_Figs_Binned+BuildGraphs output for data')
    p.add_argument('sim_files', nargs='*',
                   help='Zero or more Main_Figs_Binned+BuildGraphs sim outputs')
    p.add_argument('--data-label', default='Data')
    p.add_argument('--data-color', default='black')
    p.add_argument('--sim-label', action='append', default=[])
    p.add_argument('--sim-color', action='append', default=[])
    p.add_argument('--out-dir', default='pdf/Q2_pmiss')
    p.add_argument('--extract-fit', action='store_true',
                   help='Run ExtractFitQuantities.C on data and sim files before plotting')
    return p.parse_args()


def _build_sims(args, f_data):
    sims = []
    for i, sim_file in enumerate(args.sim_files):
        label = args.sim_label[i] if i < len(args.sim_label) else 'Sim %d' % (i + 1)
        color = args.sim_color[i] if i < len(args.sim_color) else \
            DEFAULT_SIM_COLORS[i % len(DEFAULT_SIM_COLORS)]
        f_sim = graph_io.open_file(sim_file)

        scale_ep, _ = normalization.scale_factor(f_sim, f_data, 'ep')
        scale_epp = scale_ep
        sims.append(Series(label, f_sim, color, 'sim', scale_ep, scale_epp))
    return sims


def _require_graphs(series_list, needed):
    """Return (ok, missing_messages) for required task/selection/axis_bin entries."""
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


def _save_call(callable_fn, out_path, missing_notes):
    try:
        callable_fn(out_path)
        print('Wrote %s' % out_path)
    except Exception as exc:
        missing_notes.append('Skipped %s: %s' % (out_path, exc))


def _add_missing_summary(missing_notes, header, miss, max_items=6):
    missing_notes.append(header)
    if not miss:
        return
    missing_notes.extend(miss[:max_items])
    if len(miss) > max_items:
        missing_notes.append('... and %d more missing-graph entries.' % (len(miss) - max_items))


def _run_extract_fit(filename):
    macro = 'ExtractFitQuantities.C("%s")' % filename
    cmd = ['root', '-l', '-b', '-q', macro]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError('ExtractFitQuantities failed for %s:\n%s'
                           % (filename, proc.stdout.decode('utf-8', errors='replace')))


def _plot_sigma_q2(pdf, data, sims, out_dir, has_sim, missing_notes):
    sigma_specs = [
        ('sigma_pcmx', r'$\sigma_{x,C.M.} [GeV]$', 'sigma_cmx_vs_q2.pdf'),
        ('sigma_pcmy', r'$\sigma_{y,C.M.} [GeV]$', 'sigma_cmy_vs_q2.pdf'),
    ]
    for task_name, ylabel, filename in sigma_specs:
        needed = [(task_name, 'epp', ())]
        ok, miss = _require_graphs([data] + sims, needed)
        if not ok:
            _add_missing_summary(
                missing_notes,
                'Skipped %s (missing required graphs).' % filename,
                miss)
            continue

        out_path = os.path.join(out_dir, filename)
        _save_call(
            lambda path: ph.plot_q2_single(
                pdf, data, sims, task_name, ylabel,
                selection='epp', xlim=(1.5, 5.0), ylim=(0.0, 0.35),
                with_ratio=has_sim, unit_scale=True, save_as=path),
            out_path, missing_notes)


def _plot_sigma_q2_pmiss_if_available(pdf, data, sims, out_dir, has_sim, missing_notes):
    # Optional support: if a file contains sigma-vs-Q2 graphs split by pMiss,
    # draw those 2x2 panels. Current standard extraction may not provide these.
    sigma_panel_specs = [
        ('sigma_pcmx_pmiss', r'$\sigma_{x,C.M.} [GeV]$', 'sigma_cmx_vs_q2_pmiss_2x2.pdf'),
        ('sigma_pcmy_pmiss', r'$\sigma_{y,C.M.} [GeV]$', 'sigma_cmy_vs_q2_pmiss_2x2.pdf'),
    ]
    for task_name, ylabel, filename in sigma_panel_specs:
        needed = [
            (task_name, 'epp', (0,)),
            (task_name, 'epp', (1,)),
            (task_name, 'epp', (2,)),
            (task_name, 'epp', (3,)),
        ]
        ok, _ = _require_graphs([data] + sims, needed)
        if not ok:
            missing_notes.append(
                'Skipped %s because task "%s" in pMiss bins was not found.'
                % (filename, task_name))
            continue

        out_path = os.path.join(out_dir, filename)
        _save_call(
            lambda path: ph.plot_q2_2x2(
                pdf, data, sims, task_name, ylabel, PMISS_LABELS,
                selection='epp', ylim=(0.0, 0.35),
                ylabel_size=13, unit_scale=True, save_as=path),
            out_path, missing_notes)

        if has_sim:
            ratio_out = os.path.join(
                out_dir, filename.replace('.pdf', '_data_over_sim.pdf'))
            _save_call(
                lambda path: ph.plot_q2_2x2_data_over_sim(
                    pdf, data, sims, task_name, PMISS_LABELS,
                    selection='epp', ylim=(0.5, 1.5),
                    ylabel_size=13, unit_scale=True, save_as=path),
                ratio_out, missing_notes)


def _plot_ratios(pdf, data, sims, out_dir, has_sim, missing_notes):
    _save_call(
        lambda path: ph.plot_overlay(
            pdf, 'pMiss_epp_over_pMiss_ep', [data] + sims, selection='ratio',
            xlabel=r'$p_{miss} [GeV]$', ylabel=r'$epp/ep$',
            xlim=(0.4, 1.0), ylim=(0.0, 0.25), with_ratio=has_sim,
            xlabel_size=25, ylabel_size=25, save_as=path),
        os.path.join(out_dir, 'epp_over_ep_vs_pmiss.pdf'), missing_notes)

    q2_ratio_needed = [
        ('Q2_epp_SRC_pmiss_over_Q2_ep_SRC_pmiss', 'ratio', (0,)),
        ('Q2_epp_SRC_pmiss_over_Q2_ep_SRC_pmiss', 'ratio', (1,)),
        ('Q2_epp_SRC_pmiss_over_Q2_ep_SRC_pmiss', 'ratio', (2,)),
        ('Q2_epp_SRC_pmiss_over_Q2_ep_SRC_pmiss', 'ratio', (3,)),
    ]
    ok, miss = _require_graphs([data] + sims, q2_ratio_needed)
    if not ok:
        _add_missing_summary(
            missing_notes,
            'Skipped epp/ep vs Q2 in pMiss bins (missing required ratio graphs).',
            miss)
        return

    _save_call(
        lambda path: ph.plot_q2_2x2_ratio(
            pdf, data, sims,
            'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', r'$epp/ep$', PMISS_LABELS,
            ylim=(0.0, 0.25), label_xy=(2.3, 0.2), ylabel_size=15,
            save_as=path),
        os.path.join(out_dir, 'epp_over_ep_vs_q2_pmiss_2x2.pdf'), missing_notes)

    if has_sim:
        _save_call(
            lambda path: ph.plot_q2_2x2_data_over_sim_ratio_only(
                pdf, data, sims,
                'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', PMISS_LABELS,
                ylim=(0.5, 1.5), label_xy=(2.3, 1.35), ylabel_size=13,
                save_as=path),
            os.path.join(out_dir, 'epp_over_ep_vs_q2_pmiss_2x2_data_over_sim.pdf'),
            missing_notes)


def _plot_emiss_stats_q2_pmiss(pdf, data, sims, out_dir, has_sim, missing_notes):
    panel_specs = [
        ('stddev_E1miss_ep_pmiss', 'ep',
         r'Std Dev $E_{1,miss} [GeV]$', (0.0, 0.15),
         'stddev_E1miss_ep_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}p)$'),
        ('stddev_E1miss_epp_pmiss', 'epp',
         r'Std Dev $E_{1,miss} [GeV]$', (0.0, 0.15),
         'stddev_E1miss_epp_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}pp)$'),
        ('stddev_E2miss_epp_pmiss', 'epp',
         r'Std Dev $E_{2,miss} [GeV]$', (0.0, 0.15),
         'stddev_E2miss_epp_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}pp)$'),
        ('mean_E1miss_ep_pmiss', 'ep',
         r'Mean $E_{1,miss} [GeV]$', (-0.02, 0.3),
         'mean_E1miss_ep_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}p)$'),
        ('mean_E1miss_epp_pmiss', 'epp',
         r'Mean $E_{1,miss} [GeV]$', (-0.02, 0.3),
         'mean_E1miss_epp_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}pp)$'),
        ('mean_E2miss_epp_pmiss', 'epp',
         r'Mean $E_{2,miss} [GeV]$', (-0.02, 0.3),
         'mean_E2miss_epp_vs_q2_pmiss_2x2.pdf', r'$(e,e^{\prime}pp)$'),
    ]

    for task_name, selection, ylabel, ylim, filename, big_label in panel_specs:
        needed = [
            (task_name, selection, (0,)),
            (task_name, selection, (1,)),
            (task_name, selection, (2,)),
            (task_name, selection, (3,)),
        ]
        ok, miss = _require_graphs([data] + sims, needed)
        if not ok:
            _add_missing_summary(
                missing_notes,
                'Skipped %s (missing required graphs).' % filename,
                miss)
            continue

        out_path = os.path.join(out_dir, filename)
        _save_call(
            lambda path: ph.plot_q2_2x2(
                pdf, data, sims, task_name, ylabel, PMISS_LABELS,
                selection=selection, ylim=ylim, ylabel_size=10,
                big_label=big_label, unit_scale=True, save_as=path),
            out_path, missing_notes)

        if has_sim:
            ratio_out = os.path.join(
                out_dir, filename.replace('.pdf', '_data_over_sim.pdf'))
            _save_call(
                lambda path: ph.plot_q2_2x2_data_over_sim(
                    pdf, data, sims, task_name, PMISS_LABELS,
                    selection=selection, ylim=(0.5, 1.5),
                    ylabel_size=13, unit_scale=True, save_as=path),
                ratio_out, missing_notes)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.extract_fit:
        _run_extract_fit(args.data_file)
        for sim_file in args.sim_files:
            _run_extract_fit(sim_file)

    f_data = graph_io.open_file(args.data_file)
    data = Series(args.data_label, f_data, args.data_color, 'data')
    sims = _build_sims(args, f_data)
    has_sim = len(sims) > 0

    pdf = _NoOpPdf()
    missing_notes = []

    _plot_sigma_q2(pdf, data, sims, args.out_dir, has_sim, missing_notes)
    _plot_sigma_q2_pmiss_if_available(pdf, data, sims, args.out_dir, has_sim, missing_notes)
    _plot_ratios(pdf, data, sims, args.out_dir, has_sim, missing_notes)
    _plot_emiss_stats_q2_pmiss(pdf, data, sims, args.out_dir, has_sim, missing_notes)

    if missing_notes:
        print('Some requested Q2/pMiss extracted-quantity plots were skipped:')
        for m in missing_notes:
            print(' - %s' % m)


if __name__ == '__main__':
    main()
