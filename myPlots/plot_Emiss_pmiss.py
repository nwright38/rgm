"""
plot_Emiss_pmiss.py

E_miss distributions in pMiss bins for the analysis note: E1miss for e'p
and e'pp side by side, and E2miss for e'pp only (no e'p counterpart, since
E2miss requires a detected recoil). Each row is one pMiss bin
(bE_pmiss = [0.4, 0.55, 0.7, 0.85, 1.0] in Main_Figs_Binned.cpp, labeled via
labels.PMISS_LABELS), with a single boxed pMiss-range label centered above
the row and a thin vertical line at the quasi-free single-nucleon-knockout
threshold E_miss = sqrt(pMiss^2 + mN^2) - mN, evaluated at the row's mean
pMiss (computed separately for the e'p and e'pp populations from the
finer-binned pMiss_ep_note/pMiss_epp_note histograms, not the bin-edge
midpoint). Written out as individual single-figure PDFs under
pdf/Emiss_pmiss/.

Usage:
    python plot_Emiss_pmiss.py Data_He.root
    python plot_Emiss_pmiss.py Data_He.root Sim_He_AV18.root --sim-label AV18
"""

import argparse
import os

import graph_io
import normalization
from plot_helpers import Series
import plot_helpers as ph
from labels import PMISS_LABELS

MN = 0.938272  # matches SRC_Cuts.cpp's mN
PMISS_EDGES = [0.4, 0.55, 0.7, 0.85, 1.0]

DEFAULT_SIM_COLORS = ['red', 'blue', 'green', 'darkorange']


def _mean_pmiss_per_bin(data, note_task_name, selection, edges):
    """Mean pMiss in each [edges[i], edges[i+1]) range, weighted by counts
    in the finer-binned pMiss_ep_note/pMiss_epp_note histograms (50 uniform
    bins over 0.4-1.2 GeV) -- used instead of the naive bin-edge midpoint
    so the E_miss threshold line reflects where the data actually sits
    within each wide (and, for the last bin, asymmetric) pMiss bin."""
    x, y, _ = ph.get_xy_err(data, note_task_name, selection)
    centers = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        num = sum(xi * yi for xi, yi in zip(x, y) if lo <= xi < hi)
        den = sum(yi for xi, yi in zip(x, y) if lo <= xi < hi)
        centers.append(num / den if den > 0 else 0.5 * (lo + hi))
    return centers


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
    p.add_argument('--out-dir', default='pdf/Emiss_pmiss')
    return p.parse_args()


def _build_sims(args, f_data):
    if len(args.sim_files) > 4:
        raise ValueError('plot_Emiss_pmiss.py supports at most 4 simulation files, got %d'
                         % len(args.sim_files))
    sims = []
    for i, sim_file in enumerate(args.sim_files):
        label = args.sim_label[i] if i < len(args.sim_label) else 'Sim %d' % (i + 1)
        color = args.sim_color[i] if i < len(args.sim_color) else \
            DEFAULT_SIM_COLORS[i % len(DEFAULT_SIM_COLORS)]
        f_sim = graph_io.open_file(sim_file)
        scale_ep, _ = normalization.scale_factor(f_sim, f_data, 'ep')
        # Same epp-scaled-by-ep-factor convention as plot_1D.py's _build_sims.
        scale_epp = scale_ep
        sims.append(Series(label, f_sim, color, 'sim', scale_ep, scale_epp))
    return sims


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    f_data = graph_io.open_file(args.data_file)
    data = Series(args.data_label, f_data, args.data_color, 'data')
    sims = _build_sims(args, f_data)
    series = [data] + sims

    pmiss_centers_ep = _mean_pmiss_per_bin(data, 'pMiss_ep_note', 'ep', PMISS_EDGES)
    pmiss_centers_epp = _mean_pmiss_per_bin(data, 'pMiss_epp_note', 'epp', PMISS_EDGES)

    pdf = _NoOpPdf()

    ph.plot_emiss_4x2_note(
        pdf, series, r'$E_{1,miss} [GeV]$',
        'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
        pmiss_centers_ep=pmiss_centers_ep, pmiss_centers_epp=pmiss_centers_epp,
        box_labels=PMISS_LABELS, mN=MN,
        save_as=os.path.join(args.out_dir, 'E1miss_4x2.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E1miss_4x2.pdf'))

    ph.plot_emiss_4x1_note(
        pdf, series, r'$E_{2,miss} [GeV]$',
        'E2miss_epp_SRC_pmiss',
        pmiss_centers=pmiss_centers_epp, box_labels=PMISS_LABELS, mN=MN,
        xlim=(-0.2, 0.4),
        save_as=os.path.join(args.out_dir, 'E2miss_4x1.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E2miss_4x1.pdf'))


if __name__ == '__main__':
    main()
