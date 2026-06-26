"""
plot_Emiss_pmiss.py

E_miss distributions in pMiss bins for the analysis note: E1miss for e'p
and e'pp side by side, and E2miss for e'pp only (no e'p counterpart, since
E2miss requires a detected recoil). Each row is one pMiss bin
(bE_pmiss = [0.4, 0.55, 0.7, 0.85, 1.0] in Main_Figs_Binned.cpp, labeled via
labels.PMISS_LABELS), with a single boxed pMiss-range label centered above
the row and a thin vertical line at the quasi-free single-nucleon-knockout
threshold E_miss = sqrt(pMiss^2 + mN^2) - mN, evaluated at the row's mean
pMiss (PMISS_CENTERS_EP/EPP below -- measured directly from per-event
kinematics by Ana/Q2_Ana/AvgPMiss.cpp, run once on data; not derived from
the pMiss_ep_note/pMiss_epp_note histograms here, since estimating the
mean from those 50-bin histograms turned out to be distorted by bin-grid
misalignment between the 50 uniform note-bins and these 4 wider edges).
Written out as individual single-figure PDFs under pdf/Emiss_pmiss/.

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

# Per-event mean pMiss in each of Main_Figs_Binned.cpp's bE_pmiss bins
# ([0.4, 0.55, 0.7, 0.85, 1.0]), measured directly from data by running
# Ana/Q2_Ana/AvgPMiss.cpp (not derived here from histograms -- see the
# module docstring for why). Re-run AvgPMiss.cpp and update these if the
# data sample, cuts, or bE_pmiss change.
PMISS_CENTERS_EP = [0.499816, 0.63018, 0.769554, 0.91613]
PMISS_CENTERS_EPP = [0.508686, 0.634642, 0.772566, 0.915738]
PMISS_BINS = [0.4, 0.55, 0.7, 0.85, 1.0]

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

    pdf = _NoOpPdf()

    ph.plot_emiss_4x2_note(
        pdf, series, r'$E_{1,miss} [GeV]$',
        'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
        pmiss_centers_ep=PMISS_CENTERS_EP, pmiss_centers_epp=PMISS_CENTERS_EPP,
        box_labels=PMISS_LABELS, mN=MN, draw_threshold_lines=False,
        save_as=os.path.join(args.out_dir, 'E1miss_4x2.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E1miss_4x2.pdf'))

    ph.plot_emiss_4x1_note(
        pdf, series, r'$E_{0,miss} [GeV]$',
        'E0miss_ep_SRC_pmiss',
        pmiss_centers=PMISS_CENTERS_EP, box_labels=PMISS_LABELS, mN=MN,
        xlim=(0.1, 0.5), selection='ep', panel_title=r'$(e,e^{\prime}p)$',
        save_as=os.path.join(args.out_dir, 'E0miss_4x1_ep.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E0miss_4x1_ep.pdf'))

    ph.plot_emiss_4x1_note(
        pdf, series, r'$E_{2,miss} [GeV]$',
        'E2miss_epp_SRC_pmiss',
        pmiss_centers=PMISS_CENTERS_EPP, box_labels=PMISS_LABELS, mN=MN,
        xlim=(-0.2, 0.4), draw_threshold_lines=False,
        save_as=os.path.join(args.out_dir, 'E2miss_4x1.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E2miss_4x1.pdf'))

    ph.plot_emiss_waterfall(
        pdf, PMISS_BINS, data, sims,
        task_name='E0miss_ep_SRC_pmiss', selection='ep',
        xlim=(0.1, 0.5), ylim=(0.0, 120.0),
        xlabel=r'$E_{0,miss} [GeV]$',
        panel_label=r'$(e,e^{\prime}p)$',
        save_as=os.path.join(args.out_dir, 'E0miss_waterfall_ep.pdf'))
    print('Wrote %s' % os.path.join(args.out_dir, 'E0miss_waterfall_ep.pdf'))


if __name__ == '__main__':
    main()
