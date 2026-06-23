"""
plot_data_vs_sim.py

He data vs the five theory models, every variable, each overlay with a
data/sim ratio subpanel underneath (plot_overlay's with_ratio=True default).

This is the data-vs-simulation half of what used to be one monolithic
makePlots.py; plot_A_dependence.py is the cross-nucleus half. Both read
from nuclei_config.py instead of hardcoding file paths/Series chains, and
both expect BuildGraphs.cpp to have already been run on every input file
(it adds the "graphs" TDirectory these scripts read from).

NOTE: the old C.M.-width-vs-Q^2 and E_miss-mean/width-vs-Q^2 plots
(g_sigma_pcmx, g_sigma_E1miss_ep_pmiss, g_mean_E1miss_ep_pmiss, ...) are not
reproduced here -- they came from Gaussian fits across the 100 toy
histograms, which is a separate fit-extraction step not yet built on top of
the persisted hists/nominal + hists/toy_NNN histograms. See the note in
plot_helpers.py.
"""

import matplotlib.backends.backend_pdf

import graph_io
import nuclei_config as cfg
import normalization
from labels import PMISS_LABELS
from plot_helpers import Series
import plot_helpers as ph

pdf = matplotlib.backends.backend_pdf.PdfPages("plots/data_vs_sim.pdf")

# ---------------------------------------------------------------------------
# Open files, build Series
# ---------------------------------------------------------------------------
he_entry = cfg.reference_nucleus()
f_he_data = graph_io.open_file(he_entry['file'])
he_data = Series(he_entry['label'], f_he_data, he_entry['color'], 'data')

sim_models = []
for m in cfg.SIM_MODELS:
    f_sim = graph_io.open_file(m['file'])
    scale_ep, _ = normalization.scale_factor(f_sim, f_he_data, 'ep')
    scale_epp, _ = normalization.scale_factor(f_sim, f_he_data, 'epp')
    sim_models.append(Series(m['label'], f_sim, m['color'], 'sim', scale_ep, scale_epp))

# AV4' appears only in the pMiss ratio overlay; the E_miss panels use the
# subset without it (matching the original figures).
sim_models_emiss = [s for s in sim_models if s.label != "AV4'"]

# ---------------------------------------------------------------------------
# C.M. projections
# ---------------------------------------------------------------------------
for axis in ('x', 'y'):
    ph.plot_overlay(
        pdf, 'pcm{0}_epp'.format(axis), [he_data] + sim_models, selection='epp',
        xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
        ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 350))

# ---------------------------------------------------------------------------
# epp/ep vs pMiss, all five models (already a ratio quantity -- no rescaling)
# ---------------------------------------------------------------------------
ph.plot_overlay(
    pdf, 'pMiss_epp_over_pMiss_ep', [he_data] + sim_models, selection='ratio',
    xlabel=r'$p_{miss} [GeV]$', ylabel=r'$epp/ep$',
    xlim=(0.4, 1.0), ylim=(0, 0.2), with_ratio=False,
    xlabel_size=25, ylabel_size=25)

# ---------------------------------------------------------------------------
# E1miss 4x2 (ep | epp), E2miss 4x1 (epp)
# ---------------------------------------------------------------------------
ph.plot_emiss_4x2(
    pdf, [he_data] + sim_models_emiss, r'$E_{1,miss} [GeV]$',
    'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
    ylim_ep=[1000, 1000, 850, 600], ylim_epp=[60, 110, 110, 60],
    labely_ep=[850, 850, 700, 500], labely_epp=[50, 95, 95, 50],
    pmiss_labels=PMISS_LABELS)

ph.plot_emiss_4x1(
    pdf, [he_data] + sim_models_emiss, r'$E_{2,miss} [GeV]$',
    'E2miss_epp_SRC_pmiss',
    ylim=[60, 110, 110, 80], labely=[50, 95, 95, 65],
    pmiss_labels=PMISS_LABELS)

# ---------------------------------------------------------------------------
# epp/ep vs Q^2 in pMiss bins (one model shown, matching the original)
# ---------------------------------------------------------------------------
gcf = sim_models[0]  # AV18, reused as the single comparison curve, as before
ph.plot_q2_2x2_ratio(
    pdf, he_data, gcf, 'Q2_epp_SRC_pmiss', 'Q2_ep_SRC_pmiss', r'$epp/ep$',
    PMISS_LABELS, ylim=(0, 0.25), label_xy=(2.3, 0.2), ylabel_size=15)

pdf.close()
