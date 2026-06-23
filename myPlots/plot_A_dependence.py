"""
plot_A_dependence.py

Cross-nucleus comparisons (He + C + Ca40 + Ca48 + Sn), kept in its own
script and separate from plot_data_vs_sim.py since these are about
A-dependence, not data/simulation agreement -- there is no single
simulation reference here, so with_ratio is off throughout.

Each nucleus's data is scaled to the He reference using the corrected
normalization (normalization.py, reading the integrated Q2 yield), not the
buggy bin-scanning getFactor_ep/getFactor_epp this replaces.
"""

import matplotlib.backends.backend_pdf

import graph_io
import nuclei_config as cfg
import normalization
from labels import PMISS_LABELS
from plot_helpers import Series
import plot_helpers as ph

pdf = matplotlib.backends.backend_pdf.PdfPages("plots/A_dependence.pdf")

he_entry = cfg.reference_nucleus()
f_he_data = graph_io.open_file(he_entry['file'])
he_data = Series(he_entry['label'], f_he_data, he_entry['color'], 'data')

nuclei = []
for n in cfg.NUCLEI:
    if n['key'] == cfg.REFERENCE_KEY:
        continue
    f_nuc = graph_io.open_file(n['file'])
    scale_ep, _ = normalization.scale_factor(f_nuc, f_he_data, 'ep')
    scale_epp, _ = normalization.scale_factor(f_nuc, f_he_data, 'epp')
    nuclei.append(Series(n['label'], f_nuc, n['color'], 'data_ex',
                         scale_ep, scale_epp, n['offset'], n['marker']))

# 48Ca omitted from the C.M. figures, matching the original.
nuclei_no48 = [n for n in nuclei if n.label != cfg.find(cfg.NUCLEI, 'Ca48')['label']]

annotations = [{'x': 0.5, 'y': 0.1 + 0.02 * i, 'text': s.label, 'color': s.color}
              for i, s in enumerate([he_data] + nuclei)]

# ---------------------------------------------------------------------------
# C.M. projections
# ---------------------------------------------------------------------------
for axis in ('x', 'y'):
    ph.plot_overlay(
        pdf, 'pcm{0}_epp'.format(axis), [he_data] + nuclei_no48, selection='epp',
        xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
        ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 400), with_ratio=False)

# ---------------------------------------------------------------------------
# epp/ep vs pMiss, all nuclei (already a ratio quantity -- no rescaling)
# ---------------------------------------------------------------------------
ph.plot_overlay(
    pdf, 'pMiss_epp_over_pMiss_ep', [he_data] + nuclei, selection='ratio',
    xlabel=r'$p_{miss}$', ylabel=r'$epp/ep$',
    xlim=(0.4, 1.05), ylim=(0, 0.2), with_ratio=False,
    xlabel_size=25, ylabel_size=25, annotations=annotations)

# ---------------------------------------------------------------------------
# E1miss 4x2, E2miss 4x1, all nuclei (narrow offsets via offset_scale=0.4)
# ---------------------------------------------------------------------------
ph.plot_emiss_4x2(
    pdf, [he_data] + nuclei, r'$E_{1,miss} [GeV]$',
    'E1miss_ep_SRC_pmiss', 'E1miss_epp_SRC_pmiss',
    ylim_ep=[1000, 1500, 1500, 1000], ylim_epp=[60, 110, 150, 100],
    labely_ep=[850, 1250, 1250, 850], labely_epp=[50, 95, 130, 85],
    pmiss_labels=PMISS_LABELS, offset_scale=0.4)

ph.plot_emiss_4x1(
    pdf, [he_data] + nuclei, r'$E_{2,miss} [GeV]$',
    'E2miss_epp_SRC_pmiss',
    ylim=[30, 110, 130, 100], labely=[25, 95, 110, 85],
    pmiss_labels=PMISS_LABELS, offset_scale=0.4)

pdf.close()
