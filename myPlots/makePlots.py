from __future__ import division
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import ROOT
import plotTOOL as pto
import plot_helpers as ph
from plot_helpers import Series

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
Q2bins = [1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0]
pMbins = [0.4, 0.55, 0.7, 0.85, 1.0]
kMbins = [0.3, 0.45, 0.6, 0.75, 0.9]

PMISS_LABELS = [
    r'$0.4<p_{Miss}<0.55$',
    r'$0.55<p_{Miss}<0.7$',
    r'$0.7<p_{Miss}<0.85$',
    r'$0.85<p_{Miss}<1.0$',
]

path = '/work/clas12/users/awild/RGM/CLAS_note/Boston2024/Main/Dec2024/'
pdf = matplotlib.backends.backend_pdf.PdfPages("plots/output.pdf")

# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------
def _open(name):
    return ROOT.TFile.Open(path + 'Sys_Err/' + name)

f_Data_He = ROOT.TFile.Open('/work/clas12/users/nwright/RGM/Q2/Sys_Err/rootFiles/Data_He.root')
#f_Data_He = _open('Data_C.root')
f_Data_C = _open('Data_C.root')
f_Data_40Ca = _open('Data_40Ca.root')
f_Data_48Ca = _open('Data_48Ca.root')
f_Data_Sn = _open('Data_Sn.root')

f_Sim_He_AV18 =  ROOT.TFile.Open('Sim_He_AV18_testSkim.root')
#f_Sim_He_AV18 =  ROOT.TFile.Open('/work/clas12/users/nwright/rgm_andrew/build/Ana/Q2_Ana/c12_reweighted_sim.root')
f_Sim_He_AV4 = _open('Sim_He_AV4.root')
f_Sim_He_N2LO10 = _open('Sim_He_N2LO10.root')
f_Sim_He_N2LO12 = _open('Sim_He_N2LO12.root')
f_Sim_He_NV = _open('Sim_He_NV.root')

# ---------------------------------------------------------------------------
# Normalisation factors (epp / ep), all relative to He data
# ---------------------------------------------------------------------------
C2epp = pto.getFactor_epp(f_Data_C, f_Data_He)
C2ep = pto.getFactor_ep(f_Data_C, f_Data_He)
C0epp = pto.getFactor_epp(f_Data_40Ca, f_Data_He)
C0ep = pto.getFactor_ep(f_Data_40Ca, f_Data_He)
C8epp = pto.getFactor_epp(f_Data_48Ca, f_Data_He)
C8ep = pto.getFactor_ep(f_Data_48Ca, f_Data_He)
Sepp = pto.getFactor_epp(f_Data_Sn, f_Data_He)
Sep = pto.getFactor_ep(f_Data_Sn, f_Data_He)

A8epp = pto.getFactor_epp(f_Sim_He_AV18, f_Data_He)
A8ep = pto.getFactor_ep(f_Sim_He_AV18, f_Data_He)
A4epp = pto.getFactor_epp(f_Sim_He_AV4, f_Data_He)
A4ep = pto.getFactor_ep(f_Sim_He_AV4, f_Data_He)
N0epp = pto.getFactor_epp(f_Sim_He_N2LO10, f_Data_He)
N0ep = pto.getFactor_ep(f_Sim_He_N2LO10, f_Data_He)
N2epp = pto.getFactor_epp(f_Sim_He_N2LO12, f_Data_He)
N2ep = pto.getFactor_ep(f_Sim_He_N2LO12, f_Data_He)
NVepp = pto.getFactor_epp(f_Sim_He_NV, f_Data_He)
NVep = pto.getFactor_ep(f_Sim_He_NV, f_Data_He)

# Per-Q^2-bin TGraph factors for the C.M. projections
N0epp_array = []
for i in range(0, 7):
    N0epp_array.append(pto.getFactorTGraph(
        f_Sim_He_N2LO10, 'h_pcmx_epp_SRC_Q2_' + str(i),
        f_Data_He, 'h_pcmx_epp_SRC_Q2_' + str(i)))

# ---------------------------------------------------------------------------
# Reusable Series definitions
# ---------------------------------------------------------------------------
he_data = Series(r'${}^{4}He$', f_Data_He, 'black', 'data')

# He simulation models. AV4' appears only in the pMiss ratio overlay, so the
# E_miss panels use the subset without it (matching the original figures).
sim_models = [
    Series("AV18",         f_Sim_He_AV18,   'red',    'sim', A8ep, A8epp),
    Series("AV4'",         f_Sim_He_AV4,    'pink',   'sim', A4ep, A4epp),
    Series("N2LO (1.0fm)", f_Sim_He_N2LO10, 'blue',   'sim', N0ep, N0epp),
    Series("N2LO (1.2fm)", f_Sim_He_N2LO12, 'green',  'sim', N2ep, N2epp),
    Series("Norfolk",      f_Sim_He_NV,     'orange', 'sim', NVep, NVepp),
]
sim_models_emiss = [s for s in sim_models if s.label != "AV4'"]

# The GCF (N2LO 1.0 fm) curve, reused as the single comparison line.
gcf = Series('GCF', f_Sim_He_AV18, 'blue', 'sim')
#gcf = Series('GCF', f_Sim_He_N2LO10, 'blue', 'sim')

# Other nuclei, scaled to He. Base offsets are the "wide" values used in the
# C.M./pMiss figures; the E_miss panels pass offset_scale=0.4 to get the
# original narrow offsets (-0.004, -0.002, 0.002, 0.004).
nuc_C = Series(r'${}^{12}C$',   f_Data_C,    'darkred',       'data_ex', C2ep, C2epp, -0.01,  '<')
nuc_Ca40 = Series(r'${}^{40}Ca$',  f_Data_40Ca, 'darkgreen',     'data_ex', C0ep, C0epp, -0.005, 'v')
nuc_Ca48 = Series(r'${}^{48}Ca$',  f_Data_48Ca, 'darkblue',      'data_ex', C8ep, C8epp,  0.005, '^')
nuc_Sn = Series(r'${}^{120}Sn$', f_Data_Sn,   'darkgoldenrod', 'data_ex', Sep,  Sepp,   0.01,  '>')
nuclei = [nuc_C, nuc_Ca40, nuc_Ca48, nuc_Sn]
nuclei_no48 = [nuc_C, nuc_Ca40, nuc_Sn]  # 48Ca omitted in the C.M. figures

# ===========================================================================
# He: data vs simulations
# ===========================================================================

# C.M. projections with the GCF curve (per-bin TGraph factor, Q^2 bin 3)
gcf_pcm = Series('GCF', f_Sim_He_N2LO10, 'blue', 'sim', scale_epp=N0epp_array[3])
for axis in ('x', 'y'):
    ph.plot_overlay(
        pdf, 'h_pcm{0}_epp'.format(axis), [he_data, gcf_pcm],
        xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
        ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 350),
        annotations=[
            {'x': -0.3, 'y': 300, 'text': 'Data', 'color': 'black', 'fontsize': 25},
            {'x': -0.3, 'y': 250, 'text': 'GCF',  'color': 'blue',  'fontsize': 25},
        ])

# epp/ep vs pMiss, all five models (ratio histogram -> unit_scale)
ph.plot_overlay(
    pdf, 'pMiss_epp', [he_data] + sim_models,
    xlabel=r'$p_{miss} [GeV]$', ylabel=r'$epp/ep$',
    xlim=(0.4, 1.0), ylim=(0, 0.2), unit_scale=True,
    xlabel_size=25, ylabel_size=25,
    annotations=[
        {'x': 0.65, 'y': 0.13,  'text': 'AV18',         'color': 'red',   'fontsize': 15},
        {'x': 0.45, 'y': 0.175, 'text': "AV4'",         'color': 'black', 'fontsize': 15,
         'bbox': dict(facecolor='pink', edgecolor='none')},
        {'x': 0.7,  'y': 0.07,  'text': 'N2LO (1.0fm)', 'color': 'blue',  'fontsize': 15},
        {'x': 0.92, 'y': 0.1,   'text': 'Norfolk',      'color': 'black', 'fontsize': 15,
         'bbox': dict(facecolor='orange', edgecolor='none')},
        {'x': 0.8,  'y': 0.05,  'text': 'N2LO (1.2fm)', 'color': 'green', 'fontsize': 15},
    ])

# E1miss 4x2 (ep | epp), data + models
ph.plot_emiss_4x2(
    pdf, [he_data] + sim_models_emiss, r'$E_{1,miss} [GeV]$',
    'h_E1miss_ep_SRC_pmiss', 'h_E1miss_epp_SRC_pmiss',
    ylim_ep=[1000, 1000, 850, 600], ylim_epp=[60, 110, 110, 60],
    labely_ep=[850, 850, 700, 500], labely_epp=[50, 95, 95, 50],
    pmiss_labels=PMISS_LABELS)

# E2miss 4x1 (epp), data + models
ph.plot_emiss_4x1(
    pdf, [he_data] + sim_models_emiss, r'$E_{2,miss} [GeV]$',
    'h_E2miss_epp_SRC_pmiss',
    ylim=[60, 110, 110, 80], labely=[50, 95, 95, 65],
    pmiss_labels=PMISS_LABELS)

# ===========================================================================
# C.M. width studies vs Q^2
# ===========================================================================
pto.plotSig(pdf, Q2bins, f_Data_He, f_Sim_He_N2LO10, 'pcmx', -0.6, 0.6,
            r'$\sigma_{x,C.M.} [GeV]$',
            r'$\vec{p}_{C.M.} \cdot \hat{v}_{x}  [GeV]$', N0epp_array)
pto.plotSig(pdf, Q2bins, f_Data_He, f_Sim_He_N2LO10, 'pcmy', -0.6, 0.6,
            r'$\sigma_{y,C.M.} [GeV]$',
            r'$\vec{p}_{C.M.} \cdot \hat{v}_{y}  [GeV]$', N0epp_array)

ph.plot_q2_single(pdf, he_data, gcf, 'g_sigma_pcmx', r'$\sigma_{x,C.M.} [GeV]$')
ph.plot_q2_single(pdf, he_data, gcf, 'g_sigma_pcmy', r'$\sigma_{y,C.M.} [GeV]$')

# epp/ep vs Q^2 in pmiss bins
ph.plot_q2_2x2(pdf, he_data, gcf, 'h_Q2_epp_SRC_pmiss', r'$epp/ep$',
               PMISS_LABELS, ylim=(0, 0.25), label_xy=(2.3, 0.2), ylabel_size=15)

# Std dev of E_miss vs Q^2 in pmiss bins
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_sigma_E1miss_ep_pmiss',
               r'Std DeV $E_{1,miss} [GeV]$', PMISS_LABELS,
               ylim=(0, 0.15), label_xy=(2.3, 0.01), ylabel_size=10,
               big_label=r'$(e,e^{\prime}p)$')
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_sigma_E1miss_epp_pmiss',
               r'Std DeV $E_{1,miss} [GeV]$', PMISS_LABELS,
               ylim=(0, 0.15), label_xy=(2.3, 0.01), ylabel_size=10,
               big_label=r'$(e,e^{\prime}pp)$')
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_sigma_E2miss_epp_pmiss',
               r'Std DeV $E_{2,miss} [GeV]$', PMISS_LABELS,
               ylim=(0, 0.15), label_xy=(2.3, 0.01), ylabel_size=10,
               big_label=r'$(e,e^{\prime}pp)$')


# Mean of E_miss vs Q^2 in pmiss bins
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_mean_E1miss_ep_pmiss',
               r'Mean $E_{1,miss} [GeV]$', PMISS_LABELS,
               ylim=(-.02, .3), label_xy=(2.3, 0.04), ylabel_size=10,
               big_label=r'$(e,e^{\prime}p)$')
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_mean_E1miss_epp_pmiss', 
                r'Mean $E_{1,miss} [GeV]$', PMISS_LABELS,
                ylim=(-.02, .3), label_xy=(2.3, 0.04), ylabel_size=10,
                big_label=r'$(e,e^{\prime}pp)$')
ph.plot_q2_2x2(pdf, he_data, gcf, 'g_mean_E2miss_epp_pmiss',
               r'Mean $E_{2,miss} [GeV]$', PMISS_LABELS,
               ylim=(-.02, .3), label_xy=(2.3, 0.04), ylabel_size=10,
               big_label=r'$(e,e^{\prime}pp)$')


# ===========================================================================
# Cross-nucleus comparisons (He reference + scaled nuclei)
# ===========================================================================

# C.M. projections (48Ca omitted, matching original); Sn label sits in the
# slot 48Ca would have occupied.
nuc_cm_annotations = [
    {'x': -0.4, 'y': 200, 'text': r'${}^{4}He$',   'color': 'black'},
    {'x': -0.4, 'y': 240, 'text': r'${}^{12}C$',   'color': 'darkred'},
    {'x': -0.4, 'y': 280, 'text': r'${}^{40}Ca$',  'color': 'darkgreen'},
    {'x': -0.4, 'y': 320, 'text': r'${}^{120}Sn$', 'color': 'darkgoldenrod'},
]
for axis in ('x', 'y'):
    ph.plot_overlay(
        pdf, 'h_pcm{0}_epp'.format(axis), [he_data] + nuclei_no48,
        xlabel=r'$\vec{p}_{C.M.} \cdot \hat{v}_{%s}  [GeV]$' % axis,
        ylabel=r'Counts', xlim=(-0.6, 0.6), ylim=(0.0, 400),
        annotations=nuc_cm_annotations)

# epp/ep vs pMiss, all nuclei (ratio histogram -> unit_scale)
ph.plot_overlay(
    pdf, 'pMiss_epp', [he_data] + nuclei,
    xlabel=r'$p_{miss}$', ylabel=r'$epp/ep$',
    xlim=(0.4, 1.05), ylim=(0, 0.2), unit_scale=True,
    xlabel_size=25, ylabel_size=25,
    annotations=[
        {'x': 0.5, 'y': 0.1,  'text': r'${}^{4}He$',   'color': 'black'},
        {'x': 0.5, 'y': 0.12, 'text': r'${}^{12}C$',   'color': 'darkred'},
        {'x': 0.5, 'y': 0.14, 'text': r'${}^{40}Ca$',  'color': 'darkgreen'},
        {'x': 0.5, 'y': 0.16, 'text': r'${}^{48}Ca$',  'color': 'darkblue'},
        {'x': 0.5, 'y': 0.18, 'text': r'${}^{120}Sn$', 'color': 'darkgoldenrod'},
    ])

# E1miss 4x2, all nuclei (narrow offsets via offset_scale=0.4)
ph.plot_emiss_4x2(
    pdf, [he_data] + nuclei, r'$E_{1,miss} [GeV]$',
    'h_E1miss_ep_SRC_pmiss', 'h_E1miss_epp_SRC_pmiss',
    ylim_ep=[1000, 1500, 1500, 1000], ylim_epp=[60, 110, 150, 100],
    labely_ep=[850, 1250, 1250, 850], labely_epp=[50, 95, 130, 85],
    pmiss_labels=PMISS_LABELS, offset_scale=0.4)

# E2miss 4x1, all nuclei
ph.plot_emiss_4x1(
    pdf, [he_data] + nuclei, r'$E_{2,miss} [GeV]$',
    'h_E2miss_epp_SRC_pmiss',
    ylim=[30, 110, 130, 100], labely=[25, 95, 110, 85],
    pmiss_labels=PMISS_LABELS, offset_scale=0.4)

pdf.close()