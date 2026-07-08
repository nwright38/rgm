# import matplotlib.pyplot as plt

# # ------------------------------------------------------------------
# # Options
# # ------------------------------------------------------------------
# SHOW_LITERATURE_DATA = True   # set to False to hide the literature comparison points

# # ------------------------------------------------------------------
# # This Work: hard-coded fit results, sorted by increasing A
# # ------------------------------------------------------------------
# A = [4, 12, 40, 48, 120]

# # p_x (mean, sigma) with errors [GeV]
# px_mu     = [0.0155, 0.0168, 0.0180, 0.0020, -0.0290]
# px_mu_err = [0.0038, 0.0036, 0.0055, 0.0125, 0.0178]
# px_sig    = [0.1345, 0.1580, 0.1583, 0.1836, 0.1766]
# px_sig_err= [0.0034, 0.0037, 0.0057, 0.0148, 0.0205]

# # p_y (mean, sigma) with errors [GeV]
# py_mu     = [0.0015, 0.0021, -0.0024, 0.0193, 0.0190]
# py_mu_err = [0.0044, 0.0054, 0.0093, 0.0201, 0.0296]
# py_sig    = [0.1473, 0.1980, 0.2143, 0.2338, 0.2366]
# py_sig_err= [0.0045, 0.0069, 0.0131, 0.0301, 0.0303]

# # ------------------------------------------------------------------
# # Literature comparison points (sigma_c.m., in MeV/c), digitized
# # by eye from the reference figure. Converted to GeV/c for overlay.
# # ------------------------------------------------------------------
# # "This Work" (red circles) from reference plot: C, Al, Fe, Pb
# ref_thiswork_A       = [12, 27, 56, 208]
# ref_thiswork_sig     = [143, 155, 160, 158]      # MeV/c
# ref_thiswork_sig_err = [6, 10, 8, 15]            # MeV/c

# # BNL (p,2pn) - blue square, at Carbon
# ref_bnl_A       = [12]
# ref_bnl_sig     = [143]
# ref_bnl_sig_err = [12]

# # Hall-A (e,e'pp) - blue inverted triangle, at Carbon
# ref_hallA_pp_A       = [12]
# ref_hallA_pp_sig     = [138]
# ref_hallA_pp_sig_err = [21]

# # Hall-A (e,e'pn) - blue triangle, at 4He
# ref_hallA_pn_A       = [4]
# ref_hallA_pn_sig     = [100]
# ref_hallA_pn_sig_err = [20]

# # Ciofi and Simula (open stars): 4He, C, Al, Fe, Pb
# ref_ciofi_A       = [4, 12, 27, 56, 208]
# ref_ciofi_sig     = [90, 138, 143, 130, 150]
# ref_ciofi_sig_err = [10, 0, 0, 0, 10]

# MEV_TO_GEV = 1.0e-3

# def to_gev(vals):
#     return [v * MEV_TO_GEV for v in vals]

# # ------------------------------------------------------------------
# # Simulation point at A = 12 (shown by default, distinct marker)
# # ------------------------------------------------------------------
# sim_A         = 12
# sim_px_mu     = 0.0153
# sim_px_mu_err = 0.0011
# sim_px_sig    = 0.1383
# sim_px_sig_err= 0.0010
# sim_py_mu     = -0.0070
# sim_py_mu_err = 0.0011
# sim_py_sig    = 0.1384
# sim_py_sig_err= 0.0010

# # ------------------------------------------------------------------
# # Plot 1: sigma_x and sigma_y vs A (points/error bars only, no lines)
# # ------------------------------------------------------------------
# plt.figure(figsize=(7, 5))

# plt.errorbar(A, px_sig, yerr=px_sig_err, fmt='o', color='black',
#              label=r'$\sigma_{p_{C.M.,x}}$ (This work)', capsize=3)
# plt.errorbar(A, py_sig, yerr=py_sig_err, fmt='s', color='royalblue',
#              label=r'$\sigma_{p_{C.M.,y}}$ (This work)', capsize=3)

# if SHOW_LITERATURE_DATA:
#     plt.errorbar(ref_thiswork_A, to_gev(ref_thiswork_sig), yerr=to_gev(ref_thiswork_sig_err),
#                  fmt='o', color='red', label=r'$\sigma_{c.m.}$ This Work (ref.)', capsize=3)
#     plt.errorbar(ref_bnl_A, to_gev(ref_bnl_sig), yerr=to_gev(ref_bnl_sig_err),
#                  fmt='s', color='green', label='BNL (p,2pn)', capsize=3)
#     plt.errorbar(ref_hallA_pp_A, to_gev(ref_hallA_pp_sig), yerr=to_gev(ref_hallA_pp_sig_err),
#                  fmt='v', color='green', label="Hall-A (e,e'pp)", capsize=3)
#     plt.errorbar(ref_hallA_pn_A, to_gev(ref_hallA_pn_sig), yerr=to_gev(ref_hallA_pn_sig_err),
#                  fmt='^', color='green', label="Hall-A (e,e'pn)", capsize=3)
#     plt.errorbar(ref_ciofi_A, to_gev(ref_ciofi_sig), yerr=to_gev(ref_ciofi_sig_err),
#                  fmt='*', color='black', markerfacecolor='none', markersize=10,
#                  label='Ciofi and Simula', capsize=3)

# # Simulation point at A = 12 (shown by default, distinct marker: diamond)
# plt.errorbar(sim_A, sim_px_sig, yerr=sim_px_sig_err, fmt='D', color='black',
#              markerfacecolor='none', markersize=8,
#              label=r'$\sigma_{p_{C.M.,x}}$ (Simulation)', capsize=3)
# plt.errorbar(sim_A, sim_py_sig, yerr=sim_py_sig_err, fmt='D', color='royalblue',
#              markerfacecolor='none', markersize=8,
#              label=r'$\sigma_{p_{C.M.,y}}$ (Simulation)', capsize=3)

# plt.xscale('log')
# plt.xlabel('A')
# plt.ylabel(r'$\sigma$ [GeV]')
# plt.title(r'$p_{C.M.}$ width vs nucleus size A')
# plt.legend(fontsize=8)
# plt.tight_layout()
# plt.savefig('pdf/scratch/sigma_vs_A.pdf')

# # ------------------------------------------------------------------
# # Plot 2: mu_x and mu_y vs A (points/error bars only, no lines)
# # ------------------------------------------------------------------
# plt.figure(figsize=(7, 5))
# plt.errorbar(A, px_mu, yerr=px_mu_err, fmt='o', color='black',
#              label=r'$\mu_{p_{C.M.,x}}$', capsize=3)
# plt.errorbar(A, py_mu, yerr=py_mu_err, fmt='s', color='royalblue',
#              label=r'$\mu_{p_{C.M.,y}}$', capsize=3)

# # Simulation point at A = 12 (shown by default, distinct marker: diamond)
# plt.errorbar(sim_A, sim_px_mu, yerr=sim_px_mu_err, fmt='D', color='black',
#              markerfacecolor='none', markersize=8,
#              label=r'$\mu_{p_{C.M.,x}}$ (Simulation)', capsize=3)
# plt.errorbar(sim_A, sim_py_mu, yerr=sim_py_mu_err, fmt='D', color='royalblue',
#              markerfacecolor='none', markersize=8,
#              label=r'$\mu_{p_{C.M.,y}}$ (Simulation)', capsize=3)

# plt.xscale('log')
# plt.xlabel('A')
# plt.ylabel(r'$\mu$ [GeV]')
# plt.title(r'$p_{C.M.}$ mean vs nucleus size A')
# plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
# plt.legend(fontsize=8)
# plt.tight_layout()
# plt.savefig('pdf/scratch/mean_vs_A.pdf')

# # plt.show()


import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Options
# ------------------------------------------------------------------
SHOW_LITERATURE_DATA = False   # set to False to hide the literature comparison points

# ------------------------------------------------------------------
# This Work: hard-coded fit results, sorted by increasing A
# ------------------------------------------------------------------
A = [4, 12, 40, 48, 120]

# p_x (mean, sigma) with errors [GeV]
px_mu     = [0.0085, 0.0279, 0.0436, 0.0287, 0.0190]
px_mu_err = [0.0018, 0.0028, 0.0037, 0.0064, 0.0105]
px_sig    = [0.1285, 0.1602, 0.1664, 0.1664, 0.1779]
px_sig_err= [0.0021, 0.0037, 0.0039, 0.0067, 0.0101]

# p_y (mean, sigma) with errors [GeV]
py_mu     = [-0.0071, -0.0122, -0.0172, -0.0110, -0.0063]
py_mu_err = [0.0016, 0.0025, 0.0037, 0.0059, 0.0103]
py_sig    = [0.1181, 0.1560, 0.1685, 0.1587, 0.1674]
py_sig_err= [0.0017, 0.0034, 0.0040, 0.0058, 0.0114]

# ------------------------------------------------------------------
# Literature comparison points (sigma_c.m., in MeV/c), digitized
# by eye from the reference figure. Converted to GeV/c for overlay.
# ------------------------------------------------------------------
# "This Work" (red circles) from reference plot: C, Al, Fe, Pb
ref_thiswork_A       = [12, 27, 56, 208]
ref_thiswork_sig     = [143, 155, 160, 158]      # MeV/c
ref_thiswork_sig_err = [6, 10, 8, 15]            # MeV/c

# BNL (p,2pn) - blue square, at Carbon
ref_bnl_A       = [12]
ref_bnl_sig     = [143]
ref_bnl_sig_err = [12]

# Hall-A (e,e'pp) - blue inverted triangle, at Carbon
ref_hallA_pp_A       = [12]
ref_hallA_pp_sig     = [138]
ref_hallA_pp_sig_err = [21]

# Hall-A (e,e'pn) - blue triangle, at 4He
ref_hallA_pn_A       = [4]
ref_hallA_pn_sig     = [100]
ref_hallA_pn_sig_err = [20]

# Ciofi and Simula (open stars): 4He, C, Al, Fe, Pb
ref_ciofi_A       = [4, 12, 27, 56, 208]
ref_ciofi_sig     = [90, 138, 143, 130, 150]
ref_ciofi_sig_err = [10, 0, 0, 0, 10]

MEV_TO_GEV = 1.0e-3

def to_gev(vals):
    return [v * MEV_TO_GEV for v in vals]

# ------------------------------------------------------------------
# Simulation point at A = 12 (shown by default, distinct marker)
# ------------------------------------------------------------------
sim_A         = 12
sim_px_mu     = 0.00985466
sim_px_mu_err = 0.00165
sim_px_sig    = 0.140432 
sim_px_sig_err= 0.00196
sim_py_mu     = -0.005986
sim_py_mu_err = 0.0015
sim_py_sig    = 0.13751
sim_py_sig_err= 0.00184

# ------------------------------------------------------------------
# Plot 1: sigma_x and sigma_y vs A (points/error bars only, no lines)
# ------------------------------------------------------------------
plt.figure(figsize=(7, 5))

plt.errorbar(A, px_sig, yerr=px_sig_err, fmt='o', color='black',
             label=r'$\sigma_{p_{C.M.,x}}$ (This work)', capsize=3)
plt.errorbar(A, py_sig, yerr=py_sig_err, fmt='s', color='royalblue',
             label=r'$\sigma_{p_{C.M.,y}}$ (This work)', capsize=3)

if SHOW_LITERATURE_DATA:
    plt.errorbar(ref_thiswork_A, to_gev(ref_thiswork_sig), yerr=to_gev(ref_thiswork_sig_err),
                 fmt='o', color='red', label=r'$\sigma_{c.m.}$ This Work (ref.)', capsize=3)
    plt.errorbar(ref_bnl_A, to_gev(ref_bnl_sig), yerr=to_gev(ref_bnl_sig_err),
                 fmt='s', color='green', label='BNL (p,2pn)', capsize=3)
    plt.errorbar(ref_hallA_pp_A, to_gev(ref_hallA_pp_sig), yerr=to_gev(ref_hallA_pp_sig_err),
                 fmt='v', color='green', label="Hall-A (e,e'pp)", capsize=3)
    plt.errorbar(ref_hallA_pn_A, to_gev(ref_hallA_pn_sig), yerr=to_gev(ref_hallA_pn_sig_err),
                 fmt='^', color='green', label="Hall-A (e,e'pn)", capsize=3)
    plt.errorbar(ref_ciofi_A, to_gev(ref_ciofi_sig), yerr=to_gev(ref_ciofi_sig_err),
                 fmt='*', color='orange', markerfacecolor='none', markersize=10,
                 label='Ciofi and Simula', capsize=3)

# Simulation point at A = 12 (shown by default, distinct marker: diamond)
plt.errorbar(sim_A, sim_px_sig, yerr=sim_px_sig_err, fmt='D', color='black',
             markerfacecolor='none', markersize=8,
             label=r'$\sigma_{p_{C.M.,x}}$ (Simulation)', capsize=3)
plt.errorbar(sim_A, sim_py_sig, yerr=sim_py_sig_err, fmt='D', color='royalblue',
             markerfacecolor='none', markersize=8,
             label=r'$\sigma_{p_{C.M.,y}}$ (Simulation)', capsize=3)

plt.xscale('log')
plt.xlabel('A')
plt.ylabel(r'$\sigma$ [GeV]')
plt.title(r'$p_{C.M.}$ width vs nucleus size A')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig('pdf/scratch/sigma_vs_A.pdf')

# ------------------------------------------------------------------
# Plot 2: mu_x and mu_y vs A (points/error bars only, no lines)
# ------------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.errorbar(A, px_mu, yerr=px_mu_err, fmt='o', color='black',
             label=r'$\mu_{p_{C.M.,x}}$', capsize=3)
plt.errorbar(A, py_mu, yerr=py_mu_err, fmt='s', color='royalblue',
             label=r'$\mu_{p_{C.M.,y}}$', capsize=3)

# Simulation point at A = 12 (shown by default, distinct marker: diamond)
plt.errorbar(sim_A, sim_px_mu, yerr=sim_px_mu_err, fmt='D', color='black',
             markerfacecolor='none', markersize=8,
             label=r'$\mu_{p_{C.M.,x}}$ (Simulation)', capsize=3)
plt.errorbar(sim_A, sim_py_mu, yerr=sim_py_mu_err, fmt='D', color='royalblue',
             markerfacecolor='none', markersize=8,
             label=r'$\mu_{p_{C.M.,y}}$ (Simulation)', capsize=3)

plt.xscale('log')
plt.xlabel('A')
plt.ylabel(r'$\mu$ [GeV]')
plt.title(r'$p_{C.M.}$ mean vs nucleus size A')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig('pdf/scratch/mean_vs_A.pdf')

# plt.show()

