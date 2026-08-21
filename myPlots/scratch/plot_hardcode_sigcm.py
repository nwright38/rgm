# import numpy as np
# import matplotlib.pyplot as plt

# # Column labels (4th column added for the weighted average)
# columns = ["FD\nRec CD", "CD\nRec FD", "CD\nRec CD", "Weighted\nAverage"]
# x = [0, 1, 2, 3]

# # Hardcoded values: (value, error)
# data_x = [(0.1608, 0.0076), (0.1683, 0.0134), (0.1729, 0.0058)]
# data_y = [(0.2009, 0.0156), (0.1820, 0.0182), (0.1726, 0.0057)]
# sim_x  = [(0.1470, 0.0067), (0.1613, 0.0076), (0.1837, 0.0075)]
# sim_y  = [(0.1448, 0.0056), (0.1503, 0.0059), (0.1568, 0.0047)]

# # Fit values only exist for the first two columns -- not included in the
# # weighted average (only 2 points, and it's a derived fit result rather
# # than an independent reconstruction situation)
# fit_x_idx = [0, 1]
# fit_y = [(0.185, 0.030), (0.147, 0.018)]

# def vals(pairs):
#     return [p[0] for p in pairs], [p[1] for p in pairs]

# def weighted_avg_with_sys(pairs):
#     """
#     Inverse-variance weighted mean of the input (value, stat_error) pairs,
#     with:
#       - sigma_stat: standard error on the weighted mean (purely statistical)
#       - sigma_sys : unbiased weighted sample variance of the points about
#                     the mean (reliability-weights form), i.e. the genuine
#                     spread between the different reconstruction situations,
#                     weighted by how precisely each situation is known.
#       - sigma_tot : the two added in quadrature.
#     """
#     v = np.array([p[0] for p in pairs], dtype=float)
#     e = np.array([p[1] for p in pairs], dtype=float)
#     w = 1.0 / e**2
#     sw = w.sum()

#     mean = (w * v).sum() / sw
#     sigma_stat = np.sqrt(1.0 / sw)

#     denom = sw**2 - (w**2).sum()
#     sigma_sys2 = (sw / denom) * (w * (v - mean)**2).sum()
#     sigma_sys = np.sqrt(sigma_sys2)

#     sigma_tot = np.sqrt(sigma_stat**2 + sigma_sys2)
#     return mean, sigma_stat, sigma_sys, sigma_tot

# # Weighted average across the three situations, per series
# mean_dx, stat_dx, sys_dx, tot_dx = weighted_avg_with_sys(data_x)
# mean_dy, stat_dy, sys_dy, tot_dy = weighted_avg_with_sys(data_y)
# mean_sx, stat_sx, sys_sx, tot_sx = weighted_avg_with_sys(sim_x)
# mean_sy, stat_sy, sys_sy, tot_sy = weighted_avg_with_sys(sim_y)

# print(f"Data x weighted avg = {mean_dx:.4f}  (stat {stat_dx:.4f}, sys {sys_dx:.4f}, total {tot_dx:.4f})")
# print(f"Data y weighted avg = {mean_dy:.4f}  (stat {stat_dy:.4f}, sys {sys_dy:.4f}, total {tot_dy:.4f})")
# print(f"Sim x  weighted avg = {mean_sx:.4f}  (stat {stat_sx:.4f}, sys {sys_sx:.4f}, total {tot_sx:.4f})")
# print(f"Sim y  weighted avg = {mean_sy:.4f}  (stat {stat_sy:.4f}, sys {sys_sy:.4f}, total {tot_sy:.4f})")

# fig, ax = plt.subplots(figsize=(8, 5))

# # small horizontal offsets so points don't overlap
# offset = 0.08
# situ_idx = [0, 1, 2]  # the three reconstruction-situation columns

# v, e = vals(data_x)
# ax.errorbar([i - 1.5*offset for i in situ_idx], v, yerr=e, fmt='o', label='Data x', color='black')

# v, e = vals(data_y)
# ax.errorbar([i - 0.5*offset for i in situ_idx], v, yerr=e, fmt='s', label='Data y', color='dimgray')

# v, e = vals(sim_x)
# ax.errorbar([i + 0.5*offset for i in situ_idx], v, yerr=e, fmt='o', label='Sim x', color='red')

# v, e = vals(sim_y)
# ax.errorbar([i + 1.5*offset for i in situ_idx], v, yerr=e, fmt='s', label='Sim y', color='salmon')

# v, e = vals(fit_y)
# ax.errorbar([i + 2.5*offset for i in fit_x_idx], v, yerr=e, fmt='^', label='Fit y', color='dodgerblue')

# # Weighted average column (total error = stat + sys added in quadrature)
# avg_idx = 3
# ax.errorbar([avg_idx - 1.5*offset], [mean_dx], yerr=[tot_dx], fmt='D', 
#             color='black', markerfacecolor='white', capsize=4)
# ax.errorbar([avg_idx - 0.5*offset], [mean_dy], yerr=[tot_dy], fmt='D',
#             color='dimgray', markerfacecolor='white', capsize=4)
# ax.errorbar([avg_idx + 0.5*offset], [mean_sx], yerr=[tot_sx], fmt='D',
#             color='red', markerfacecolor='white', capsize=4)
# ax.errorbar([avg_idx + 1.5*offset], [mean_sy], yerr=[tot_sy], fmt='D',
#             color='salmon', markerfacecolor='white', capsize=4)

# ax.set_xticks(x)
# ax.set_xticklabels(columns)
# ax.set_ylabel(r'$\sigma_{CM}$')
# ax.set_title(r'$\sigma_{CM}$ values with error bars')
# ax.legend()
# ax.grid(True, alpha=0.3)

# plt.tight_layout()
# plt.savefig('sigma_cm_plot.pdf', dpi=150)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Column labels (4th column added for the weighted average)
columns = ["Lead FD\nRec CD", "Lead CD\nRec FD", "Lead CD\nRec CD", "Weighted\nAverage"]
x = [0, 1, 2, 3]

# Hardcoded values: (value, error)
data_x = [(0.1303, 0.0093), (0.1427, 0.0158), (0.1288, 0.0040)]
data_y = [(0.1106, 0.0061), (0.1253, 0.0111), (0.1159, 0.0031)]

# Simulation values only exist for the first column (Lead FD / Rec CD)
sim_x_idx = [0]
sim_x = [(0.1385, 0.0031)]
sim_y_idx = [0]
sim_y = [(0.1359, 0.0029)]

def vals(pairs):
    return [p[0] for p in pairs], [p[1] for p in pairs]

def weighted_avg_with_sys(pairs):
    """
    Inverse-variance weighted mean of the input (value, stat_error) pairs,
    with:
      - sigma_stat: standard error on the weighted mean (purely statistical)
      - sigma_sys : unbiased weighted sample variance of the points about
                    the mean (reliability-weights form), i.e. the genuine
                    spread between the different reconstruction situations,
                    weighted by how precisely each situation is known.
      - sigma_tot : the two added in quadrature.
    """
    v = np.array([p[0] for p in pairs], dtype=float)
    e = np.array([p[1] for p in pairs], dtype=float)
    w = 1.0 / e**2
    sw = w.sum()

    mean = (w * v).sum() / sw
    sigma_stat = np.sqrt(1.0 / sw)

    denom = sw**2 - (w**2).sum()
    sigma_sys2 = (sw / denom) * (w * (v - mean)**2).sum()
    sigma_sys = np.sqrt(sigma_sys2)

    sigma_tot = np.sqrt(sigma_stat**2 + sigma_sys2)
    return mean, sigma_stat, sigma_sys, sigma_tot

# Weighted average of the three Data points, per axis
mean_x, stat_x, sys_x, tot_x = weighted_avg_with_sys(data_x)
mean_y, stat_y, sys_y, tot_y = weighted_avg_with_sys(data_y)

print(f"sigma_CM,x weighted avg = {mean_x:.4f}  (stat {stat_x:.4f}, sys {sys_x:.4f}, total {tot_x:.4f})")
print(f"sigma_CM,y weighted avg = {mean_y:.4f}  (stat {stat_y:.4f}, sys {sys_y:.4f}, total {tot_y:.4f})")

fig, ax = plt.subplots(figsize=(8, 5))

# small horizontal offsets so points don't overlap
offset = 0.08
data_idx = [0, 1, 2]  # only the first three columns have per-situation Data points

v, e = vals(data_x)
ax.errorbar([i - 0.5*offset for i in data_idx], v, yerr=e, fmt='o', label='Data x', color='black')

v, e = vals(data_y)
ax.errorbar([i + 0.5*offset for i in data_idx], v, yerr=e, fmt='s', label='Data y', color='dimgray')

v, e = vals(sim_x)
ax.errorbar([i - 0.5*offset for i in sim_x_idx], v, yerr=e, fmt='o', label='Sim x', color='red')

v, e = vals(sim_y)
ax.errorbar([i + 0.5*offset for i in sim_y_idx], v, yerr=e, fmt='s', label='Sim y', color='salmon')

# Weighted average column (total error = stat + sys added in quadrature)
avg_idx = 3
ax.errorbar([avg_idx - 0.5*offset], [mean_x], yerr=[tot_x], fmt='D',
            color='black', markerfacecolor='white', capsize=4)
ax.errorbar([avg_idx + 0.5*offset], [mean_y], yerr=[tot_y], fmt='D',
            color='dimgray', markerfacecolor='white', capsize=4)

ax.set_xticks(x)
ax.set_xticklabels(columns)
ax.set_ylabel(r'$\sigma_{CM}$')
ax.set_title(r'$\sigma_{CM}$ values with error bars (Lead)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('sigma_cm_plot.pdf', dpi=150)
# plt.show()