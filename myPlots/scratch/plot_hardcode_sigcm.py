import matplotlib.pyplot as plt

# Column labels
columns = ["Lead FD\nRec CD", "Lead CD\nRec FD", "Lead CD\nRec CD"]
x = [0, 1, 2]

# Hardcoded values: (value, error)
data_x = [(0.1608, 0.0076), (0.1683, 0.0134), (0.1729, 0.0058)]
data_y = [(0.2009, 0.0156), (0.1820, 0.0182), (0.1726, 0.0057)]
sim_x  = [(0.1470, 0.0067), (0.1613, 0.0076), (0.1837, 0.0075)]
sim_y  = [(0.1448, 0.0056), (0.1503, 0.0059), (0.1568, 0.0047)]

# Fit values only exist for first two columns
fit_x = [0, 1]
fit_y = [(0.185, 0.030), (0.147, 0.018)]

def vals(pairs):
    return [p[0] for p in pairs], [p[1] for p in pairs]

fig, ax = plt.subplots(figsize=(8, 5))

# small horizontal offsets so points don't overlap
offset = 0.08
v, e = vals(data_x)
ax.errorbar([i - 1.5*offset for i in x], v, yerr=e, fmt='o', label='Data x', color='black')

v, e = vals(data_y)
ax.errorbar([i - 0.5*offset for i in x], v, yerr=e, fmt='s', label='Data y', color='dimgray')

v, e = vals(sim_x)
ax.errorbar([i + 0.5*offset for i in x], v, yerr=e, fmt='o', label='Sim x', color='red')

v, e = vals(sim_y)
ax.errorbar([i + 1.5*offset for i in x], v, yerr=e, fmt='s', label='Sim y', color='salmon')

v, e = vals(fit_y)
ax.errorbar([i + 2.5*offset for i in fit_x], v, yerr=e, fmt='^', label='Fit y', color='dodgerblue')

ax.set_xticks(x)
ax.set_xticklabels(columns)
ax.set_ylabel(r'$\sigma_{CM}$')
ax.set_title(r'$\sigma_{CM}$ values with error bars')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/sigma_cm_plot.png', dpi=150)
plt.show()