"""
plot_primitives.py

The matplotlib drawing half of what plotTOOL.py's plotTGEStep / plotTGELine
/ plotTGEStepEx / plotTGEStepQ2 / plotTGELineQ2 used to do in one function
each (open file, unpack TGraph, plot). Now that graph_io.read_graph() does
all the ROOT unpacking, these functions only take plain (x, y, yerr) lists
and never touch a ROOT object -- any scaling/x-shift is applied by the
caller (plot_helpers.draw) before calling these, so these stay simple.

plotTOOL.py itself is left alone (Sys_Err.py still uses it); this is a
parallel, smaller module for the new driver scripts.
"""


def step_with_error(ax, x, y, yerr_low, yerr_high, color):
    """Data-style: step line + plain error bars (no marker)."""
    ax.errorbar(x, y, yerr=[yerr_low, yerr_high], color=color, linestyle='')
    ax.step(x, y, color=color, where='mid')


def errorbar_marker(ax, x, y, yerr_low, yerr_high, color, marker):
    """Scaled/offset data-style point set with a marker (no step line)."""
    ax.errorbar(x, y, yerr=[yerr_low, yerr_high], color=color,
                linestyle='', marker=marker)


def line_with_band(ax, x, y, yerr_low, yerr_high, color, alpha=0.3):
    """Simulation-style: a line with a shaded error band."""
    ylow = [yi - ei for yi, ei in zip(y, yerr_low)]
    yhigh = [yi + ei for yi, ei in zip(y, yerr_high)]
    ax.plot(x, y, color=color)
    ax.fill_between(x, ylow, yhigh, alpha=alpha, facecolor=color)


def _pad_for_step(x, y, yerr_low, yerr_high, xmin=0.0, xmax=10.0):
    """Q2-panel convention: repeat the first/last point out to wide x limits
    so the step line/band fills the whole panel instead of stopping short."""
    xp = [xmin] + list(x) + [xmax]
    yp = [y[0]] + list(y) + [y[-1]]
    elow = [yerr_low[0]] + list(yerr_low) + [yerr_low[-1]]
    ehigh = [yerr_high[0]] + list(yerr_high) + [yerr_high[-1]]
    return xp, yp, elow, ehigh


def step_with_error_q2(ax, x, y, yerr_low, yerr_high, color):
    xp, yp, elow, ehigh = _pad_for_step(x, y, yerr_low, yerr_high)
    ax.errorbar(xp, yp, yerr=[elow, ehigh], color=color, linestyle='')
    ax.step(xp, yp, color=color, where='mid')


def line_with_band_q2(ax, x, y, yerr_low, yerr_high, color, alpha=0.3):
    xp, yp, elow, ehigh = _pad_for_step(x, y, yerr_low, yerr_high)
    line_with_band(ax, xp, yp, elow, ehigh, color, alpha=alpha)
