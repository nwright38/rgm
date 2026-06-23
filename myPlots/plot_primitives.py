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


def step_with_error(ax, x, y, yerr, color):
    """Data-style: step line + plain error bars (no marker)."""
    ax.errorbar(x, y, yerr, color=color, linestyle='')
    ax.step(x, y, color=color, where='mid')


def errorbar_marker(ax, x, y, yerr, color, marker):
    """Scaled/offset data-style point set with a marker (no step line)."""
    ax.errorbar(x, y, yerr, color=color, linestyle='', marker=marker)


def line_with_band(ax, x, y, yerr, color, alpha=0.3):
    """Simulation-style: a line with a shaded +/- error band."""
    ylow = [yi - ei for yi, ei in zip(y, yerr)]
    yhigh = [yi + ei for yi, ei in zip(y, yerr)]
    ax.plot(x, y, color=color)
    ax.fill_between(x, ylow, yhigh, alpha=alpha, facecolor=color)


def _pad_for_step(x, y, yerr, xmin=0.0, xmax=10.0):
    """Q2-panel convention: repeat the first/last point out to wide x limits
    so the step line/band fills the whole panel instead of stopping short."""
    xp = [xmin] + list(x) + [xmax]
    yp = [y[0]] + list(y) + [y[-1]]
    ep = [yerr[0]] + list(yerr) + [yerr[-1]]
    return xp, yp, ep


def step_with_error_q2(ax, x, y, yerr, color):
    xp, yp, ep = _pad_for_step(x, y, yerr)
    ax.errorbar(xp, yp, ep, color=color, linestyle='')
    ax.step(xp, yp, color=color, where='mid')


def line_with_band_q2(ax, x, y, yerr, color, alpha=0.3):
    xp, yp, ep = _pad_for_step(x, y, yerr)
    line_with_band(ax, xp, yp, ep, color, alpha=alpha)
