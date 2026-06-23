"""
graph_names.py

The one place that knows the graph-naming convention BuildGraphs.cpp uses
when it writes into the "graphs" TDirectory. Driver scripts should always
go through these functions instead of hand-building name strings, so a
naming-convention change only has to happen in one place (here and in
BuildGraphs.cpp).

Convention (must match BuildGraphs.cpp exactly):
  differential graph : "<task_name>|<selection>|<b0>_<b1>_..."   (or "|none" with no selector axes)
  integrated graph    : "<task_name>|<selection>|<pattern>"

`axis_bin` order matches the order the axes were listed in that task's
FillTask definition in Main_Figs_Binned.cpp (e.g. for
"E1miss_ep_SRC_pmiss_Q2", axis_bin[0] is the pMiss bin, axis_bin[1] is the
Q2 bin) -- this module doesn't re-derive that, callers need to know it from
the task they're asking for, same as they would have needed to know which
histogram array index meant what in the old pipeline.
"""

PATTERN_NONE = "collapse:none"   # no axis collapsed -- one point per selector bin
PATTERN_ALL = "collapse:all"     # every selector axis collapsed -- single grand-total point


def pattern_collapse_axis(axis_name):
    """Pattern for collapsing exactly one named axis (e.g. 'pMiss' or 'Q2')."""
    return "collapse:" + axis_name


def diff_graph_name(task_name, selection, axis_bin):
    """axis_bin: list of selector-axis bin indices, e.g. [2] or [1, 0], or
    [] for a task with no selector axes at all."""
    suffix = "_".join(str(b) for b in axis_bin) if axis_bin else "none"
    return "%s|%s|%s" % (task_name, selection, suffix)


def integrated_graph_name(task_name, selection, pattern):
    return "%s|%s|%s" % (task_name, selection, pattern)


def ratio_diff_graph_name(numerator_task, denominator_task, axis_bin):
    """Matches buildRatioDiffRows() in BinnedHistStore.h, which sets
    task_name = "<numerator>_over_<denominator>" and selection = "ratio"."""
    ratio_task = "%s_over_%s" % (numerator_task, denominator_task)
    return diff_graph_name(ratio_task, "ratio", axis_bin)
