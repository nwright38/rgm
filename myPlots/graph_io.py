"""
graph_io.py

The only module that reads ROOT files directly. Every other plotting
module works with plain python lists, never with a TGraph/TFile object
directly. This is what replaces the ~15 near-duplicate "open file, loop
over TGraph::GetX()/GetY()/GetErrorY(), then plot" functions in the old
plotTOOL.py: the ROOT-unpacking happens once, here, and everything
downstream just gets plain python arrays.

Uses uproot, not PyROOT -- PyROOT isn't installed in this environment
(confirmed: `import ROOT` fails, but `import uproot` works, matching the
BAND correlated_bg plotting scripts at
/Users/nataliewright/Desktop/MIT/BAND/plotting/correlated_bg, which read
ROOT files the same way). uproot is pure-python, so these scripts can run
locally without a ROOT installation at all -- only the C++ side
(Main_Figs_Binned.cpp, BuildGraphs.cpp) needs real ROOT, on the cluster.
"""

import uproot


def open_file(path):
    """Thin wrapper so driver scripts don't import uproot directly."""
    try:
        return uproot.open(path)
    except Exception as exc:
        raise IOError("Could not open ROOT file: %s (%s)" % (path, exc))


def read_graph_asymm(root_file, name):
    """Reads a TGraphErrors/TGraphAsymmErrors written by BuildGraphs.cpp.

    Returns (x, y, yerr_low, yerr_high). Symmetric graphs are promoted to
    asymmetric form with yerr_low == yerr_high.

    Raises KeyError if the named graph isn't in the file -- callers should
    let this surface rather than silently plotting nothing, since a missing
    graph almost always means a naming mismatch with graph_names.py.
    """
    try:
        graphs_dir = root_file["graphs"]
    except Exception:
        raise KeyError("No 'graphs' directory in file %s -- did you run "
                        "BuildGraphs on it?" % root_file.file_path)

    try:
        graph = graphs_dir[name]
    except Exception:
        raise KeyError("Graph '%s' not found in %s" % (name, root_file.file_path))

    # uproot's graph models don't all share the same convenience helpers, so
    # read the stored members directly.
    x = graph.member("fX")
    y = graph.member("fY")
    try:
        yerr_low = graph.member("fEYlow")
        yerr_high = graph.member("fEYhigh")
    except Exception:
        yerr = graph.member("fEY")
        yerr_low = yerr
        yerr_high = yerr
    return list(x), list(y), list(yerr_low), list(yerr_high)


def read_graph(root_file, name):
    """Returns (x, y, yerr) for callers that still need one symmetric error.

    For asymmetric graphs, yerr is max(yerr_low, yerr_high) per point so the
    fallback never understates the plotted uncertainty.
    """
    x, y, yerr_low, yerr_high = read_graph_asymm(root_file, name)
    yerr = [max(lo, hi) for lo, hi in zip(yerr_low, yerr_high)]
    return x, y, yerr
