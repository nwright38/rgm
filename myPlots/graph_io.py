"""
graph_io.py

The only module that reads ROOT files directly. Every other plotting
module works with plain python lists, never with a TGraph/TFile object
directly. This is what replaces the ~15 near-duplicate "open file, loop
over TGraph::GetX()/GetY()/GetErrorY(), then plot" functions in the old
plotTOOL.py: the ROOT-unpacking happens once, here, and everything
downstream just gets (x, y, yerr) lists.

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


def read_graph(root_file, name):
    """Reads a TGraphErrors written by BuildGraphs.cpp, from the "graphs"
    TDirectory it created.

    Returns (x, y, yerr) as plain python lists. yerr is whatever error was
    stored on the graph (BuildGraphs.cpp combines stat+sys in quadrature
    before writing, so this is already the combined error).

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

    # TGraphErrors' uproot model doesn't mix in TGraph's values() helper
    # (confirmed: Model_TGraphErrors_v3 only inherits the empty TGraphErrors
    # behavior class, not TGraph's), so read the members directly. member()
    # searches the embedded TGraph base for fX/fY by default.
    x = graph.member("fX")
    y = graph.member("fY")
    yerr = graph.member("fEY")
    return list(x), list(y), list(yerr)
