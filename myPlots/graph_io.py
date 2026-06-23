"""
graph_io.py

The only module that talks directly to ROOT objects. Every other plotting
module works with plain python lists, never with a TGraph/TFile directly.
This is what replaces the ~15 near-duplicate "open file, loop over
TGraph::GetX()/GetY()/GetErrorY(), then plot" functions in plotTOOL.py: the
ROOT-unpacking happens once, here, and everything downstream just gets
(x, y, yerr) lists.
"""

import ROOT


def open_file(path):
    """Thin wrapper so driver scripts don't import ROOT directly."""
    f = ROOT.TFile.Open(path)
    if (not f) or f.IsZombie():
        raise IOError("Could not open ROOT file: %s" % path)
    return f


def read_graph(root_file, name):
    """Reads a TGraphErrors written by BuildGraphs.cpp.

    Returns (x, y, yerr) as plain python lists. yerr is whatever error was
    stored on the graph (BuildGraphs.cpp combines stat+sys in quadrature
    before writing, so this is already the combined error).

    Raises KeyError if the named graph isn't in the file -- callers should
    let this surface rather than silently plotting nothing, since a missing
    graph almost always means a naming mismatch with graph_names.py.
    """
    graph_dir = root_file.Get("graphs")
    if not graph_dir:
        raise KeyError("No 'graphs' directory in file %s -- did you run "
                        "BuildGraphs on it?" % root_file.GetName())

    graph = graph_dir.Get(name)
    if not graph:
        raise KeyError("Graph '%s' not found in %s" % (name, root_file.GetName()))

    x, y, yerr = [], [], []
    xvals = graph.GetX()
    yvals = graph.GetY()
    for i in range(graph.GetN()):
        x.append(xvals[i])
        y.append(yvals[i])
        yerr.append(graph.GetErrorY(i))
    return x, y, yerr
