"""
normalization.py

Computes the data/He (or sim/He) normalization scale factors used to put
different files on the same display scale in the overlay plots.

Replaces plotTOOL.py's getFactor_ep/getFactor_epp, which had a bug: their
loop over the 7 Q2 bins overwrote an accumulator each iteration instead of
summing, so the result only ever reflected the last bin. Here we instead
read the already-correctly-integrated Q2 yield straight from the
integratedTable graphs BuildGraphs.cpp wrote (which sums over all bins by
construction), via graph_io + ratio.py.
"""

import graph_io
import graph_names
import ratio


def _integrated_yield(root_file, task_name, selection):
    """Reads the single grand-total point for a one-selector-axis task whose
    selector axis is Q2 itself (e.g. 'Q2_ep_SRC_Q2' / 'Q2_epp_SRC_Q2')."""
    pattern = graph_names.pattern_collapse_axis('Q2')
    name = graph_names.integrated_graph_name(task_name, selection, pattern)
    x, y, yerr = graph_io.read_graph(root_file, name)
    if not y:
        raise ValueError("No integrated yield found for %s/%s in %s" %
                          (task_name, selection, root_file.file_path))
    return y[0], yerr[0]


def _legacy_last_q2_slice_yield(root_file, task_name, selection):
    """Emulates plotTOOL.getFactor_ep/getFactor_epp legacy behavior.

    The old helpers effectively used only the LAST Q2 selector slice (bin 6)
    when computing normalization factors.
    """
    # Legacy compatibility: depending on producer mode, the EP-gated
    # Q2_ep_SRC_Q2 task can be stored under either selection='ep'
    # (modern-consistent) or selection='epp' (legacy-mirroring quirk).
    # Try requested selection first, then fallback to the alternate key.
    candidates = [selection]
    if selection == 'ep':
        candidates.append('epp')
    elif selection == 'epp':
        candidates.append('ep')

    last_exc = None
    for sel in candidates:
        name = graph_names.diff_graph_name(task_name, sel, [6])
        try:
            _x, y, yerr = graph_io.read_graph(root_file, name)
            if y:
                return sum(y), (sum(e * e for e in yerr)) ** 0.5
        except Exception as exc:
            last_exc = exc

    raise KeyError("Legacy normalization graph not found for %s/%s in %s (tried selections: %s). Last error: %s" %
                   (task_name, selection, root_file.file_path, candidates, last_exc))


def scale_factor(target_file, reference_file, selection, mode='integrated'):
    """Returns (factor, error) such that target * factor displays on the
    same scale as reference, for the given selection ('ep' or 'epp')."""
    task_name = 'Q2_%s_SRC_Q2' % selection
    if mode == 'legacy-last-q2':
        target_y, target_err = _legacy_last_q2_slice_yield(target_file, task_name, selection)
        ref_y, ref_err = _legacy_last_q2_slice_yield(reference_file, task_name, selection)
    else:
        target_y, target_err = _integrated_yield(target_file, task_name, selection)
        ref_y, ref_err = _integrated_yield(reference_file, task_name, selection)
    return ratio.ratio_with_error(ref_y, ref_err, target_y, target_err)
