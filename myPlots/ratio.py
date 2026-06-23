"""
ratio.py

One small, pure function for dividing two independent measurements and
propagating their error. Used for:
  - the data/sim ratio subpanel under every overlay plot
  - the A-dependence normalization factor (nucleus yield / He yield)

Data and simulation (and different nuclei) are independent samples, not one
nested inside the other -- unlike the C++ side's epp/ep ratio (where every
epp event is also an ep event), so the standard independent-quantities
formula applies directly; no toy-correlation handling is needed here.
"""


def ratio_with_error(y1, e1, y2, e2):
    """Returns (y1/y2, error) for one point. error is 0 if either value is 0
    (matches the convention already used elsewhere in this codebase for
    empty bins, rather than raising or returning NaN)."""
    if y1 == 0 or y2 == 0:
        return 0.0, 0.0
    r = y1 / y2
    rel = (e1 / y1) ** 2 + (e2 / y2) ** 2
    return r, abs(r) * (rel ** 0.5)


def ratio_series_with_error(y1_list, e1_list, y2_list, e2_list):
    """Same as ratio_with_error, applied bin-by-bin to two equal-length
    series (e.g. two graphs sharing the same x values)."""
    if len(y1_list) != len(y2_list):
        raise ValueError("ratio_series_with_error: series must be the same length "
                          "(got %d and %d) -- are these graphs binned the same way?"
                          % (len(y1_list), len(y2_list)))
    ratios, errors = [], []
    for y1, e1, y2, e2 in zip(y1_list, e1_list, y2_list, e2_list):
        r, e = ratio_with_error(y1, e1, y2, e2)
        ratios.append(r)
        errors.append(e)
    return ratios, errors
