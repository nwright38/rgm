"""
plot_1D_overlay.py

Unit-normalized point overlays for the same 1D distributions used by
plot_1D.py. Unlike plot_1D.py, every input file is treated generically: data,
simulation, alternate data releases, or systematic variants can all be
overlaid together as shifted points with error bars.

Examples:
    python plot_1D_overlay.py Data_He.root Sim_AV18.root Sim_N2LO.root \\
        --label Data --label AV18 --label N2LO --with-ratio

    python plot_1D_overlay.py Data_new.root Data_old.root \\
        --label "new data" --label "reference data" \\
        --with-ratio --ratio-reference-index 1

By default the ratio panel shows every non-reference file divided by the
reference input file (input/reference). Use --ratio-mode reference-over-inputs
for reference/input, or --ratio-mode pairwise for fully manual
numerator/denominator choices.
"""

from __future__ import division

import argparse
import os
import tempfile

os.environ.setdefault('MPLCONFIGDIR', os.path.join(tempfile.gettempdir(), 'matplotlib'))
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import graph_io
import graph_names
import ratio


# (output key, task_name, selection, xlabel, xlim)
VARIABLES = [
    ('pMiss_ep', 'pMiss_ep_note', 'ep', r'$p_{miss} [GeV]$', (0.4, 1.)),
    ('pMiss_epp', 'pMiss_epp_note', 'epp', r'$p_{miss} [GeV]$', (0.4, 1.)),
    ('kMiss_ep', 'kMiss_ep_note', 'ep', r'$k_{miss} [GeV]$', (0.3, 1.)),
    ('kMiss_epp', 'kMiss_epp_note', 'epp', r'$k_{miss} [GeV]$', (0.3, 1.)),
    ('pRel_epp', 'pRel_epp', 'epp', r'$p_{rel} [GeV]$', (0.15, 1.0)),
    ('pcmz_epp', 'pcmz_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{z} [GeV]$', (-0.75, 0.75)),
    ('pcmx_epp', 'pcmx_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{x} [GeV]$', (-0.75, 0.75)),
    ('pcmy_epp', 'pcmy_epp', 'epp', r'$\vec{p}_{C.M.} \cdot \hat{v}_{y} [GeV]$', (-0.75, 0.75)),
    ('theta_pmiss_ep', 'theta_pmiss_ep', 'ep', r'$\theta_{p_{miss},q} [deg]$', (100, 180)),
    ('theta_pmiss_epp', 'theta_pmiss_epp', 'epp', r'$\theta_{p_{miss},q} [deg]$', (100, 180)),
    ('q_ep', 'q_ep', 'ep', r'$|\vec{q}| [GeV]$', (0, 4)),
    ('q_epp', 'q_epp', 'epp', r'$|\vec{q}| [GeV]$', (0, 4)),
    ('theta_pLeadq_ep', 'theta_pLeadq_ep', 'ep', r'$\theta_{p_{lead},q} [deg]$', (0, 40)),
    ('theta_pLeadq_epp', 'theta_pLeadq_epp', 'epp', r'$\theta_{p_{lead},q} [deg]$', (0, 40)),
]

DEFAULT_COLORS = [
    'black', 'red', 'blue', 'green', 'darkorange', 'purple',
    'darkcyan', 'magenta', 'saddlebrown', 'gray',
]
DEFAULT_MARKERS = ['o', 's', '^', 'v', 'D', '<', '>', 'P', 'X', '*']
CORNER_LABELS = {'ep': r"$(e,e^{\prime}p)$", 'epp': r"$(e,e^{\prime}pp)$"}


class OverlayInput(object):
    def __init__(self, path, label, color, marker, root_file):
        self.path = path
        self.label = label
        self.color = color
        self.marker = marker
        self.file = root_file


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_files', nargs='+',
                   help='ROOT files after Main_Figs_Binned + BuildGraphs')
    p.add_argument('--label', action='append', default=[],
                   help='Legend label for an input file; repeat in input-file order')
    p.add_argument('--color', action='append', default=[],
                   help='Matplotlib color for an input file; repeat in input-file order')
    p.add_argument('--marker', action='append', default=[],
                   help='Matplotlib marker for an input file; repeat in input-file order')
    p.add_argument('--out-dir', default='pdf/1D_overlay')
    p.add_argument('--out-format', default='pdf', choices=['pdf', 'png'])
    p.add_argument('--variable', action='append', default=[],
                   help='Only plot this output key from VARIABLES; repeat as needed')
    p.add_argument('--with-ratio', action='store_true',
                   help='Add a bottom ratio panel')
    p.add_argument('--ratio-mode', default='inputs-over-reference',
                   choices=['inputs-over-reference', 'reference-over-inputs', 'pairwise'],
                   help='Ratio panel convention')
    p.add_argument('--ratio-reference-index', type=int, default=0,
                   help='Reference input-file index for inputs-over-reference or reference-over-inputs')
    p.add_argument('--ratio-numerator-index', type=int, default=0,
                   help='Numerator input-file index for pairwise ratio mode')
    p.add_argument('--ratio-denominator-index', type=int, action='append', default=[],
                   help='Denominator/reference input-file index; repeat for pairwise overlays')
    p.add_argument('--ratio-ylabel', default=None,
                   help='Override the bottom-panel y-axis label')
    p.add_argument('--ratio-ylim', nargs=2, type=float, default=(0.5, 1.5),
                   metavar=('YMIN', 'YMAX'))
    p.add_argument('--shift-fraction', type=float, default=0.12,
                   help='Total x-shift span as a fraction of median bin spacing')
    p.add_argument('--no-legend', action='store_true')
    return p.parse_args()


def _check_index(index, n_inputs, option_name):
    if index < 0 or index >= n_inputs:
        raise ValueError('%s=%d is out of range for %d input files'
                         % (option_name, index, n_inputs))


def _build_inputs(args):
    inputs = []
    for i, path in enumerate(args.input_files):
        label = args.label[i] if i < len(args.label) else os.path.basename(path)
        color = args.color[i] if i < len(args.color) else DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
        marker = args.marker[i] if i < len(args.marker) else DEFAULT_MARKERS[i % len(DEFAULT_MARKERS)]
        inputs.append(OverlayInput(path, label, color, marker, graph_io.open_file(path)))
    return inputs


def _selected_variables(args):
    if not args.variable:
        return VARIABLES
    wanted = set(args.variable)
    variables = [v for v in VARIABLES if v[0] in wanted]
    missing = sorted(wanted - set(v[0] for v in variables))
    if missing:
        raise ValueError('Unknown --variable key(s): %s' % ', '.join(missing))
    return variables


def _read_normalized(inp, task_name, selection):
    name = graph_names.diff_graph_name(task_name, selection, [])
    x, y, yerr_low, yerr_high = graph_io.read_graph_asymm(inp.file, name)
    integral = sum(y)
    if integral == 0.0:
        raise ValueError('Integral is zero for graph %s in %s' % (name, inp.path))
    y_norm = [yi / integral for yi in y]
    err_low_norm = [ei / integral for ei in yerr_low]
    err_high_norm = [ei / integral for ei in yerr_high]
    err_sym_norm = [max(lo, hi) for lo, hi in zip(err_low_norm, err_high_norm)]
    return x, y_norm, err_low_norm, err_high_norm, err_sym_norm


def _median_spacing(x):
    spacings = sorted(abs(b - a) for a, b in zip(x[:-1], x[1:]) if abs(b - a) > 0.0)
    if not spacings:
        return 0.0
    mid = len(spacings) // 2
    if len(spacings) % 2:
        return spacings[mid]
    return 0.5 * (spacings[mid - 1] + spacings[mid])


def _offsets(n_inputs, spacing, shift_fraction):
    if n_inputs <= 1 or spacing == 0.0 or shift_fraction == 0.0:
        return [0.0] * n_inputs
    span = spacing * shift_fraction
    if n_inputs == 2:
        return [-0.5 * span, 0.5 * span]
    return [span * (i / float(n_inputs - 1) - 0.5) for i in range(n_inputs)]


def _draw_points(ax, x, y, yerr_low, yerr_high, offset, inp, markersize=4):
    shifted_x = [xi + offset for xi in x]
    ax.errorbar(shifted_x, y, yerr=[yerr_low, yerr_high],
                color=inp.color, marker=inp.marker, linestyle='',
                markersize=markersize, label=inp.label)
    return shifted_x


def _manual_ratio_denominators(args, n_inputs):
    _check_index(args.ratio_numerator_index, n_inputs, '--ratio-numerator-index')
    if args.ratio_denominator_index:
        denominators = list(args.ratio_denominator_index)
    else:
        denominators = [i for i in range(n_inputs) if i != args.ratio_numerator_index]
    for idx in denominators:
        _check_index(idx, n_inputs, '--ratio-denominator-index')
    if args.ratio_numerator_index in denominators:
        raise ValueError('Ratio numerator cannot also be a denominator')
    return denominators


def _ratio_specs(args, n_inputs):
    if args.ratio_mode == 'pairwise':
        return [(args.ratio_numerator_index, den_idx)
                for den_idx in _manual_ratio_denominators(args, n_inputs)]

    _check_index(args.ratio_reference_index, n_inputs, '--ratio-reference-index')
    others = [i for i in range(n_inputs) if i != args.ratio_reference_index]
    if args.ratio_denominator_index:
        others = list(args.ratio_denominator_index)
        for idx in others:
            _check_index(idx, n_inputs, '--ratio-denominator-index')
        if args.ratio_reference_index in others:
            raise ValueError('Reference input cannot also be listed as a ratio partner')

    if args.ratio_mode == 'inputs-over-reference':
        return [(idx, args.ratio_reference_index) for idx in others]
    return [(args.ratio_reference_index, idx) for idx in others]


def _plot_one(inputs, variable, args, out_path):
    key, task_name, selection, xlabel, xlim = variable
    normalized = [_read_normalized(inp, task_name, selection) for inp in inputs]
    spacing = _median_spacing(normalized[0][0])
    offsets = _offsets(len(inputs), spacing, args.shift_fraction)

    show_ratio = args.with_ratio and len(inputs) > 1
    if show_ratio:
        ratio_specs = _ratio_specs(args, len(inputs))
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    else:
        ratio_specs = []
        fig, ax = plt.subplots(1, 1)
        ax_ratio = None

    ymax = 0.0
    for inp, offset, values in zip(inputs, offsets, normalized):
        x, y, yerr_low, yerr_high, _ = values
        _draw_points(ax, x, y, yerr_low, yerr_high, offset, inp)
        local_max = max([yi + ei for yi, ei in zip(y, yerr_high)] or [0.0])
        ymax = max(ymax, local_max)

    ax.set_ylabel('Unit-normalized counts')
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, ymax * 1.2 if ymax > 0.0 else 1.0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.text(0.04, 0.92, CORNER_LABELS.get(selection, ''),
            transform=ax.transAxes, fontsize=18, ha='left', va='top')
    if not args.no_legend:
        ax.legend(frameon=False, fontsize=9)

    if show_ratio:
        ratio_offsets = _offsets(len(ratio_specs), spacing, args.shift_fraction)
        for j, (num_idx, den_idx) in enumerate(ratio_specs):
            num_x, num_y, _, _, num_err = normalized[num_idx]
            den_x, den_y, _, _, den_err = normalized[den_idx]
            if len(num_x) != len(den_x):
                raise ValueError('Cannot ratio %s and %s for %s: different point counts'
                                 % (inputs[num_idx].label, inputs[den_idx].label, key))
            r, e = ratio.ratio_series_with_error(num_y, num_err, den_y, den_err)
            shifted_x = [xi + ratio_offsets[j] for xi in num_x]
            label = '%s / %s' % (inputs[num_idx].label, inputs[den_idx].label)
            ax_ratio.errorbar(shifted_x, r, e, color=inputs[num_idx].color,
                              marker=inputs[num_idx].marker, linestyle='',
                              markersize=4, label=label)
        ax_ratio.axhline(1.0, color='gray', linewidth=0.7, linestyle='--')
        ylabel = args.ratio_ylabel
        if ylabel is None:
            if args.ratio_mode == 'inputs-over-reference':
                ylabel = 'input / %s' % inputs[args.ratio_reference_index].label
            elif args.ratio_mode == 'reference-over-inputs':
                ylabel = '%s / input' % inputs[args.ratio_reference_index].label
            else:
                ylabel = 'ratio'
        ax_ratio.set_ylabel(ylabel)
        ax_ratio.set_xlabel(xlabel)
        ax_ratio.set_ylim(*args.ratio_ylim)
        ax_ratio.yaxis.set_major_locator(MaxNLocator(nbins=4))
        if not args.no_legend:
            ax_ratio.legend(frameon=False, fontsize=8)
    else:
        ax.set_xlabel(xlabel)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    inputs = _build_inputs(args)
    variables = _selected_variables(args)

    missing = []
    for variable in variables:
        key = variable[0]
        out_path = os.path.join(args.out_dir, key + '.' + args.out_format)
        try:
            _plot_one(inputs, variable, args, out_path)
        except Exception as exc:
            missing.append('%s: %s' % (key, exc))
            continue
        print('Wrote %s' % out_path)

    if missing:
        print('Skipped variables with missing or incompatible graphs:')
        for msg in missing:
            print(' - %s' % msg)


if __name__ == '__main__':
    main()
