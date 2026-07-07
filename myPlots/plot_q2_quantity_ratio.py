"""
plot_q2_quantity_ratio.py

Ratio plots for graph quantities between two Main_Figs_Binned output ROOT
files. Defaults to sigma_cm_x vs Q2, sigma_cm_y vs Q2, and the 1D
theta_pLead,q distributions. The same machinery can plot any graph stored
with the BuildGraphs/ExtractFitQuantities name convention:

    "<task_name>|<selection>|<axis-bin-suffix>"

Examples:
    python plot_q2_quantity_ratio.py numerator.root denominator.root
    python plot_q2_quantity_ratio.py stat.root nominal.root --quantity sigma_pcmx
    python plot_q2_quantity_ratio.py a.root b.root \\
        --quantity mean_E1miss_ep_pmiss:ep:0:'Mean E1miss, pMiss bin 0'
"""

from __future__ import division

import argparse
from collections import OrderedDict, namedtuple
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

import graph_io
import graph_names
import ratio


QuantitySpec = namedtuple(
    'QuantitySpec',
    ['task_name', 'selection', 'axis_bin', 'label', 'xlabel', 'xlim', 'ylim'])


DEFAULT_QUANTITIES = OrderedDict([
    ('sigma_pcmx', QuantitySpec(
        'sigma_pcmx', 'epp', (), r'$\sigma_{x,C.M.}$',
        r'$Q^{2}$', (1.5, 5.0), None)),
    ('sigma_pcmy', QuantitySpec(
        'sigma_pcmy', 'epp', (), r'$\sigma_{y,C.M.}$',
        r'$Q^{2}$', (1.5, 5.0), None)),
    ('theta_pLeadq_ep', QuantitySpec(
        'theta_pLeadq_ep', 'ep', (), r'$\theta_{p_{lead},q}$ $(e,e^{\prime}p)$',
        r'$\theta_{p_{lead},q} [deg]$', (0.0, 40.0), None)),
    ('theta_pLeadq_epp', QuantitySpec(
        'theta_pLeadq_epp', 'epp', (), r'$\theta_{p_{lead},q}$ $(e,e^{\prime}pp)$',
        r'$\theta_{p_{lead},q} [deg]$', (0.0, 40.0), None)),
])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('numerator_file')
    p.add_argument('denominator_file')
    p.add_argument('--numerator-label', default='Numerator')
    p.add_argument('--denominator-label', default='Denominator')
    p.add_argument('--quantity', action='append', default=[],
                   help=("Quantity to plot. Use a predefined key like sigma_pcmx "
                         "or custom TASK[:SELECTION[:AXIS_BINS[:LABEL]]]. "
                         "AXIS_BINS is comma-separated, e.g. 2 or 1,3."))
    p.add_argument('--out-file', default='pdf/quantity_ratios.pdf',
                   help='Multi-page PDF output path.')
    p.add_argument('--xlim', nargs=2, type=float, default=None,
                   help='Override x limits for every page.')
    p.add_argument('--ylim', nargs=2, type=float, default=None)
    p.add_argument('--marker', default='o')
    p.add_argument('--color', default='black')
    return p.parse_args()


def _parse_axis_bin(text):
    if text in ('', 'none', 'None'):
        return ()
    return tuple(int(x) for x in text.split(',') if x != '')


def quantity_from_token(token):
    if token in DEFAULT_QUANTITIES:
        return DEFAULT_QUANTITIES[token]

    parts = token.split(':', 3)
    task_name = parts[0]
    selection = parts[1] if len(parts) > 1 and parts[1] else 'epp'
    axis_bin = _parse_axis_bin(parts[2]) if len(parts) > 2 else ()
    label = parts[3] if len(parts) > 3 and parts[3] else task_name
    return QuantitySpec(
        task_name, selection, axis_bin, label,
        r'$Q^{2}$', (1.5, 5.0), None)


def requested_quantities(args):
    tokens = args.quantity or list(DEFAULT_QUANTITIES.keys())
    return [quantity_from_token(t) for t in tokens]


def read_quantity(root_file, spec):
    graph_name = graph_names.diff_graph_name(
        spec.task_name, spec.selection, list(spec.axis_bin))
    return graph_io.read_graph(root_file, graph_name)


def check_matching_x(x_num, x_den, spec, tolerance=1e-6):
    if len(x_num) != len(x_den):
        raise ValueError(
            '%s has different point counts in numerator and denominator (%d vs %d)'
            % (spec.task_name, len(x_num), len(x_den)))
    mismatches = [
        (a, b) for a, b in zip(x_num, x_den)
        if abs(a - b) > tolerance
    ]
    if mismatches:
        first = mismatches[0]
        print('Warning: %s x values differ; plotting ratio at numerator x. '
              'First mismatch: %.6g vs %.6g'
              % (spec.task_name, first[0], first[1]))


def ratio_for_quantity(num_file, den_file, spec):
    x_num, y_num, err_num = read_quantity(num_file, spec)
    x_den, y_den, err_den = read_quantity(den_file, spec)
    check_matching_x(x_num, x_den, spec)
    y_ratio, err_ratio = ratio.ratio_series_with_error(
        y_num, err_num, y_den, err_den)
    return x_num, y_ratio, err_ratio


def plot_quantity_ratio(pdf, x, y, yerr, spec, args):
    fig, ax = plt.subplots(1, 1)
    ax.errorbar(x, y, yerr=yerr, color=args.color,
                linestyle='', marker=args.marker)
    ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')
    xlim = tuple(args.xlim) if args.xlim is not None else spec.xlim
    if xlim is not None:
        ax.set_xlim(*xlim)
    if args.ylim is not None:
        ax.set_ylim(*args.ylim)
    elif spec.ylim is not None:
        ax.set_ylim(*spec.ylim)
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel('%s / %s' % (args.numerator_label, args.denominator_label))
    ax.set_title(spec.label)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    args = parse_args()
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    num_file = graph_io.open_file(args.numerator_file)
    den_file = graph_io.open_file(args.denominator_file)

    missing = []
    n_written = 0
    with PdfPages(args.out_file) as pdf:
        for spec in requested_quantities(args):
            try:
                x, y, yerr = ratio_for_quantity(num_file, den_file, spec)
                plot_quantity_ratio(pdf, x, y, yerr, spec, args)
                n_written += 1
            except Exception as exc:
                missing.append('Skipped %s: %s' % (spec.task_name, exc))

    print('Wrote %d page(s) to %s' % (n_written, args.out_file))

    if missing:
        print('Skipped quantities:')
        for msg in missing:
            print('  ' + msg)


if __name__ == '__main__':
    main()
