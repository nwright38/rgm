"""
plot_src_overlay_hists.py

Overlay plotting for histograms produced by
save_overlay_default_multi_hists.C.

The input ROOT file is expected to contain histograms under:
  <variable>/<selection>/<group>/<topology>/h_<...>_s<sampleIndex>

Examples:
  python plot_src_overlay_hists.py scratch/overlay_src_tools/overlay_default_multi_hists.root

  python plot_src_overlay_hists.py scratch/overlay_src_tools/overlay_default_multi_hists.root \
      --selection epp --topology-group drawn --topology lead_fd_rec_cd --topology lead_cd_rec_cd

  python plot_src_overlay_hists.py scratch/overlay_src_tools/overlay_default_multi_hists.root \
      --output-mode per-plot --out-dir pdf/src_overlay_hists_per_plot
"""

from __future__ import division

import argparse
import math
import os
import re
import tempfile

os.environ.setdefault('MPLCONFIGDIR', os.path.join(tempfile.gettempdir(), 'matplotlib'))
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

import uproot
import ratio


DEFAULT_COLORS = [
    'black', 'red', 'blue', 'green', 'darkorange', 'purple',
    'darkcyan', 'magenta', 'saddlebrown', 'gray',
]
DEFAULT_MARKERS = ['o', 's', '^', 'v', 'D', '<', '>', 'P', 'X', '*']

CORNER_LABELS = {'ep': r"$(e,e^{\prime}p)$", 'epp': r"$(e,e^{\prime}pp)$"}

SAMPLE_RE = re.compile(r'_s(\d+)$')


class OverlayHist(object):
    def __init__(self, variable, selection, group, topology, sample_idx,
                 xlabel, x, xerr_low, xerr_high, y, yerr_low, yerr_high):
        self.variable = variable
        self.selection = selection
        self.group = group
        self.topology = topology
        self.sample_idx = sample_idx
        self.xlabel = xlabel
        self.x = x
        self.xerr_low = xerr_low
        self.xerr_high = xerr_high
        self.y = y
        self.yerr_low = yerr_low
        self.yerr_high = yerr_high


class SampleStyle(object):
    def __init__(self, idx, label, color, marker):
        self.idx = idx
        self.label = label
        self.color = color
        self.marker = marker


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_file',
                   help='ROOT file produced by save_overlay_default_multi_hists.C')
    p.add_argument('--label', action='append', default=[],
                   help='Legend label for sample index order; repeat as needed')
    p.add_argument('--color', action='append', default=[],
                   help='Plot color for sample index order; repeat as needed')
    p.add_argument('--marker', action='append', default=[],
                   help='Plot marker for sample index order; repeat as needed')
    p.add_argument('--selection', choices=['all', 'ep', 'epp'], default='all',
                   help='Selection family to draw')
    p.add_argument('--topology-group', action='append',
                   choices=['drawn', 'lead_only', 'inclusive'], default=[],
                   help='Restrict to topology groups; repeat and/or comma-separate values')
    p.add_argument('--topology', action='append', default=[],
                   help='Restrict to specific topology keys; repeat and/or comma-separate values')
    p.add_argument('--variable', action='append', default=[],
                   help='Restrict to specific variable names; repeat and/or comma-separate values')
    p.add_argument('--output-mode', choices=['book', 'per-plot'], default='book',
                   help='book: one multi-page PDF, per-plot: one PDF per overlay figure')
    p.add_argument('--out-file', default='pdf/src_overlay_hists/overlay_src_hists_book.pdf',
                   help='Output PDF path when --output-mode book')
    p.add_argument('--out-dir', default='pdf/src_overlay_hists',
                   help='Output directory for --output-mode per-plot')
    p.add_argument('--renormalize', action='store_true',
                   help='Normalize each sample histogram to unity before plotting')
    p.add_argument('--no-legend', action='store_true')
    p.add_argument('--ylabel', default='Normalized counts',
                   help='Y axis label')
    p.add_argument('--shift-fraction', type=float, default=0.12,
                   help='Total x-shift span as a fraction of median bin spacing')
    p.add_argument('--run-tag', default='',
                   help='Optional top-level run tag directory inside the ROOT file')
    p.add_argument('--with-ratio', dest='with_ratio', action='store_true',
                   help='Draw a bottom ratio panel when at least two samples are present')
    p.add_argument('--no-ratio', dest='with_ratio', action='store_false',
                   help='Disable the bottom ratio panel')
    p.set_defaults(with_ratio=True)
    p.add_argument('--ratio-reference-index', type=int, default=0,
                   help='Sample index to use as denominator in the ratio panel')
    p.add_argument('--ratio-ylabel', default=None,
                   help='Override the ratio-panel y-axis label')
    p.add_argument('--ratio-ylim', nargs=2, type=float, default=(0.5, 1.5),
                   metavar=('YMIN', 'YMAX'))
    return p.parse_args()


def _strip_cycle(name):
    return name.split(';')[0]


def _sqrt_or_zero(x):
    if x is None or x <= 0.0:
        return 0.0
    return math.sqrt(x)


def _axis_label_or_fallback(hist, variable):
    axis = hist.axes[0]

    # uproot API differs across versions/models; try common accessors first.
    label = getattr(axis, 'label', None)
    if isinstance(label, str) and label.strip():
        return label

    title = getattr(axis, 'title', None)
    if isinstance(title, str) and title.strip():
        return title

    try:
        member_title = axis.member('fTitle')
    except Exception:
        member_title = None
    if isinstance(member_title, str) and member_title.strip():
        return member_title

    return variable


def _root_label_to_mathtext(label):
    if not isinstance(label, str):
        return label
    s = label.strip()
    if not s:
        return s

    # If the label already contains TeX-like backslashes, keep it as-is.
    # Some ROOT titles already use math-like fragments that matplotlib's
    # mathtext parser can misread when we re-wrap them.
    if '\\' in s:
        return s

    root_tokens = {
        '#theta': r'\theta',
        '#phi': r'\phi',
        '#omega': r'\omega',
        '#alpha': r'\alpha',
        '#pi': r'\pi',
        '#vec': r'\vec',
        '#hat': r'\hat',
        '#bar': r'\bar',
        '#cdot': r'\cdot',
    }

    converted = False
    for src, dst in root_tokens.items():
        if src in s:
            s = s.replace(src, dst)
            converted = True

    # Common function names in titles.
    if 'Cos(' in s:
        s = s.replace('Cos(', r'\cos(')
        converted = True
    if 'Sin(' in s:
        s = s.replace('Sin(', r'\sin(')
        converted = True

    if converted:
        return '$%s$' % s
    return s


def _normalize_filter_values(raw_values):
    out = []
    for raw in (raw_values or []):
        for item in str(raw).split(','):
            t = item.strip()
            if t:
                out.append(t)
    return out


def _read_meta_labels(root_file, run_tag=''):
    labels = {}
    try:
        base = root_file[run_tag] if run_tag else root_file
        meta = base['meta']
    except Exception:
        return labels

    for k in meta.keys(cycle=False):
        obj = meta[k]
        if not hasattr(obj, 'title'):
            continue
        name = _strip_cycle(k)
        if not name.startswith('sample_'):
            continue
        try:
            idx = int(name.split('_')[1])
        except Exception:
            continue
        labels[idx] = str(obj.title)
    return labels


def _read_overlay_hists(path, run_tag=''):
    root_file = uproot.open(path)
    class_map = root_file.classnames(recursive=True)

    records = {}
    sample_indices = set()

    for key, class_name in class_map.items():
        if not class_name.startswith('TH1'):
            continue

        clean_parts = [_strip_cycle(p) for p in key.split('/')]
        if run_tag:
            if not clean_parts or clean_parts[0] != run_tag:
                continue
            clean_parts = clean_parts[1:]
        if len(clean_parts) < 5:
            continue
        if clean_parts[0] == 'meta':
            continue

        variable, selection, group, topology, hist_name = clean_parts[:5]
        m = SAMPLE_RE.search(hist_name)
        if not m:
            continue
        sample_idx = int(m.group(1))

        hist = root_file[key]
        values = list(hist.values(flow=False))
        variances = hist.variances(flow=False)
        variances = list(variances) if variances is not None else [0.0] * len(values)
        edges = list(hist.axes[0].edges())

        x = []
        xerr_low = []
        xerr_high = []
        for i in range(len(values)):
            lo = edges[i]
            hi = edges[i + 1]
            center = 0.5 * (lo + hi)
            half = 0.5 * (hi - lo)
            x.append(center)
            xerr_low.append(half)
            xerr_high.append(half)

        y = [float(v) for v in values]
        yerr = [_sqrt_or_zero(vv) for vv in variances]
        xlabel = _axis_label_or_fallback(hist, variable)

        rec = OverlayHist(
            variable=variable,
            selection=selection,
            group=group,
            topology=topology,
            sample_idx=sample_idx,
            xlabel=xlabel,
            x=x,
            xerr_low=xerr_low,
            xerr_high=xerr_high,
            y=y,
            yerr_low=list(yerr),
            yerr_high=list(yerr),
        )

        plot_key = (variable, selection, group, topology)
        records.setdefault(plot_key, {})[sample_idx] = rec
        sample_indices.add(sample_idx)

    return root_file, records, sorted(sample_indices)


def _median_spacing(x):
    spacings = sorted(abs(b - a) for a, b in zip(x[:-1], x[1:]) if abs(b - a) > 0.0)
    if not spacings:
        return 0.0
    mid = len(spacings) // 2
    if len(spacings) % 2:
        return spacings[mid]
    return 0.5 * (spacings[mid - 1] + spacings[mid])


def _offsets(n, spacing, shift_fraction):
    if n <= 1 or spacing == 0.0 or shift_fraction == 0.0:
        return [0.0] * n
    span = spacing * shift_fraction
    if n == 2:
        return [-0.5 * span, 0.5 * span]
    return [span * (i / float(n - 1) - 0.5) for i in range(n)]


def _build_styles(sample_indices, args, meta_labels):
    styles = {}
    for i, sample_idx in enumerate(sample_indices):
        if i < len(args.label):
            label = args.label[i]
        elif sample_idx in meta_labels:
            label = meta_labels[sample_idx]
        else:
            label = 'Sample %d' % (sample_idx + 1)

        color = args.color[i] if i < len(args.color) else DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
        marker = args.marker[i] if i < len(args.marker) else DEFAULT_MARKERS[i % len(DEFAULT_MARKERS)]
        styles[sample_idx] = SampleStyle(sample_idx, label, color, marker)
    return styles


def _normalize(rec):
    integral = sum(rec.y)
    if integral <= 0.0:
        return
    rec.y = [yi / integral for yi in rec.y]
    rec.yerr_low = [ei / integral for ei in rec.yerr_low]
    rec.yerr_high = [ei / integral for ei in rec.yerr_high]


def _passes_filters(plot_key, args):
    variable, selection, group, topology = plot_key
    if args.selection != 'all' and selection != args.selection:
        return False
    if args.topology_group and group not in args.topology_group_set:
        return False
    if args.topology and topology not in args.topology_set:
        return False
    if args.variable and variable not in args.variable_set:
        return False
    return True


def _title_from_key(plot_key):
    variable, selection, group, topology = plot_key
    return '%s | %s | %s | %s' % (variable, selection, group, topology)


def _plot_one(plot_key, sample_map, ordered_sample_indices, styles, args):
    present = [idx for idx in ordered_sample_indices if idx in sample_map]
    if not present:
        raise ValueError('No samples available for key=%s' % (plot_key,))

    first = sample_map[present[0]]
    spacing = _median_spacing(first.x)
    shifts = _offsets(len(present), spacing, args.shift_fraction)

    show_ratio = args.with_ratio and len(present) > 1
    if show_ratio and args.ratio_reference_index not in present:
      raise ValueError('ratio reference sample index %d is not present for key=%s'
                       % (args.ratio_reference_index, plot_key))

    if show_ratio:
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    else:
        fig, ax = plt.subplots(1, 1)
        ax_ratio = None

    ymax = 0.0
    drawn = {}

    x_min = min((x - e) for x, e in zip(first.x, first.xerr_low))
    x_max = max((x + e) for x, e in zip(first.x, first.xerr_high))

    for offset, sample_idx in zip(shifts, present):
        rec = sample_map[sample_idx]
        if args.renormalize:
            _normalize(rec)

        style = styles[sample_idx]
        shifted_x = [x + offset for x in rec.x]
        ax.errorbar(shifted_x, rec.y,
                    xerr=[rec.xerr_low, rec.xerr_high],
                    yerr=[rec.yerr_low, rec.yerr_high],
                    color=style.color, marker=style.marker,
                    linestyle='', markersize=4, label=style.label)
        drawn[sample_idx] = {
            'x': list(rec.x),
            'y': list(rec.y),
            'err': [max(lo, hi) for lo, hi in zip(rec.yerr_low, rec.yerr_high)],
        }
        local_max = max([yy + ee for yy, ee in zip(rec.y, rec.yerr_high)] or [0.0])
        ymax = max(ymax, local_max)

    corner = CORNER_LABELS.get(plot_key[1], plot_key[1])
    ax.text(0.04, 0.92, corner, transform=ax.transAxes, fontsize=16,
            ha='left', va='top')

    ax.set_title(_title_from_key(plot_key))
    ax.set_xlabel(_root_label_to_mathtext(first.xlabel))
    ax.set_ylabel(args.ylabel)
    pad = max([abs(s) for s in shifts] or [0.0])
    ax.set_xlim(x_min - pad, x_max + pad)
    ax.set_ylim(0.0, ymax * 1.2 if ymax > 0.0 else 1.0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    if not args.no_legend:
        ax.legend(frameon=False, fontsize=9)

    if show_ratio:
        ref_idx = args.ratio_reference_index
        ref = drawn[ref_idx]
        ratio_indices = [idx for idx in present if idx != ref_idx]
        ratio_offsets = _offsets(len(ratio_indices), spacing, args.shift_fraction)

        for r_off, sample_idx in zip(ratio_offsets, ratio_indices):
            values = drawn[sample_idx]
            if len(values['x']) != len(ref['x']):
                raise ValueError('Cannot ratio sample %d to %d for key=%s: different bin counts'
                                 % (sample_idx, ref_idx, plot_key))

            ratio_y, ratio_err = ratio.ratio_series_with_error(
                values['y'], values['err'], ref['y'], ref['err'])
            shifted_x = [x + r_off for x in values['x']]

            style = styles[sample_idx]
            ref_label = styles[ref_idx].label
            ax_ratio.errorbar(
                shifted_x, ratio_y, ratio_err,
                color=style.color, marker=style.marker, linestyle='', markersize=4,
                label='%s / %s' % (style.label, ref_label))

        ax_ratio.axhline(1.0, color='gray', linewidth=0.7, linestyle='--')
        ax_ratio.set_ylim(*args.ratio_ylim)
        ax_ratio.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ratio_ylabel = args.ratio_ylabel
        if ratio_ylabel is None:
            ratio_ylabel = 'sample / %s' % styles[ref_idx].label
        ax_ratio.set_ylabel(ratio_ylabel)
        ax_ratio.set_xlabel(_root_label_to_mathtext(first.xlabel))
        if (not args.no_legend) and ratio_indices:
            ax_ratio.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    return fig


def _plot_path(out_dir, plot_key):
    variable, selection, group, topology = plot_key
    name = '%s_%s_%s_%s.pdf' % (variable, selection, group, topology)
    return os.path.join(out_dir, name)


def main():
    args = parse_args()
    args.topology_group = _normalize_filter_values(args.topology_group)
    args.topology = _normalize_filter_values(args.topology)
    args.variable = _normalize_filter_values(args.variable)
    args.topology_group_set = set(args.topology_group)
    args.topology_set = set(args.topology)
    args.variable_set = set(args.variable)

    root_file, records, sample_indices = _read_overlay_hists(args.input_file, args.run_tag)
    if not records:
        if args.run_tag:
            raise ValueError('No overlay histograms found in %s for runTag=%s'
                             % (args.input_file, args.run_tag))
        raise ValueError('No overlay histograms found in %s' % args.input_file)

    meta_labels = _read_meta_labels(root_file, args.run_tag)
    styles = _build_styles(sample_indices, args, meta_labels)

    keys = [k for k in sorted(records.keys()) if _passes_filters(k, args)]
    if not keys:
        raise ValueError('No histograms matched the requested filters')

    if args.output_mode == 'book':
        out_parent = os.path.dirname(args.out_file)
        if out_parent:
            os.makedirs(out_parent, exist_ok=True)
        with PdfPages(args.out_file) as pdf:
            for key in keys:
                fig = _plot_one(key, records[key], sample_indices, styles, args)
                pdf.savefig(fig)
                plt.close(fig)
                print('Added page for %s' % (_title_from_key(key),))
        print('Wrote %s' % args.out_file)
    else:
        os.makedirs(args.out_dir, exist_ok=True)
        for key in keys:
            fig = _plot_one(key, records[key], sample_indices, styles, args)
            out_path = _plot_path(args.out_dir, key)
            fig.savefig(out_path)
            plt.close(fig)
            print('Wrote %s' % out_path)


if __name__ == '__main__':
    main()
