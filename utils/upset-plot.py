#!/usr/bin/env python3
"""
Generate an UpSet plot (Lex et al. 2014) showing the exact intersection
structure of the source sets contributing genes to a class (curated or
candidate), replacing the 3-set-capped Venn diagram in utils/venn-plot.py.

Unlike a Venn diagram, an UpSet plot is not limited to 3 sets and shows
every observed intersection combination explicitly -- e.g. it can show
the Orphanet neurological/developmental AND relationship that defines
curated eligibility directly, instead of pre-collapsing them into a
single "orpha" union set.

Input format is identical to utils/venn-plot.py: a 4-column TSV
(gene, class, score, pipe-separated subgroups), e.g. data/gene_all_score.txt.
"""

import sys
import os
from collections import defaultdict


def _apply_upsetplot_pandas3_patch(_up_plotting, pd, np):
    """
    upsetplot 0.9.0 (the latest release as of this writing) has two rendering
    bugs surfaced by pandas>=3 (mandatory copy-on-write, which turns its
    `.fillna(..., inplace=True)` calls into no-ops) and by newer numpy (which
    no longer silently converts a 1-element ndarray to a Python scalar where
    matplotlib's Text artist expects one). Both are patched in-process here
    rather than pinning an older/incompatible pandas or numpy in the shared
    conda environment. Safe to drop once upstream upsetplot ships a fix.
    """
    from upsetplot.plotting import util as _util

    def _patched_plot_matrix(self, ax):
        ax = self._reorient(ax)
        data = self.intersections
        n_cats = data.index.nlevels
        inclusion = data.index.to_frame().values
        styles = [
            [
                self.subset_styles[i]
                if inclusion[i, j]
                else {"facecolor": self._other_dots_color, "linewidth": 0}
                for j in range(n_cats)
            ]
            for i in range(len(data))
        ]
        styles = sum(styles, [])
        style_columns = {
            "facecolor": "facecolors",
            "edgecolor": "edgecolors",
            "linewidth": "linewidths",
            "linestyle": "linestyles",
            "hatch": "hatch",
        }
        styles = (
            pd.DataFrame(styles)
            .reindex(columns=style_columns.keys())
            .astype(
                {
                    "facecolor": "O",
                    "edgecolor": "O",
                    "linewidth": float,
                    "linestyle": "O",
                    "hatch": "O",
                }
            )
        )
        # non-inplace .fillna assignment -- the pandas-3-safe form
        styles["linewidth"] = styles["linewidth"].fillna(1)
        styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
        styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
        styles["linestyle"] = styles["linestyle"].fillna("solid")
        del styles["hatch"]

        x = np.repeat(np.arange(len(data)), n_cats)
        y = np.tile(np.arange(n_cats), len(data))
        s = (self._element_size * 0.35) ** 2 if self._element_size is not None else 200
        ax.scatter(*self._swapaxes(x, y), s=s, zorder=10, **styles.rename(columns=style_columns))

        if self._with_lines:
            idx = np.flatnonzero(inclusion)
            line_data = pd.Series(y[idx], index=x[idx]).groupby(level=0).aggregate(["min", "max"])
            colors = pd.Series(
                [style.get("edgecolor", style.get("facecolor", self._facecolor)) for style in self.subset_styles],
                name="color",
            )
            line_data = line_data.join(colors)
            ax.vlines(line_data.index.values, line_data["min"], line_data["max"], lw=2,
                      colors=line_data["color"], zorder=5)

        tick_axis = ax.yaxis
        tick_axis.set_ticks(np.arange(n_cats))
        tick_axis.set_ticklabels(data.index.names, rotation=0 if self._horizontal else -90)
        ax.xaxis.set_visible(False)
        ax.tick_params(axis="both", which="both", length=0)
        if not self._horizontal:
            ax.yaxis.set_ticks_position("top")
        ax.set_frame_on(False)
        ax.set_xlim(-0.5, x[-1] + 0.5, auto=False)
        ax.grid(False)

    def _patched_label_sizes(self, ax, rects, where):
        if not self._show_counts and not self._show_percentages:
            return
        if self._show_counts is True:
            count_fmt = "{:.0f}"
        else:
            count_fmt = self._show_counts
            if "{" not in count_fmt:
                count_fmt = _util.to_new_pos_format(count_fmt)
        pct_fmt = "{:.1%}" if self._show_percentages is True else self._show_percentages

        if count_fmt and pct_fmt:
            fmt = f"{count_fmt}\n({pct_fmt})" if where == "top" else f"{count_fmt} ({pct_fmt})"
            def make_args(val):
                return val, val / self.total
        elif count_fmt:
            fmt = count_fmt
            def make_args(val):
                return (val,)
        else:
            fmt = pct_fmt
            def make_args(val):
                return (val / self.total,)

        # float(...) -- collapses the 1-element ndarray np.diff() returns to a
        # Python scalar, which matplotlib's Text artist requires
        if where == "right":
            margin = float(0.01 * abs(np.diff(ax.get_xlim()))[0])
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(width + margin, rect.get_y() + rect.get_height() * 0.5,
                        fmt.format(*make_args(width)), ha="left", va="center")
        elif where == "left":
            margin = float(0.01 * abs(np.diff(ax.get_xlim()))[0])
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(width + margin, rect.get_y() + rect.get_height() * 0.5,
                        fmt.format(*make_args(width)), ha="right", va="center")
        elif where == "top":
            margin = float(0.01 * abs(np.diff(ax.get_ylim()))[0])
            for rect in rects:
                height = rect.get_height() + rect.get_y()
                ax.text(rect.get_x() + rect.get_width() * 0.5, height + margin,
                        fmt.format(*make_args(height)), ha="center", va="bottom")
        else:
            raise NotImplementedError("unhandled where: %r" % where)

    _up_plotting.UpSet.plot_matrix = _patched_plot_matrix
    _up_plotting.UpSet._label_sizes = _patched_label_sizes


def main():
    if len(sys.argv) < 3:
        print("Usage: python upset-plot.py <input.tsv> <class_name> [subgroups] [output_file]")
        print()
        print("Arguments:")
        print("  input.tsv    - File with 4 columns: gene, class, score, pipe-separated subgroups")
        print("  class_name   - Class to filter (e.g., curated, candidate, no_evidence)")
        print("  subgroups    - Comma-separated subgroups to compare (e.g., sanchis,hc,neuro,develop,sfari+)")
        print("                 Note: 'orpha' is a special subgroup = union of 'neuro' and 'develop'")
        print("                 If not provided, uses all subgroups found in the class")
        print("  output_file  - Output PNG path (default: upset_<class_name>.png)")
        print()
        print("Note: Lines starting with # are treated as comments and skipped")
        sys.exit(1)

    input_file = sys.argv[1]
    class_name = sys.argv[2]
    subgroups_arg = sys.argv[3] if len(sys.argv) > 3 else None
    output_file = sys.argv[4] if len(sys.argv) > 4 else f"upset_{class_name}.png"

    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    specified_subgroups = []
    if subgroups_arg:
        specified_subgroups = [s.strip().lower() for s in subgroups_arg.split(',') if s.strip()]

    # Read gene list and organize by subgroup (same parsing as venn-plot.py)
    gene_subgroups = defaultdict(set)  # subgroup -> set of genes
    class_genes = set()

    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 4:
                print(f"Warning: Line {line_num} has fewer than 4 columns, skipping", file=sys.stderr)
                continue

            gene, gene_class, score, evidence = (f.strip() for f in fields[:4])

            if not gene or gene_class != class_name:
                continue

            class_genes.add(gene)

            if evidence.lower() != 'none' and evidence:
                for subgroup in (s.strip().lower() for s in evidence.split('|') if s.strip()):
                    gene_subgroups[subgroup].add(gene)

    if not class_genes:
        print(f"Error: No genes found for class '{class_name}'", file=sys.stderr)
        sys.exit(1)

    # "orpha" special subgroup = union of 'neuro' and 'develop' (kept for parity with venn-plot.py;
    # by default we keep neuro/develop separate instead, since that's exactly the intersection
    # detail an UpSet plot -- unlike a Venn -- can show without collapsing sets)
    if 'neuro' in gene_subgroups or 'develop' in gene_subgroups:
        orpha_genes = set()
        orpha_genes.update(gene_subgroups.get('neuro', set()))
        orpha_genes.update(gene_subgroups.get('develop', set()))
        if orpha_genes:
            gene_subgroups['orpha'] = orpha_genes

    if specified_subgroups:
        subgroups_to_plot = []
        missing = []
        for subgroup in specified_subgroups:
            if subgroup in gene_subgroups:
                subgroups_to_plot.append(subgroup)
            else:
                missing.append(subgroup)
        if missing:
            print(f"Warning: Subgroups not found: {', '.join(missing)}", file=sys.stderr)
    else:
        subgroups_to_plot = sorted(gene_subgroups.keys())

    if len(subgroups_to_plot) < 2:
        print("Error: Need at least 2 subgroups to create an UpSet plot", file=sys.stderr)
        print(f"Available subgroups: {', '.join(sorted(gene_subgroups.keys()))}", file=sys.stderr)
        sys.exit(1)

    try:
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from upsetplot import UpSet, from_contents
        import upsetplot.plotting as _up_plotting
    except ImportError as e:
        print(f"Error: Required libraries not found: {e}", file=sys.stderr)
        print("Please install: pip install matplotlib upsetplot", file=sys.stderr)
        sys.exit(1)

    _apply_upsetplot_pandas3_patch(_up_plotting, pd, np)

    contents = {sg: gene_subgroups[sg] for sg in subgroups_to_plot}
    membership = from_contents(contents)

    fig = plt.figure(figsize=(max(8, 1.5 * len(subgroups_to_plot) + 4.5), 6))
    upset = UpSet(
        membership,
        subset_size='count',
        show_counts=True,
        sort_by='cardinality',
        # Drop the left-hand per-category totals bar chart entirely -- with it removed,
        # the intersection-size bar chart (top) and the membership matrix (bottom) share
        # the same column grid and therefore render at exactly the same width. The set
        # totals it would have shown are added back as right-hand row labels below instead.
        totals_plot_elements=0,
    )
    axes = upset.plot(fig=fig)
    # Style the intersection-size bars to match the site's accent color (avoids passing
    # facecolor to UpSet() directly, which triggers a pandas-version styling bug in
    # upsetplot 0.9.0's fillna(..., inplace=True) chained-assignment path).
    if "intersections" in axes:
        inter_ax = axes["intersections"]
        for bar in inter_ax.patches:
            bar.set_facecolor('#4689a3')
        # upsetplot only draws the left spine by default; add the remaining three so the
        # intersection-size bar chart reads as a fully bordered plot, matching its width
        # to the membership matrix below it.
        for spine in inter_ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#333333')
            spine.set_linewidth(1.0)

    # Total set size per category, printed to the right of each matrix row (in place of
    # the removed left-hand totals bar chart).
    if "matrix" in axes:
        mat_ax = axes["matrix"]
        xlim = mat_ax.get_xlim()
        margin = 0.03 * abs(xlim[1] - xlim[0])
        row_labels = [t.get_text() for t in mat_ax.get_yticklabels()]
        for row_idx, category in enumerate(row_labels):
            total = len(gene_subgroups.get(category, set()))
            mat_ax.text(
                xlim[1] + margin, row_idx, f"{total:,}",
                ha="left", va="center", fontsize=9, color="#333333",
            )
        # keep the data xlim unchanged so the added labels sit in the margin outside it
        mat_ax.set_xlim(xlim)

    fig.suptitle(f"UpSet plot \u2014 {class_name.capitalize()} genes "
                 f"(total: {len(class_genes):,})", fontsize=13, fontweight='bold')

    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"UpSet plot saved to: {output_file}")

    # Print statistics, mirroring venn-plot.py's reporting style
    total_genes = len(class_genes)
    print(f"\nStatistics for {class_name} genes:")
    print(f"Total genes in class: {total_genes:,}")
    print()
    for subgroup in subgroups_to_plot:
        genes = gene_subgroups[subgroup]
        pct = (len(genes) / total_genes * 100) if total_genes > 0 else 0
        if subgroup == 'orpha':
            print(f"  {subgroup} (union of neuro + develop): {len(genes):,} genes ({pct:.1f}%)")
        else:
            print(f"  {subgroup}: {len(genes):,} genes ({pct:.1f}%)")

    # Genes contributing to the class via exactly one source (i.e. would be lost if that
    # source were dropped) -- the per-source unique-contribution figure referenced in the
    # reviewer response (R2.3).
    print("\nUnique-contribution genes (belong to exactly one listed source):")
    for subgroup in subgroups_to_plot:
        others = set()
        for other in subgroups_to_plot:
            if other != subgroup:
                others |= gene_subgroups[other]
        unique = gene_subgroups[subgroup] - others
        print(f"  {subgroup} only: {len(unique):,} genes")

    covered = set()
    for sg in subgroups_to_plot:
        covered |= gene_subgroups[sg]
    uncovered = class_genes - covered
    if uncovered:
        print(f"\nWarning: {len(uncovered):,} {class_name} genes are not covered by any of the "
              f"selected subgroups ({', '.join(subgroups_to_plot)}) -- likely captured by a "
              f"source flag not included in this plot.", file=sys.stderr)


if __name__ == "__main__":
    main()
