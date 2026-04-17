#!/usr/bin/env python3
"""
Create a stacked horizontal bar chart from multiple TSV taxonomy tables.

Features
--------
- Reads all .tsv files from a folder
- Ignores comment lines starting with '#'
- Ignores lines starting with 'Unassigned'
- Aggregates repeated entries of the same selected rank within each file
- Normalizes each file to percentages
- Automatically creates an "Other" category:
    all groups with less than X percent within a given bar are merged into "Other"
- Uses a fixed global category order across all bars
- Uses a fixed color per category across all bars
- Optional filtering:
    only include rows whose taxonomy contains a particular name at a particular rank
    using --filter-rank and --filter-name together
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict

import pandas as pd
import matplotlib.pyplot as plt


VALID_RANKS = ("d", "p", "c", "o", "f", "g")
RANK_LABELS = {
    "d": "Domain",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
}


def parse_taxonomy_field(taxonomy: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}

    for part in taxonomy.strip().split(";"):
        part = part.strip()
        if "__" not in part:
            continue

        rank, value = part.split("__", 1)
        rank = rank.strip()
        value = value.strip()

        if rank in VALID_RANKS:
            parsed[rank] = value

    return parsed


def normalize_label(label: str) -> str:
    if not label or label == "__":
        return "Unclassified"
    return label


def passes_filter(
    taxonomy_map: Dict[str, str],
    filter_rank: str | None,
    filter_name: str | None,
) -> bool:
    """
    Return True if the row passes the optional taxonomy filter.

    If no filter is specified, all rows pass.
    If both filter_rank and filter_name are specified, only rows where
    taxonomy_map[filter_rank] == filter_name pass.
    """
    if filter_rank is None and filter_name is None:
        return True

    if filter_rank is None or filter_name is None:
        return False

    return taxonomy_map.get(filter_rank, "") == filter_name


def read_and_aggregate_tsv(
    file_path: Path,
    rank: str,
    filter_rank: str | None = None,
    filter_name: str | None = None,
) -> Dict[str, float]:
    aggregated: Dict[str, float] = defaultdict(float)

    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("Unassigned"):
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"{file_path}: line {line_number} malformed (expected at least 2 columns)"
                )

            taxonomy = parts[0].strip()
            value_str = parts[1].strip()

            try:
                value = float(value_str)
            except ValueError as exc:
                raise ValueError(
                    f"{file_path}: line {line_number} has invalid numeric value {value_str!r}"
                ) from exc

            taxonomy_map = parse_taxonomy_field(taxonomy)

            if not passes_filter(taxonomy_map, filter_rank, filter_name):
                continue

            label = normalize_label(taxonomy_map.get(rank, "Unclassified"))
            aggregated[label] += value

    return dict(aggregated)


def collapse_small_categories(
    percentages: Dict[str, float],
    threshold: float,
) -> Dict[str, float]:
    if threshold <= 0:
        return dict(percentages)

    collapsed: Dict[str, float] = {}
    other_total = 0.0

    for category, pct in percentages.items():
        if pct < threshold:
            other_total += pct
        else:
            collapsed[category] = pct

    if other_total > 0:
        collapsed["Other"] = other_total

    return collapsed


def collect_tables(
    folder: Path,
    rank: str,
    threshold: float,
    filter_rank: str | None = None,
    filter_name: str | None = None,
) -> pd.DataFrame:
    tsv_files = sorted(folder.glob("*.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No .tsv files found in {folder}")

    all_tables: Dict[str, Dict[str, float]] = {}

    for file_path in tsv_files:
        counts = read_and_aggregate_tsv(
            file_path=file_path,
            rank=rank,
            filter_rank=filter_rank,
            filter_name=filter_name,
        )
        total = sum(counts.values())

        if total == 0:
            percentages = {}
        else:
            percentages = {
                k: (v / total) * 100.0
                for k, v in counts.items()
            }

        percentages = collapse_small_categories(percentages, threshold)
        all_tables[file_path.stem] = percentages

    df = pd.DataFrame.from_dict(all_tables, orient="index").fillna(0.0)

    if df.empty:
        raise ValueError("No valid data found.")

    # remove samples that have no remaining signal after filtering
    df = df.loc[df.sum(axis=1) > 0]

    if df.empty:
        raise ValueError("Filtering removed all data; no samples remain to plot.")

    order = list(df.sum(axis=0).sort_values(ascending=False).index)

    if "Other" in order:
        order = [c for c in order if c != "Other"] + ["Other"]

    return df[order]


def build_color_map(columns):
    cmap = plt.get_cmap("tab20")
    color_map = {}

    non_other = [c for c in columns if c != "Other"]

    for i, col in enumerate(non_other):
        color_map[col] = cmap(i % cmap.N)

    if "Other" in columns:
        color_map["Other"] = (0.8, 0.8, 0.8, 1.0)

    return color_map


def plot_stacked_barh(
    df: pd.DataFrame,
    rank: str,
    output_prefix: Path,
    threshold: float,
    filter_rank: str | None = None,
    filter_name: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(14, max(4, 0.6 * len(df))))

    color_map = build_color_map(df.columns)
    left = pd.Series(0.0, index=df.index)

    for category in df.columns:
        values = df[category]
        ax.barh(
            df.index,
            values,
            left=left,
            label=category,
            color=color_map[category],
            edgecolor="white",
            linewidth=0.5,
        )
        left += values

    ax.set_xlabel("Relative abundance (%)")
    ax.set_xlim(0, 100)

    title = f"{RANK_LABELS[rank]} composition (threshold < {threshold:g}% → Other)"
    if filter_rank is not None and filter_name is not None:
        title += f" | filtered to {RANK_LABELS[filter_rank]} = {filter_name}"
    ax.set_title(title)

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()

    fig.savefig(output_prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("-r", "--rank", required=True, choices=VALID_RANKS)
    parser.add_argument("-o", "--output-prefix", required=True, type=Path)
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.0,
        help='Merge categories < threshold (%) into "Other"',
    )
    parser.add_argument(
        "--filter-rank",
        choices=VALID_RANKS,
        default=None,
        help="Optional taxonomy rank used to filter rows. Must be used together with --filter-name.",
    )
    parser.add_argument(
        "--filter-name",
        default=None,
        help="Optional taxonomy name used to filter rows. Must be used together with --filter-rank.",
    )

    args = parser.parse_args()

    if not (0 <= args.threshold <= 100):
        raise ValueError("Threshold must be between 0 and 100.")

    if (args.filter_rank is None) != (args.filter_name is None):
        raise ValueError(
            "--filter-rank and --filter-name must be provided together."
        )

    df = collect_tables(
        folder=args.folder,
        rank=args.rank,
        threshold=args.threshold,
        filter_rank=args.filter_rank,
        filter_name=args.filter_name,
    )

    plot_stacked_barh(
        df=df,
        rank=args.rank,
        output_prefix=args.output_prefix,
        threshold=args.threshold,
        filter_rank=args.filter_rank,
        filter_name=args.filter_name,
    )


if __name__ == "__main__":
    main()
