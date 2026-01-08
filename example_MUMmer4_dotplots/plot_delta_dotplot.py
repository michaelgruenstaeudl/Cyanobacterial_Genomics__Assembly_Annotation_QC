#!/usr/bin/env python3
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt

COLOR_MAP = {
    "+": "blue",
    "-": "red",
}

LINE_WIDTH = 2.4
POINT_SIZE = 6.0

def usage():
    sys.stderr.write(
        "Usage: plot_delta_dotplot.py <filtered.delta> <output_basename>\n"
    )
    sys.exit(1)

def read_delta(delta_path: str) -> pd.DataFrame:
    rows = []
    ref_id = None
    qry_id = None

    with open(delta_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    if not lines:
        raise ValueError("Delta file is empty")

    lines = lines[1:]

    header_re = re.compile(r"^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)")

    for line in lines:
        line = line.strip()
        if not line or line == "0":
            continue

        m = header_re.match(line)
        if m:
            ref_id = m.group(1)
            qry_id = m.group(2)
            continue

        parts = line.split()
        if len(parts) == 7 and all(p.lstrip("-").isdigit() for p in parts):
            rs, re_, qs, qe, err = map(int, parts[:5])
            rows.append(
                {
                    "rs": rs,
                    "re": re_,
                    "qs": qs,
                    "qe": qe,
                    "error": err,
                    "rid": ref_id,
                    "qid": qry_id,
                    "strand": "+" if (qe - qs) > 0 else "-",
                }
            )

    if not rows:
        raise ValueError("No alignment records parsed from delta file")

    return pd.DataFrame(rows)

def dotplot(df: pd.DataFrame, out_base: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))

    for strand, g in df.groupby("strand", sort=True):
        color = COLOR_MAP.get(strand, "black")
        for _, r in g.iterrows():
            ax.plot(
                [r["rs"], r["re"]],
                [r["qs"], r["qe"]],
                linewidth=LINE_WIDTH,
                color=color,
            )
        ax.scatter(g["rs"], g["qs"], s=POINT_SIZE, alpha=0.5, color=color)
        ax.scatter(g["re"], g["qe"], s=POINT_SIZE, alpha=0.5, color=color)

    ax.set_xlabel("reference position (bp)")
    ax.set_ylabel("query position (bp)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(f"{out_base}.png", dpi=300)
    fig.savefig(f"{out_base}.svg")
    plt.close(fig)

def main():
    if len(sys.argv) != 3:
        usage()

    delta_file = sys.argv[1]
    out_base = sys.argv[2]

    df = read_delta(delta_file)
    dotplot(df, out_base)

if __name__ == "__main__":
    main()
