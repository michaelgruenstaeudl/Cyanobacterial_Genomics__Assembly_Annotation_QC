#!/usr/bin/env python3
import sys
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
        "Usage: plot_coords_dotplot.py <coords.tsv> <output_basename>\n"
    )
    sys.exit(1)

def read_coords_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)

    if df.shape[1] < 4:
        raise ValueError("Input file does not appear to be a valid show-coords TSV")

    df = df.rename(columns={0: "rs", 1: "re", 2: "qs", 3: "qe"}).copy()
    for c in ["rs", "re", "qs", "qe"]:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(int)

    df["strand"] = df.apply(
        lambda r: "+" if (r["qe"] - r["qs"]) > 0 else "-", axis=1
    )
    return df

def dotplot_segments(df: pd.DataFrame, out_base: str) -> None:
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

    coords_tsv = sys.argv[1]
    out_base = sys.argv[2]

    df = read_coords_tsv(coords_tsv)
    dotplot_segments(df, out_base)

if __name__ == "__main__":
    main()
