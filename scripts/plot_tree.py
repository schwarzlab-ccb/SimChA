#!/usr/bin/env python3

import argparse
from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd


REQUIRED_COLUMNS = ["sample_id", "parent_id", "dist", "fitness"]


def get_samples_data(input_file):
    samples_df = pd.read_csv(input_file, sep="\t")
    missing = [col for col in REQUIRED_COLUMNS if col not in samples_df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s) in {input_file}: {', '.join(missing)}"
        )

    samples_df = samples_df.copy()
    samples_df["sample_id"] = samples_df["sample_id"].astype(str)
    samples_df["parent_id"] = samples_df["parent_id"].astype(str)
    samples_df["dist"] = pd.to_numeric(samples_df["dist"], errors="raise")
    samples_df["fitness"] = pd.to_numeric(samples_df["fitness"], errors="raise")
    return samples_df


def add_cumulative_distance(sample_df):
    if sample_df["sample_id"].duplicated().any():
        raise ValueError("sample_id must be unique to compute cumulative distances")

    indexed = sample_df.set_index("sample_id", drop=False)
    cache = {}
    visiting = set()

    def dist_to_root(node_id):
        if node_id in cache:
            return cache[node_id]
        if node_id in visiting:
            raise ValueError(f"Cycle detected in tree at node {node_id}")

        visiting.add(node_id)
        row = indexed.loc[node_id]
        edge_dist = float(row["dist"])
        parent_id = row["parent_id"]

        if parent_id == node_id or parent_id not in indexed.index:
            total = edge_dist
        else:
            total = dist_to_root(parent_id) + edge_dist

        cache[node_id] = total
        visiting.remove(node_id)
        return total

    out_df = sample_df.copy()
    out_df["dist_cumulative"] = out_df["sample_id"].map(dist_to_root)
    return out_df


def plot_tree(sample_df, title="", dpi=100):
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=dpi)

    sample_df = add_cumulative_distance(sample_df)

    indexed = sample_df.set_index("sample_id", drop=False)

    # Draw branch segments from parent to child where parent exists in the table.
    for _, row in sample_df.iterrows():
        parent_id = row["parent_id"]
        if parent_id == row["sample_id"] or parent_id not in indexed.index:
            continue

        parent = indexed.loc[parent_id]
        ax.plot(
            [parent["dist_cumulative"], row["dist_cumulative"]],
            [parent["fitness"], row["fitness"]],
            color="#444444",
            linewidth=1.1,
            alpha=0.8,
            zorder=1,
        )

    ax.scatter(
        sample_df["dist_cumulative"],
        sample_df["fitness"],
        s=24,
        c="#1f77b4",
        edgecolors="white",
        linewidths=0.5,
        alpha=0.95,
        zorder=2,
    )

    ax.set_xlabel("mutation count")
    ax.set_ylabel("fitness")
    ax.set_title(title if title else "Evolutionary tree (time vs. fitness)")
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    return fig, ax


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot tree from samples.tsv using dist (x) and fitness (y)"
    )
    parser.add_argument("-I", "--input", required=True, help="Path to samples.tsv")
    parser.add_argument("-O", "--output", default="", help="Output image path")
    parser.add_argument("--dpi", default=100, type=int, help="Output DPI")
    args = parser.parse_args()

    sample_df = get_samples_data(args.input)
    fig, ax = plot_tree(sample_df, title=f"Tree: {args.input}", dpi=args.dpi)

    out_path = (
        args.output
        if args.output
        else join(abspath(dirname(args.input)), "tree_dist_fitness.png")
    )
    print(f"Saving tree plot to {out_path}")
    fig.savefig(out_path, dpi=args.dpi)
