#!/usr/bin/env python3
"""Portable entry point for the preserved historical score-analysis functions."""
import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd
from score_analysis_original import (
    analyze_score_distribution, identify_top_binders, save_top_peptides,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--top", type=int, default=100)
    args = parser.parse_args()
    df = pd.read_csv(args.input.resolve())
    required = {"Ligand_ID", "Ligand", "Root", "Part", "Score_KPNA4", "Score_MC1R"}
    if required - set(df.columns):
        parser.error(f"Missing columns: {sorted(required - set(df.columns))}")
    if df.empty or args.top < 1:
        parser.error("Input must contain records and --top must be positive")
    scores = df[["Score_KPNA4", "Score_MC1R"]].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(scores.to_numpy()).all():
        parser.error("Scores must be finite and non-missing")
    if df["Ligand_ID"].isna().any() or df["Ligand_ID"].duplicated().any():
        parser.error("Ligand_ID must be present and unique")
    df[scores.columns] = scores
    args.output.mkdir(parents=True, exist_ok=True)
    os.chdir(args.output.resolve())
    scores.describe().to_csv("score_summary.csv")
    df.groupby("Root", dropna=False).size().to_csv("source_counts.csv", header=["count"])
    analyze_score_distribution(df)
    save_top_peptides(*identify_top_binders(df, min(args.top, len(df))))
    print(f"Analyzed {len(df):,} rows; outputs: {Path.cwd()}")


if __name__ == "__main__":
    main()
