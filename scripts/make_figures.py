#!/usr/bin/env python3
"""Peak-level figures from MACS narrowPeak files -> results/figures/.

Usage: python scripts/make_figures.py   (run from the repo root after the pipeline)
"""
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PEAKDIR = Path("results/macs3"); F = Path("results/figures"); F.mkdir(exist_ok=True)
COLS = ["chrom", "start", "end", "name", "score", "strand", "fc", "neglog10p", "neglog10q", "summit"]
peaks = {p.name.replace("_peaks.narrowPeak", ""): pd.read_csv(p, sep="\t", header=None, names=COLS)
         for p in sorted(PEAKDIR.glob("*_peaks.narrowPeak"))}

def overlap(a, b):
    n = 0
    for c, g in a.groupby("chrom"):
        o = b[b.chrom == c]
        s, e = o.start.values, o.end.values
        n += sum(((s < r.end) & (e > r.start)).any() for r in g.itertuples())
    return n

names = list(peaks)
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
axs[0].bar(names, [len(p) for p in peaks.values()], color="#4c72b0")
axs[0].set_ylabel("peaks"); axs[0].set_title("Peaks called")
for n, p in peaks.items():
    axs[1].hist(p.end - p.start, bins=40, histtype="step", label=n, lw=1.3)
    axs[2].hist(p.fc, bins=40, histtype="step", label=n, lw=1.3)
axs[1].set_xlabel("peak width (bp)"); axs[1].set_title("Peak width"); axs[1].legend(frameon=False, fontsize=7)
axs[2].set_xlabel("fold enrichment at summit"); axs[2].set_title("Peak strength")
for a in axs: a.tick_params(axis="x", labelrotation=30)
fig.tight_layout(); fig.savefig(F / "peak_summary.png", dpi=130); plt.close(fig)

# pairwise overlap (fraction of row's peaks overlapping column's peaks)
m = np.array([[overlap(peaks[a], peaks[b]) / len(peaks[a]) * 100 for b in names] for a in names])
fig, ax = plt.subplots(figsize=(1.2 * len(names) + 2.5, 1.0 * len(names) + 2))
im = ax.imshow(m, cmap="viridis", vmin=0, vmax=100)
ax.set_xticks(range(len(names)), names, rotation=30, ha="right"); ax.set_yticks(range(len(names)), names)
for i in range(len(names)):
    for j in range(len(names)):
        ax.text(j, i, f"{m[i,j]:.0f}", ha="center", va="center", color="w" if m[i,j] < 60 else "k")
ax.set_title("% of row's peaks overlapping column's peaks"); fig.colorbar(im, ax=ax, fraction=.046)
fig.tight_layout(); fig.savefig(F / "peak_overlap.png", dpi=130); plt.close(fig)

mt = pd.read_csv("results/qc/mito_fractions.tsv", sep="	")
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar(mt["sample"], mt.pct_mito, color="#dd8452"); ax.set_ylabel("% mitochondrial reads")
ax.set_title("Mitochondrial read fraction"); plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
fig.tight_layout(); fig.savefig(F / "mito_fraction.png", dpi=130); plt.close(fig)

print({n: len(p) for n, p in peaks.items()})
