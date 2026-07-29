#!/usr/bin/env python3
"""Snakemake script: fragment-size distribution from a paired-end BAM --
the classic ATAC-seq QC signature (nucleosome-free fragments <100bp,
periodic mono-/di-nucleosome peaks) computed from properly-paired,
first-in-pair TLEN values."""
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bam = snakemake.input[0]
sample = snakemake.wildcards.sample

# -f 0x42: paired + proper pair + first in pair (each fragment counted once)
out = subprocess.run(
    ["samtools", "view", "-f", "0x42", bam],
    capture_output=True, text=True, check=True,
).stdout

sizes = []
for line in out.splitlines():
    fields = line.split("\t")
    tlen = abs(int(fields[8]))
    if 0 < tlen < 1000:
        sizes.append(tlen)

n = len(sizes)
sizes_sorted = sorted(sizes)
median = sizes_sorted[n // 2]
pct_nfr = 100 * sum(1 for v in sizes if v < 100) / n
pct_mono = 100 * sum(1 for v in sizes if 147 <= v <= 294) / n

with open(snakemake.output.stats, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"n_fragments\t{n}\n")
    f.write(f"median_fragment_size\t{median}\n")
    f.write(f"pct_nucleosome_free_lt100bp\t{pct_nfr:.1f}\n")
    f.write(f"pct_mononucleosome_147_294bp\t{pct_mono:.1f}\n")

print(f"{sample}: n={n}, median={median}bp, "
      f"NFR(<100bp)={pct_nfr:.1f}%, mono-nuc(147-294bp)={pct_mono:.1f}%")

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(sizes, bins=range(0, 800, 5), color="#2a6f77")
ax.axvspan(0, 100, alpha=0.15, color="green", label="nucleosome-free (<100bp)")
ax.axvspan(147, 294, alpha=0.15, color="orange", label="mono-nucleosome (147-294bp)")
ax.set_xlabel("Fragment size (bp)")
ax.set_ylabel("Count")
ax.set_title(f"{sample} — ATAC-seq fragment size distribution")
ax.legend()
fig.tight_layout()
fig.savefig(snakemake.output.plot, dpi=150)
