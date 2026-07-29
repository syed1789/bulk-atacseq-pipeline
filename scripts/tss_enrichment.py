#!/usr/bin/env python3
"""Snakemake script: aggregate TSS enrichment profile from Tn5-shifted
insertion sites -- the standard ENCODE-style ATAC-seq QC metric. Uses the
full yeast annotation (7,126 transcripts), unlike the sparse 300-gene test
annotation in the bulk-chipseq-pipeline repo, so this one is a real,
informative enrichment test rather than a chance-level artifact."""
import bisect
import subprocess
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bam = snakemake.input.bam
gtf = snakemake.input.gtf
sample = snakemake.wildcards.sample

WINDOW = 1000
BIN = 20

# TSS per transcript: start coord for +, end coord for - (GTF is 1-based).
tss_by_chrom = defaultdict(list)
with open(gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[2] != "transcript":
            continue
        chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
        tss = start if strand == "+" else end
        tss_by_chrom[chrom].append((tss, strand))

# Tn5 insertion sites = 5' end of each read in the shifted BAM.
insertions_by_chrom = defaultdict(list)
out = subprocess.run(
    ["samtools", "view", "-F", "0x104", bam],
    capture_output=True, text=True, check=True,
).stdout
for line in out.splitlines():
    f = line.split("\t")
    chrom = f[2]
    flag = int(f[1])
    pos = int(f[3])
    is_reverse = bool(flag & 0x10)
    cigar = f[5]
    if is_reverse:
        # 5' end of a reverse-strand read is at its rightmost mapped base.
        length = sum(int(n) for n, op in __import__("re").findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")
        insertions_by_chrom[chrom].append(pos + length - 1)
    else:
        insertions_by_chrom[chrom].append(pos)

for c in insertions_by_chrom:
    insertions_by_chrom[c].sort()

n_bins = (2 * WINDOW) // BIN
profile = [0] * n_bins

for chrom, tss_list in tss_by_chrom.items():
    positions = insertions_by_chrom.get(chrom, [])
    if not positions:
        continue
    for tss, strand in tss_list:
        lo, hi = tss - WINDOW, tss + WINDOW
        i = bisect.bisect_left(positions, lo)
        j = bisect.bisect_right(positions, hi)
        for pos in positions[i:j]:
            rel = pos - tss if strand == "+" else tss - pos
            b = (rel + WINDOW) // BIN
            if 0 <= b < n_bins:
                profile[b] += 1

flank_n = 100 // BIN  # outer 100bp on each side of the 2kb window
flank_bins = profile[:flank_n] + profile[-flank_n:]
flank_mean = sum(flank_bins) / len(flank_bins) if flank_bins else 0
center_bins = profile[n_bins // 2 - 5 : n_bins // 2 + 5]
center_mean = sum(center_bins) / len(center_bins) if center_bins else 0
enrichment = center_mean / flank_mean if flank_mean > 0 else float("nan")

with open(snakemake.output.stats, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"n_transcripts\t{len(sum(tss_by_chrom.values(), []))}\n")
    f.write(f"tss_enrichment_score\t{enrichment:.2f}\n")

print(f"{sample}: TSS enrichment score = {enrichment:.2f} "
      f"(center signal / flank signal, ENCODE-style; >5-7 is typically considered good)")

x = [(-WINDOW + i * BIN) for i in range(n_bins)]
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(x, profile, color="#c0453a")
ax.axvline(0, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("Distance from TSS (bp)")
ax.set_ylabel("Aggregate Tn5 insertions")
ax.set_title(f"{sample} — TSS enrichment (score={enrichment:.2f})")
fig.tight_layout()
fig.savefig(snakemake.output.plot, dpi=150)
