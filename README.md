# Bulk ATAC-seq: Alignment → Tn5 Shift → QC → Peak Calling

[![CI](https://github.com/syed1789/bulk-atacseq-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/syed1789/bulk-atacseq-pipeline/actions/workflows/ci.yml)

**Author:** Syed Nurul Hasan, Ph.D.
**Contact:** nayanchoton1789@gmail.com | [ORCID: 0000-0002-4564-872X](https://orcid.org/0000-0002-4564-872X)

---

## Why this project

The last public gap from my Ph.D.'s genome-wide profiling work (ChIP-seq +
RNA-seq + ATAC-seq in Tregs, POSTECH): bulk ATAC-seq specifically, distinct
from the single-cell ATAC-seq already in `bioinformatics-portfolio` and
`scATAC-WNN-multimodal`, and from `cutandrun-pipeline`'s CUT&Tag/CUT&RUN.
Bulk ATAC-seq has its own QC conventions — mitochondrial-read filtering,
Tn5 insertion-site shift, fragment-size periodicity, TSS enrichment — none
of which the other repos exercise, so this one builds the real thing rather
than relabeling an adjacent pipeline.

## Pipeline

```
Paired-end FASTQ (osmotic-stress time course, 100k read pairs x4 samples)
        │
        ▼
1. bowtie2_index / bowtie2_align  (Docker: biocontainers/bowtie2)
   --very-sensitive -X 2000 (wide fragment range for nucleosome-spanning fragments)
        │
        ▼
2. filter_mito  (samtools) -- drop reads mapping to the mitochondrial genome
        │
        ▼
3. dedup  (samtools fixmate -> markdup -r) -- PCR duplicate removal
        │
        ▼
4. atac_shift  (Docker: biocontainers/deeptools, alignmentSieve --ATACshift)
   +4/-5bp Tn5 insertion-site correction
        │
        ├──▶ fragment_size_qc  -- nucleosome-free / mono-nucleosome periodicity
        ├──▶ tss_enrichment    -- ENCODE-style aggregate TSS enrichment score
        └──▶ macs3_callpeak    -- no Input control (ATAC compares against its
                                   own background, unlike ChIP-seq)
```

## Dataset

Paired-end FASTQ (downsampled to 100k read pairs/sample, osmotic-stress
time course: T0 vs. T15 minutes, 2 replicates each) + the full *S.
cerevisiae* reference genome and annotation (7,126 transcripts), from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets)
(`atacseq` branch), derived from **GSE66386** — Schep AN, Buenrostro JD,
Denny SK, et al. "Structured nucleosome fingerprints enable high-resolution
mapping of chromatin architecture within regulatory regions." *Genome
Research*. 2015 — one of the foundational ATAC-seq/NucleoATAC papers.

## Results (this run)

- Alignment: 73.6–84.2% overall rate (whole-genome reference, real numbers).
- **Mitochondrial contamination: 0.2–0.9%** across all 4 samples — low, as
  expected for yeast (whose mito genome is far smaller relative to the
  nuclear genome than in mammalian samples, where 10–50%+ is common).
  [`results/qc/mito_fractions.tsv`](results/qc/mito_fractions.tsv).
- **Fragment-size periodicity**, the classic ATAC-seq signature: 14–18% of
  fragments in the nucleosome-free range (<100bp), 37–42% in the
  mono-nucleosome range (147–294bp), median fragment size 176–189bp across
  samples. [`results/qc/*_fragment_sizes.png`](results/qc/).
- **TSS enrichment (ENCODE-style score): 1.16–1.27** across samples —
  computed exactly as ENCODE defines it, and reported honestly even though
  it's well below the 5–15+ scores typical in mammalian ATAC-seq. This
  isn't a failed experiment: yeast's genome is ~70% coding and gene-dense,
  so the "flanking" regions 500–1000bp from a TSS that this metric uses as
  background are themselves often inside or near another gene's regulatory
  region — the sparse-background assumption the ENCODE metric was
  calibrated for doesn't hold in a compact genome. The profile itself still
  shows real, above-flank signal centered on the TSS; the *ratio* just
  reads lower here for a structural reason, not a QC failure.
  [`results/qc/*_tss_enrichment.png`](results/qc/).
- MACS3 peak calling (`-q 0.05`, BAMPE, no control): **290–512 peaks** per
  sample — same order of magnitude as nf-core/atacseq's own documented
  expected output for this exact dataset (627–1135 peaks with their more
  elaborate merged-replicate broadPeak pipeline; my narrower single-sample
  narrowPeak approach naturally calls fewer).

## Running it

```bash
bash scripts/download_data.sh

python3 -m venv .venv && source .venv/bin/activate
pip install snakemake matplotlib

snakemake -s workflow/Snakefile --cores 4 -p

# SLURM cluster execution
snakemake -s workflow/Snakefile --workflow-profile workflow/profiles/slurm
```

## Key Features

- Real bowtie2 alignment + MACS3 peak calling, containerized
- Full ATAC-seq-specific QC: mito filtering, duplicate removal, Tn5 shift
  (deeptools), fragment-size periodicity, TSS enrichment — not a relabeled
  ChIP-seq pipeline
- A QC metric (TSS enrichment) reported honestly with the structural reason
  it reads lower than mammalian benchmarks, rather than omitted
- An actual Snakemake DAG (35 rule invocations) with a SLURM executor profile

## References

1. Schep AN, Buenrostro JD, Denny SK, Schwartz K, Sherlock G, Greenleaf WJ.
   Structured nucleosome fingerprints enable high-resolution mapping of
   chromatin architecture within regulatory regions. *Genome Research*.
   2015;25(11):1757-1770.
2. Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ.
   Transposition of native chromatin for fast and sensitive epigenomic
   profiling of open chromatin, DNA-binding proteins and nucleosome
   position. *Nature Methods*. 2013;10(12):1213-1218.
3. Ramírez F, Ryan DP, Grüning B, et al. deepTools2: a next generation web
   server for deep-sequencing data analysis. *Nucleic Acids Research*.
   2016;44(W1):W160-W165.

## License

MIT — see [LICENSE](LICENSE).

## Running on HPC (SLURM) and AWS

Includes `CONTAINER_RUNTIME=apptainer` switch for the Snakefile, a SLURM head-job script, and a one-shot EC2 launcher. See [HPC_AWS.md](HPC_AWS.md) for steps and what was verified.
