#!/bin/bash
# Fetches paired-end ATAC-seq FASTQ (osmotic-stress time course, 100k read
# pairs/sample) and the S. cerevisiae reference genome + full annotation
# from nf-core/test-datasets (atacseq branch), derived from Schep et al.,
# Genome Research 2015 (GSE66386).
set -euo pipefail
cd "$(dirname "$0")/.."

BASE="https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/testdata"
REFBASE="https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/reference"

mkdir -p data/fastq data/reference

for f in SRR1822153_1 SRR1822153_2 SRR1822154_1 SRR1822154_2 \
         SRR1822157_1 SRR1822157_2 SRR1822158_1 SRR1822158_2; do
  out="data/fastq/${f}.fastq.gz"
  [ -f "$out" ] && { echo "Already downloaded: $out"; continue; }
  echo "Downloading $out"
  curl -sL --fail -o "$out" "${BASE}/${f}.fastq.gz"
done

if [ ! -f data/reference/genome.fa ]; then
  echo "Downloading S. cerevisiae reference genome"
  curl -sL --fail -o data/reference/genome.fa "${REFBASE}/genome.fa"
fi
if [ ! -f data/reference/genes.gtf ]; then
  echo "Downloading full annotation (7,126 transcripts)"
  curl -sL --fail -o data/reference/genes.gtf "${REFBASE}/genes.gtf"
fi

echo "Done."
