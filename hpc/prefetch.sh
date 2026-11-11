#!/usr/bin/env bash
# Run on a login node (internet access) so compute nodes can use the cache.
set -euo pipefail
cd "$(dirname "$0")/.."
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRATCH:-$HOME}/.apptainer}"
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
for img in \
  "quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0" \
  "quay.io/biocontainers/macs3:3.0.4--py313hda738de_0" \
  "quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0"
do
  echo "Caching $img"
  apptainer exec "docker://$img" true
done
bash scripts/download_data.sh
