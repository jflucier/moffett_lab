#!/bin/bash
#SBATCH --job-name=at_screen
#SBATCH --account=def-moffettp   # Replace with your allocation account
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=screen_%j_%a.out

# Usage (with job array, one task per shard):
#   sbatch --gres=gpu:<mig_type>:1 --array=0-N slurm_esmfold_screen.sh <BAIT_SEQ> <SHARD_DIR> <OUTDIR> [CHUNK_SIZE]
#
# Usage (single job, no sharding):
#   sbatch --gres=gpu:<mig_type>:1 slurm_esmfold_screen.sh <BAIT_SEQ> <SHARD_DIR> <OUTDIR> [CHUNK_SIZE]
#   (picks shard_0000.fasta when SLURM_ARRAY_TASK_ID is unset)
#
# Positional arguments:
#   $1  BAIT_SEQ    Amino acid sequence of the bait protein
#   $2  SHARD_DIR   Directory containing shard_0000.fasta, shard_0001.fasta, ... (see shard_fasta.py)
#   $3  OUTDIR      Directory to save generated PDB files (shared across all array tasks; filenames are per-target so no collisions)
#   $4  CHUNK_SIZE  (optional) axial attention chunk size, omit for no chunking

set -euo pipefail

BAIT_SEQ="$1"
SHARD_DIR="$2"
OUTDIR="$3"
CHUNK_SIZE="${4:-}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
SHARD_ID=$(printf '%04d' "$TASK_ID")
FASTA_PATH="${SHARD_DIR}/shard_${SHARD_ID}.fasta"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/../python/esmfold_screening.py"

if [[ ! -f "$FASTA_PATH" ]]; then
    echo "ERROR: shard file not found: $FASTA_PATH" >&2
    exit 1
fi

# Load environment
module load StdEnv/2023 python/3.10 cuda/12.2 gcc/12.3
source /home/jflucier/links/scratch/programs/esm/venv/bin/activate

CHUNK_ARGS=()
if [[ -n "$CHUNK_SIZE" ]]; then
    CHUNK_ARGS=(--chunk-size "$CHUNK_SIZE")
fi

echo "Array task ${TASK_ID}: processing ${FASTA_PATH}"

# Execute with arguments
python "$PYTHON_SCRIPT" \
    --bait "$BAIT_SEQ" \
    --fasta "$FASTA_PATH" \
    --outdir "$OUTDIR" \
    "${CHUNK_ARGS[@]}"
