#!/bin/bash
#SBATCH --job-name=esm_bins
#SBATCH --account=def-moffettp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=/home/jflucier/links/scratch/20260825_folds_arabidopsis/bin_testset/logs_bin_%A_%a.out


# Set a uniform resource ceiling that safely covers the largest 2501-5000aa bins
#SBATCH --gpus=h100:1
#SBATCH --mem=62G



# Define the absolute directory mapping based on your path
DATA_DIR="/home/jflucier/links/scratch/20260825_folds_arabidopsis/bin_testset"
PYTHON_SCRIPT="/home/jflucier/links/scratch/programs/moffett_lab/python/esmfold_screening.py"
MY_BAIT="MTTSRFATFDIESETGLTPGAYPAPLPTLEQQLHDRNAILAAIPGLARTKLDAPTLKRAFANFLLTLGMVGTTSKGSYEELIIPPVKGMGSSTGFRARELVQIITSSPAPPGFDGNQTLRQFARPYAPQVQNMIAQGKFKTNLYDKYGKSVGAPPHVCIDFNDAMDLQMFHSTAEFESAHKVRELAIAEAAARDNAPRPAANPRAAKPVIGQTAPAFHSGDAAKQSGGQPVNIKPSLAQSAGFDSHRPPPETPPRASTPSSQKSGQSGQTIIQPPASHGILSSALGSHKSTPHASPQQTPKK"

CHUNK_ARGS=""

case $SLURM_ARRAY_TASK_ID in
    1) FASTA="bin_1-500aa.testset.fasta" ;;
    2) FASTA="bin_501-1000aa.testset.fasta";   CHUNK_ARGS="--chunk-size 128" ;;
    3) FASTA="bin_1001-1500aa.testset.fasta";  CHUNK_ARGS="--chunk-size 64"  ;; # Smaller chunk for larger sequences
    4) FASTA="bin_1501-2000aa.testset.fasta";  CHUNK_ARGS="--chunk-size 64"  ;;
    5) FASTA="bin_2001-2500aa.testset.fasta";  CHUNK_ARGS="--chunk-size 64"  ;;
    6) FASTA="bin_2501-5000aa.testset.fasta";  CHUNK_ARGS="--chunk-size 32"  ;; # Aggressive chunking for extreme targets
esac

FULL_FASTA_PATH="${DATA_DIR}/${FASTA}"
RESULTS_DIR="${DATA_DIR}/screening_results_${FASTA%.testset.fasta}"

# Load your cluster runtime environment
module load StdEnv/2023 python/3.10 cuda/12.2 gcc/12.3
source /home/jflucier/links/scratch/programs/esm/venv/bin/activate

export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"


# Define parameter variables
echo "Processing Task ID: $SLURM_ARRAY_TASK_ID"
echo "Target File: $FULL_FASTA_PATH"
echo "Allocated Resource Profile: $MIG_PROFILE"

# Execute your python screening script
python "$PYTHON_SCRIPT" \
    --bait "$MY_BAIT" \
    --fasta "$FULL_FASTA_PATH" \
    --outdir "$RESULTS_DIR" \
    $CHUNK_ARGS
