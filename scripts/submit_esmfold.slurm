#!/bin/bash
#SBATCH --job-name=esm_bins
#SBATCH --account=def-moffettp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=/home/jflucier/links/scratch/20260825_folds_arabidopsis/arabidopsis_openprot_bins_500aa_v2/logs/esm_bins_%A_%a.out
#SBATCH --mem=31G

########
# To submit use the following call:
# sbatch --array=0-xx --gpus=[H100-1g.10gb|H100-2g.20gb|H100-3g.40gb|h100:1] submit_esmfold.slurm -i /pth/to/fasta_dir -o /path/to/result_dir -c chunk_size

# Define the absolute directory mapping based on your path
PYTHON_SCRIPT="/home/jflucier/links/scratch/programs/moffett_lab/python/esmfold_screening.py"
MY_BAIT="MTTSRFATFDIESETGLTPGAYPAPLPTLEQQLHDRNAILAAIPGLARTKLDAPTLKRAFANFLLTLGMVGTTSKGSYEELIIPPVKGMGSSTGFRARELVQIITSSPAPPGFDGNQTLRQFARPYAPQVQNMIAQGKFKTNLYDKYGKSVGAPPHVCIDFNDAMDLQMFHSTAEFESAHKVRELAIAEAAARDNAPRPAANPRAAKPVIGQTAPAFHSGDAAKQSGGQPVNIKPSLAQSAGFDSHRPPPETPPRASTPSSQKSGQSGQTIIQPPASHGILSSALGSHKSTPHASPQQTPKK"

CHUNK=""
FASTA_DIR=""
RESULTS_DIR=""

TEMP=$(getopt -o i:o:c: --long in:,out:,chunk: -n 'submit_esmfold.slurm' -- "$@")
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -i|--in) FASTA_DIR="$2" ; shift 2 ;;
        -o|--out) RESULTS_DIR="$2" ; shift 2 ;;
        -c|--chunk) CHUNK="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) echo "Invalid Option" ; exit 1 ;;
    esac
done

if [ -z "$FASTA_DIR" ]; then
    echo "Error: fasta dir must be provided."
    exit 1
fi
if [ -z "$RESULTS_DIR" ]; then
    RESULTS_DIR="${FASTA_DIR}/esmfold_results"
fi

if [ -n "$CHUNK" ]; then
    if [ "$CHUNK" -eq 32 ] || [ "$CHUNK" -eq 64 ] || [ "$CHUNK" -eq 128 ]; then
        CHUNK_ARGS="--chunk-size $CHUNK"
    else
        echo "Error: Invalid chunk size '$CHUNK'. Allowed values are 32, 64, or 128."
        exit 1
    fi
else
    CHUNK_ARGS=""
fi

# need to fetch fasta based on $SLURM_ARRAY_INDEX
shopt -s nullglob
FASTA_FILES=($(ls -1 "$FASTA_DIR"/*.fasta | sort))

FULL_FASTA_PATH="${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}"

# Safety check if the target index is out of bounds for the number of files available
if [ -z "$FULL_FASTA_PATH" ]; then
    echo "Error: No FASTA file found for Slurm Task ID $SLURM_ARRAY_TASK_ID in $FASTA_DIR"
    exit 1
fi

mkdir -p ${RESULTS_DIR}

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


