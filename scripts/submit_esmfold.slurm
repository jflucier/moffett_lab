#!/bin/bash
#SBATCH --job-name=esm_bins
#SBATCH --account=def-moffettp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=/home/jflucier/links/scratch/20260825_folds_arabidopsis/bin_testset/logs_bin_%A_%a.out
# Set a uniform resource ceiling that safely covers the largest 2501-5000aa bins
#SBATCH --gpus=h100:1
#SBATCH --mem=62G


########
# To submit use the following call:
# sbatch submit_esmfold.slurm --array=0-xx --gpus=[H100-1g.10gb|H100-2g.20gb|H100-3g.40gb|h100:1] -i /pth/to/fasta_dir -o /path/to/result_dir -c chunk_size

# Define the absolute directory mapping based on your path
PYTHON_SCRIPT="/home/jflucier/links/scratch/programs/moffett_lab/python/esmfold_screening.py"
MY_BAIT="MTTSRFATFDIESETGLTPGAYPAPLPTLEQQLHDRNAILAAIPGLARTKLDAPTLKRAFANFLLTLGMVGTTSKGSYEELIIPPVKGMGSSTGFRARELVQIITSSPAPPGFDGNQTLRQFARPYAPQVQNMIAQGKFKTNLYDKYGKSVGAPPHVCIDFNDAMDLQMFHSTAEFESAHKVRELAIAEAAARDNAPRPAANPRAAKPVIGQTAPAFHSGDAAKQSGGQPVNIKPSLAQSAGFDSHRPPPETPPRASTPSSQKSGQSGQTIIQPPASHGILSSALGSHKSTPHASPQQTPKK"

CHUNK_ARGS=""
FASTA_DIR=""
RESULTS_DIR=""

TEMP=$(getopt -o i:o:c: --long in:,out:,chunk: -n 'submit_esmfold.slurm' -- "$@")
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -i|--in) FASTA_DIR="$2" ; shift 2 ;;
        -o|--out) RESULTS_DIR="$2" ; shift 2 ;;
        -c|--chunk) CHUNK_ARGS="$2" ; shift 2 ;;
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

# need to fetch fasta based on $SLURM_ARRAY_INDEX
#FULL_FASTA_PATH="${DATA_DIR}/${FASTA}"
FULL_FASTA_PATH=$(ls ${FASTA_DIR}/*.fasta | awk ...)
RESULTS_DIR="${DATA_DIR}/esmfold_results"
mkdir -p

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


