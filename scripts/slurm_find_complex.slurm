#!/bin/bash
#SBATCH --job-name=esmfold_screen
#SBATCH --account=def-moffettp
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/screen_%j.out

# --- Default Fallback Parameters ---
INPUT_DIR=""
OUTPUT_FILE=""

# --- Parse Arguments Passed via sbatch ---
# Usage: sbatch submit_screen.sh -i <input_dir> -o <output_file>
while getopts "i:o:" opt; do
  case ${opt} in
    i ) INPUT_DIR=$OPTARG ;;
    o ) OUTPUT_FILE=$OPTARG ;;
    \? ) echo "Usage: sbatch submit_screen.sh [-i input_dir] [-o output_file]" && exit 1 ;;
  esac
done

if [ -z "$INPUT_DIR" ]; then
    echo "ERROR: Input directory (-i) is required but was not provided." >&2
    echo "Usage: sbatch submit_screen.sh -i <input_dir> -o <output_file>" >&2
    exit 1
fi

if [ -z "$OUTPUT_FILE" ]; then
    echo "ERROR: Output file path (-o) is required but was not provided." >&2
    echo "Usage: sbatch submit_screen.sh -i <input_dir> -o <output_file>" >&2
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Provided input directory does not exist: $INPUT_DIR" >&2
    exit 1
fi

# 1. Create a logs directory if it doesn't exist
mkdir -p logs

# 2. Load required cluster modules
module load scipy-stack                 # Pre-loaded scientific stack (highly optimized numpy/scipy)

# 3. Activate your existing virtual environment
ESM_PATH=/home/jflucier/links/scratch/programs/esm
source ${ESM_PATH}/venv/bin/activate

# 5. Execute the parallelized screening script
# We pass $SLURM_CPUS_PER_TASK directly to tell the script how many parallel workers to open
python /home/jflucier/links/scratch/programs/moffett_lab/python/find_complexes.py \
    --input_dir $INPUT_DIR \
    --output $OUTPUT_FILE \
    --workers $SLURM_CPUS_PER_TASK

echo "Slurm Job $SLURM_JOB_ID complete."
