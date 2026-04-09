#!/bin/bash
# filepath: /vast/scratch/users/putri.g/submit_apptainer_jobs.sh

#SBATCH --job-name=apptainer-pull_array
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=regular
#SBATCH --time=01:00:00
#SBATCH --array=1-31
#SBATCH --output=/vast/scratch/users/putri.g/slurm_log/apptainer/pull_%a.out
#SBATCH --error=/vast/scratch/users/putri.g/slurm_log/apptainer/pull_%a.err

echo "Running Array Job ID="$SLURM_ARRAY_JOB_ID", Array Index="$SLURM_ARRAY_TASK_ID

module load apptainer

IMAGES_FILE="/vast/scratch/users/putri.g/cytobenchmark/apptainer_images/images.txt"
WORK_DIR="/vast/scratch/users/putri.g/nextflow/apptainer_cache"

# Each job gets its own isolated cache
export APPTAINER_CACHEDIR="$WORK_DIR/cache_$SLURM_ARRAY_JOB_ID_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$APPTAINER_CACHEDIR"

# Read image for this task
image_url=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IMAGES_FILE")
if [ -z "$image_url" ]; then
    echo "Error: No image at line $SLURM_ARRAY_TASK_ID in $IMAGES_FILE" >&2
    exit 1
fi

echo "Job $SLURM_ARRAY_TASK_ID: Building $image_url"
echo "Cache directory: $APPTAINER_CACHEDIR"

final_dir="$WORK_DIR/images"
mkdir -p "$final_dir"

# Extract full name from URL: docker://ghcr.io/path/to/image:tag -> ghcr.io-path-to-image-tag
image_name=$(echo "$image_url" | sed 's|docker://||' | sed 's|/|-|g' | sed 's|:|-|g')
output_file="$final_dir/${image_name}.img"

echo "Output file will be: $output_file"

# Build directly to .img format
if apptainer build "$output_file" "$image_url"; then
    echo "Successfully built: $output_file"
else
    echo "Error: Failed to build $image_url" >&2
    exit 1
fi
