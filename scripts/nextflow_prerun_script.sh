# paste me as pre-run script in Tower if setting up workflow run in WEHI HPC.
# load module first so the variables don't get overwritten
module load nextflow/25.04.2

# Tower pre-run script
export SHARED_SCRATCH="/vast/scratch/users/putri.g/nextflow"

export NXF_APPTAINER_CACHEDIR="$SHARED_SCRATCH/apptainer_cache"
export APPTAINER_CACHEDIR="$SHARED_SCRATCH/apptainer_cache"
export APPTAINER_TMPDIR="$SHARED_SCRATCH/apptainer_tmp"
export APPTAINER_LIBRARYDIR="$SHARED_SCRATCH/apptainer_library"
export SINGULARITY_CACHEDIR="$SHARED_SCRATCH/apptainer_cache"
export SINGULARITY_TMPDIR="$SHARED_SCRATCH/apptainer_tmp"
export TMPDIR="$SHARED_SCRATCH/apptainer_tmp"
export NXF_HOME="$SHARED_SCRATCH/nxf_home"
export NXF_TEMP="$SHARED_SCRATCH/nxf_tmp"
export HOME="$SHARED_SCRATCH/home"

mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$APPTAINER_LIBRARYDIR" "$NXF_HOME" "$NXF_TEMP" "$HOME"