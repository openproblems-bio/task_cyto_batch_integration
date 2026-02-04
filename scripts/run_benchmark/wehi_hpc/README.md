# Running the benchmark on WEHI hpc

This directory contains all the scripts needed to submit the benchmark job on WEHI hpc.

The pipeline will use the config file in `scripts/labels_tw_wehi.config`.

Few things to note if you want to replicate this on different slurm system:

1. Before submitting the job, make sure to run `scripts/run_benchmark/wehi_hpc/build_apptainer_images.sh` script first. See sections below for more info.
2. Make sure you have setup seqera tower's compute environment to point to the HPC cluster. 
   You must select Slurm as the as environment. Have a look at the existing one and clone the relevant settings.
3. Change the workspace ID to point to the new environment created in step 2. 
4. Create all the necessary temp folders. See the config file.
5. Run tower agent in the background, using tmux or something.


## Apptainr image build

The `scripts/run_benchmark/wehi_hpc/build_apptainer_images.sh` will pull the docker images
from url in `scripts/run_benchmark/wehi_hpc/images.txt` and build them in parallel (one job = 1 image).
It will automatically rename the images to something that the pipeline can understand.
So something like: `ghcr.io-openproblems-bio-openproblems-utils-extract_uns_metadata-build_main.img`.

This script must be run before you run the full benchmark as otherwise the head job will be overwhelmed with
the number of images it has to build and just die.

Why not get each worker to build its own image? It is because the apptainer cache directory
is shared between different jobs. 
So there will be locking contention. 
This may be fixable by setting each job to use its own cache, but i don't have time to test it.

If there is a new method or metric, make sure that is added to the `scripts/run_benchmark/wehi_hpc/images.txt` file.
If we have remove a method or metric, make sure you remove it from the file as well.

If you only update one or two methods or metrics, there is no needed to pull and rebuild everything.
Just build those that were updated.
Create a new text file with the url of the images you need to build and run the script.
Make sure you change the filename of the `scripts/run_benchmark/wehi_hpc/images.txt`.

The apptainer images must be stored directly in the `APPTAINER_CACHE` directory.
If you want to test whether the job uses the images properly, run the `scripts/run_benchmark/wehi_hpc/run_subset_hpc.sh`
script which by default will only run one metric and one method.
If you see `cache` folder or `*.lock` file in the `APPTAINER_CACHE` directory, then you know
the head job is repulling and rebuilding the images.

