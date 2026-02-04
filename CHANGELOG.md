# task_cyto_batch_integration x.y.z

## BREAKING CHANGES

* Updated the file schemas (PR #94):
  - `src/api`: The common dataset is now split into two censored datasets: `censored_split1` and `censored_split2`. The method component should now be run twice, once for each split. The control methods return both splits at once. The metrics now receive both splits as input.
  - `src/control_methods`: Control methods now need to return both `integrated_split1` and `integrated_split2`.
  - `src/metrics`: Metrics now receive both `integrated_split1` and `integrated_split2` as input.
  - `src/workflows`: Update the workflows accordingly.

## NEW FUNCTIONALITY

* Initialised repository with API files, readme, and a few test resources.

* Added `methods/harmonypy` component (PR #3).

* Added `methods/limma_remove_batch_effect` component (PR #7).

* Added three negative control methods (PR #8):
  - `control_methods/shuffle_integration`
  - `control_methods/shuffle_integration_by_batch`
  - `control_methods/shuffle_integration_by_cell_type`

* Added `metrics/emd_per_samples` component (PR #9).

* Added `methods/combat` component (PR #25).

* Added `control_methods/no_integration ` (PR #26).

* Added `metrics/n_inconsistent_peaks` + some new helper functions (PR #31).

* Added `methods/cycombine_nocontrol` (PR #35).

* Added `metrics/average_batch_r2` + helper function (PR #36).

* Added `methods/gaussnorm` (PR #45).

* Added `methods/cytornom_controls` (PR #50).

* Added processing scripts for datasets (PR #55).

* Added `metrics/flowsom_mapping_similarity` (PR #59).

* Added `methods/mnn` (PR #75).

* Updated cyCombine (PR #78):
  * Added cyCombine with all controls (`methods/cycombine_all_controls`) 
  * Added cyCombine with one control - samples from only one condition (`methods/cycombine_one_control`).
  * Added parameters to tune cyCombine.

* Added `metrics/cms` (PR #79).

* Added `methods/batchadjust_one_control` and `methods/batchadjust_all_controls` (PR #82).

* Updated CytoNorm methods (PR #84):
  * Added CytoNorm with one control - samples from only one condition (`methods/cytonorm_one_control`)
  * Added CytoNorm with aggregate of samples as controls (`methods/cytonorm_no_controls`).
  * Added parameters to tune CytoNorm.

* Added CytoNorm correction to a goal batch (PR #92).
* Added cyCombine correction to a reference batch (PR #90).
* Added `metrics/bras` (PR #91).

* Added Seurat rPCA (PR #95).

* Added processing scripts for CLL dataset (PR #106).

* Added new metric `ratio_inconsistent_peaks` (PR #114).

* Added processing scripts for Lille dataset and remove ones for CLL dataset (PR #118).

* Added config and run scripts for running the benchmark on WEHI HPC (PR #119). 
  Refer to the pull request on Github to see what needs to be set up.
  TLDR; Setup the compute environment and the caching directory, then run the warmup job
  to pull and create the apptainer images (one job = one method + one metric).
  After that, run the main benchmarking job with the config file for the HPC system.

## MAJOR CHANGES

* Updated file schema (PR #18): 
  * Add is_control obs to indicate whether a cell should be used as control when correcting batch effect.
  * Removed donor_id obs from unintegrated censored.
  * Removed to_correct var from everything except common_dataset. 
  All datasets now will only contain markers that need to be corrected.

* Reupdated the file schema (PR #19):
  * Included changes in PR #21: data Processor component partitions cells between unintegrated(censored) 
  and validation.
  * Add back to_correct var to every file except integrated to reflect the real world 
  batch correction workflow better.
  * Reverted PR #18 to retain only the 1st two changes (add is_control and remove 
  donor_id from unintegrated_censored).

* Changed emd to emd_mean and emd_max (PR #27).

* Updated output anndata for methods to return all vars - corrected or not (PR #28).

* Added positive control (PR #30).

* Rewrote EMD metric so it no longer relies on implementation in cytonormpy (PR #48):
  * EMD is now calculated for all cell types as well.
  * Returned values are average and max across all cell types (exclude the one calculated agnostic of cell types),
  average and max across all donors (mean and max of values computed agnostic of cell types).

* Bump Viash version to 0.9.4 (PR #61, PR #62).

* Added EMD vertical global metric and split perfect integration into horizontal and vertical 
  for computing horizontal and vertical metrics (PR #63).

* Fix problems identified during a full run (PR #99).

* Update CytoVI (PR #114).

* Update CytoVI to normalise using minmax scaler fitted on batch 1 post correction (PR #119).

* Update batchadjust, cytonorm to use HPC temp dir if the environment variable is set or else
  default to what is set by viash. This is to prevent collision in temp files when the jobs are running (PR #119).  

* Update ratio inconsistent peaks to handle edge cases where methods return only zero values
  for a marker/cell type/donor combination, causing sd to be zero and division by zero (PR #119).

## MINOR CHANGES

* Enabled unit tests (PR #2).

* Added integrated test resource (PR #5).

* Updated file description in yaml file (PR #15).

* Updated scripts to enable running benchmark on seqera (PR #29).

* Updated project description (PR #42).

* Implemented function for writing FCS files from an anndata object (PR #49).

* Changed cytonorm and cycombine clustering to use lineage markers (PR #54).

* The metric `flowsom_mapping_similarity` now works at the cluster level (PR #68).

* Updated vertical EMD (PR #70):
  * Metric is computed in group (condition) specific manner. 
  * Split the metric into global and per cell type.
  * Refactored horizontal EMD.

* Schema defined in `src\api` has been modified to include dataset specific parameters (PR #71).

* Disabled parallelization of `metrics/cms` + cms distributions are now kept in the output AnnData object (PR #83).

* Add arguments for including/excluding methods and metrics in the benchmarking workflow (PR #100).

* Removed EMD max from calculation (PR #113).

* Tune the resource requirement for each method (PR #119).
  * Low time, mem, cpu for control methods.
  * Mid time, mem, cpu for most methods, except below.
  * High (or very high) time, mem, cpu for computationally ones like rPCA.


## BUG FIXES

* Change n_inconsistent_peaks output to float and add R2 to main.nf (PR #40).

* Changed naming 'gaussNorm' to 'gaussnorm' (PR #47).

* Fixed cyCombine so it now batch correct unnormalised data rather than normalised data (PR #58). 

* Fixed FlowSOM mapping similarity metric (PR #64).

* Fixed get_obs_var_for_integrated function in helper.R giving out error in mac (PR #65).

* Remove unlabelled cells when computing n_inconsistent_peaks metric (PR #69).

* Fix vertical EMD (PR #76):
  * Return NA if there are less than 2 samples per group in the data.
  * Refactoring and introduced "global" variable for output.

* Add MNN to run script (PR #77).

* Added `argument_groups` field in `methods/mnn` config file (PR #81).

* Fix unparseable characters in EMD metrics (PR #87).

* Fix missing anndata in yaml file and set the base_r docker image version to 1 instead of 1.0.0 (PR #89).

* Fix bug in control methods (PR #107, #108, #109).
  * All control methods are updated to cater the new schema.
  * All control methods are re-enabled. Selectively disable them when running the pipeline using method exclude.

* Fix bug in EMD where nan cannot be written out and added sklearn dependency for cytovi (PR #110). 

* Fix bug in EMD vertical where sample combination was malformed (PR #113)

* Fix lisi inconsistent naming (PR #117) for issue #116.

* Fix bug in perfect integration where if batch is str (not int), it only returns control samples (PR #119).

* Fix bug in batchadjust needing "Batch_" in the sample names for non-control samples (PR #119). 

* Fix bug in cytonorm to mid where recompute was set to FALSE. It is now set to TRUE (PR #119).

* Remove transpose in harmonypy as new updates to harmonypy no longer need the transpose (PR #119). 

* Fix bug in get_obs_var_for_integrated to handle the cases where batch column in obs is str
  and thus can't be directly overriden (new values given by get_donor_batch_map is int) (PR #119).