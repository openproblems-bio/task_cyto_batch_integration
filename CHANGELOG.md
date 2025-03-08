# task_cyto_batch_integration x.y.z

## BREAKING CHANGES

<!-- * Restructured `src` directory (PR #3). -->

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

* Added `metrics/n_inconsistent_peaks` + some new helper functions (PR #31)

* Added `methods/cycombine_nocontrol` (PR #35).

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

* Change emd to emd_mean and emd_max (PR #27).

* Update output anndata for methods to return all vars - corrected or not (PR #28).

* Add positive control (PR #30).

## MINOR CHANGES

* Enabled unit tests (PR #2).

* Added integrated test resource (PR #5).

* Updated file description in yaml file (PR #15).

* Updated scripts to enable running benchmark on seqera (PR #29).

## BUGFIXES

