# NOTE: These helper functions are ports of the original Python functions in 'helper_functions.py'

library(dplyr)
requireNamespace("anndataR", quietly = TRUE)

#' Adds annotations (.var and .obs) from the unintegrated data to the
#' integrated dataset.
#' In the case of the control method "perfect_integration",
#' the function will fetch the batch label from the unintegrated data
#' based on the split.
#' i.e., if in split 1, donor 3-5 is from batch 2, then the batch label for that split
#' will be changed from batch 1 to batch 2.
#' Note: this implementation slightly differs from the `get_obs_var_for_integrated` python version as it has to be run each time for each split.
#'
#' @param i_adata AnnData object, integrated data
#' @param u_adata AnnData object, unintegrated dataset
#' @param split_id numeric, split id of the integrated data
#' @return AnnData object with .var and .obs added
#'
get_obs_var_for_integrated <- function(i_adata, u_adata, split_id) {

    i_adata$obs <- u_adata$obs[i_adata$obs_names, ]
    i_adata$var <- u_adata$var[i_adata$var_names, ]

    # if integrated data came from perfect integration, change the batch labels of the samples
    # everything is from batch 1, but some samples need to be labelled to come from batch 2
    if (i_adata$uns["method_id"] == "perfect_integration") {
        cat(
            "Control method 'perfect_integration' detected. Changing batch labels.\n"
        )
        cat("Computing new batch labels\n")
        # mutate is needed as donors that are used for controls, we won't have the mapping
        i_adata_new_batch_labels <- get_batch_label_perfect_integration(
            u_adata = u_adata,
            i_adata = i_adata,
            split_id = split_id
        )

        # safeguard
        if (! all(i_adata_new_batch_labels$donor == i_adata$obs$donor)) {
            stop(
                "Donor labels do not match between new batch labels and integrated data. This should not happen!"
            )
        }

        cat("Attaching new batch labels\n")
        i_adata$obs$batch <- i_adata_new_batch_labels$new_batch_label
    }

    return(i_adata)
}

#' Helper function to get the batch label for perfect integration.
#' First, get donor batch map for a given split.
#' Then apply the map to the integrated data.
#'
#' @param u_adata AnnData object, unintegrated dataset.
#' @param i_adata AnnData object, integrated data.
#' @param split_id numeric, split id of the integrated data.
#'
#' @return a dataframe with donor and new batch label
#'
get_batch_label_perfect_integration <- function(u_adata, i_adata, split_id) {
    # this return which batch sample we used for a donor for a given split
    actual_donor_batch_map <- unique(
        u_adata$obs[(u_adata$obs$split == split_id), c("donor", "batch")]
    )
    # mutate is needed as donors that are used for controls, won't have batch_new as 
    # the split id is 0
    i_adata_new_batch_labels <- i_adata$obs[, c("donor", "batch")] %>%
        left_join(actual_donor_batch_map, by = "donor", suffix = c("_old", "_new")) %>%
        mutate(new_batch_label = ifelse(is.na(batch_new), batch_old, batch_new)) %>%
        select(donor, new_batch_label)

    return(i_adata_new_batch_labels)
}

#' Subsets the anndata object to remove the control cells.
#' Control cells are all the entries in adata.obs where is_control != 0.
#'
#' @param adata AnnData object
#' @return AnnData object with cells from control samples removed
subset_nocontrols <- function(adata) {
    if (!"is_control" %in% colnames(adata$obs)) {
        stop("The column 'is_control' is not present in the adata object.")
    }

    # Subset the adata to remove cells where is_control != 0
    adata[adata$obs$is_control == 0, ]$copy()
}


#' Subset AnnData Object to Single Control
#'
#' @description
#' Extracts no control samples and just one control sample from an AnnData object based on the specified index.
#'
#' @param adata An AnnData object containing multiple samples with control groups.
#' @param which_control Integer specifying which control sample to extract (default: 1).
#'
#' @return
#' An AnnData object subset to contain only the specified control sample.
#'
#' @examples
#' \dontrun{
#' control_sample <- subset_onecontrol(adata, which_control = 1)
#' }
#'
subset_onecontrol <- function(adata, which_control = 1) {
    if (!"is_control" %in% colnames(adata$obs)) {
        stop("The column 'is_control' is not present in the adata object.")
    }

    # Subset the adata to keep cells where is_control == which_control
    adata[adata$obs$is_control %in% c(which_control, 0), ]$copy()
}

#' Subsets the anndata object to only include markers that need to be
#' (or have been) corrected.
#' These markers are all the entries in adata.var where to_correct == TRUE.
#'
#' @param adata AnnData object
#' @return AnnData object with only the markers to correct
subset_markers_tocorrect <- function(adata) {
    adata[, adata$var$to_correct]
}

#' Subsets the anndata object to remove all cells where the marker is not
#' labelled. This is determined by the column "cell_type" in adata.obs.
#' Particularly useful when dealing with cell type specific metrics.
#'
#' @param adata AnnData object
#' @return AnnData object with only the labeled cells
remove_unlabelled <- function(adata) {
    if (!"cell_type" %in% colnames(adata$obs)) {
        stop("The column 'cell_type' is not present in the adata object.")
    }

    # Convert to lowercase and filter out "unlabelled" and "unlabeled"
    is_unlabelled <- tolower(adata$obs$cell_type) %in%
        c("unlabelled", "unlabeled")
    adata[!is_unlabelled, ]
}

#' Subsets the anndata object in a stratified manner
#' with 'cell type' and 'sample' as strata.
#'
#' @param adata AnnData object
#' @param frac numeric, fraction of cells to keep for each cell type
#' @param seed numeric, seed for reproducibility
#' @param anndatar logical, whether the input is anndataR object or not
#' @return AnnData object with only the markers to correct
subset_by_celltype <- function(adata, frac = 0.5, seed = 1, anndatar = TRUE) {
  set.seed(seed)
  
  obs <- adata$obs
  obs$cell_id <- rownames(obs)
  obs$.row <- seq_len(nrow(obs))   # original order
  
  keep_ids <- obs %>%
    group_by(cell_type, sample) %>%
    slice_sample(prop = frac) %>%
    ungroup() %>%
    arrange(.row) %>%     # restore original order
    pull(cell_id)
  
  if (anndatar == TRUE){
    keep_idx <- match(keep_ids, adata$obs_names)
    
    adata_sub <- anndataR::AnnData(
      X   = NULL,
      obs = adata$obs[keep_idx, , drop = FALSE],
      var = adata$var,
      uns = adata$uns,
      layers = list(
        "integrated" = adata$layers$integrated[keep_idx, , drop = FALSE]
      )
    )
  } else{
    adata_sub <- adata[keep_ids, ]
  }
}

