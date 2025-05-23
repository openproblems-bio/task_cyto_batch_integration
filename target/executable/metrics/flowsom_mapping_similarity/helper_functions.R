# NOTE: These helper functions are ports of the original Python functions in 'helper_functions.py'

#' Adds annotations (.var and .obs) from the unintegrated dataset to the
#' integrated dataset. In the case of the control method "perfect_integration",
#' the function will fetch annotations from the validation dataset instead.
#'
#' @param i_adata AnnData object, batch-integrated dataset
#' @param v_adata AnnData object, validation dataset
#' @param u_adata AnnData object, unintegrated dataset
#' @return AnnData object with .var and .obs added
get_obs_var_for_integrated <- function(i_adata, v_adata, u_adata) {
  if (i_adata$uns$method_id == "perfect_integration_horizontal") {
    if (i_adata$n_obs != v_adata$n_obs) {
      stop(
        "The number of cells in the integrated (perfect_integration_horizontal) ",
        "and validation datasets do not match"
      )
    }
    i_adata$obs <- v_adata$obs[rownames(i_adata), , drop = FALSE]
    i_adata$var <- v_adata$var[colnames(i_adata), , drop = FALSE]
  } else if (i_adata$uns$method_id == "perfect_integration_vertical") {
    comb_adata <- anndata::concat(list(v_adata, u_adata))
    # subset to just batch 1
    # Check if 'batch' column exists
    if (!"batch" %in% colnames(comb_adata$obs)) {
      stop(
        "Column 'batch' not found in comb_adata$obs for ",
        "perfect_integration_vertical."
      )
    }
    comb_adata <- comb_adata[comb_adata$obs$batch == 1, ]

    if (i_adata$n_obs != comb_adata$n_obs) {
      stop(
        "The number of cells in the integrated (perfect_integration_vertical) ",
        "and validation + unintegrated datasets do not match."
      )
    }
    i_adata$obs <- comb_adata$obs[rownames(i_adata), , drop = FALSE]
    i_adata$var <- v_adata$var[colnames(i_adata), , drop = FALSE]
  } else {
    if (i_adata$n_obs != u_adata$n_obs) {
      stop(
        "The number of cells in the integrated and unintegrated datasets do not match"
      )
    }
    # Compare obs_names for ordering
    if (!all(rownames(i_adata) == rownames(u_adata))) {
      warning(
        "The cell ordering in the integrated and unintegrated datasets do not match"
      )
    }

    i_adata$obs <- u_adata$obs[rownames(i_adata), , drop = FALSE]
    i_adata$var <- u_adata$var[colnames(i_adata), , drop = FALSE]
  }

  i_adata
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
  adata[adata$obs$is_control == 0, ]
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
