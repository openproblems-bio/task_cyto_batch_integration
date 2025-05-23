library(flowCore)
library(anndata)
library(FlowSOM)

#' @title Create a flowframe and a cell type annotation list from anndata object
#'
#' @param adata Anndata object
#' @param donor_name a field present in adata$obs$donor
#' @param layer name of the layer to which expression values are fetched for the creation of the FlowSOM tree
get_ff_annotations <- function(adata, donor_name, layer_name) {
  if (!is.element(donor_name, unique(adata$obs$donor))) {
    stop(paste0(donor_name, " not found in the anndata object"))
  }
  #extract expression matrix + annotations from donor-specific cells
  expression_matrix <- adata$layers[[layer_name]][
    adata$obs$donor == donor_name,
  ]
  ct_annotations <- adata$obs$cell_type[adata$obs$donor == donor_name]
  # Create flowframe
  colnames(expression_matrix) <- adata$var_names
  ff <- flowFrame(expression_matrix)
  list(
    flowframe = ff,
    ct_annotations = ct_annotations
  )
}

#' @title Create a dataframe with proportions of cell types for each metacluster
#'
#' @param vector_meta  factor containing the metacluster annotations for each cell
#' @param vector_celltype factor (or any vector) containing the metacluster annotations for each (integrated) cell
#'
#' @note it is important that vector_meta is a factor
compute_meta_pct <- function(vector_meta, vector_celltype) {
  tab <- table(vector_meta, vector_celltype)
  pct_tab <- prop.table(tab, margin = 1) * 100
  pct_tab <- round(pct_tab, 2)
  pct_tab[is.na(pct_tab)] <- 0 #Nans as zeroes
  pct_tab
}

#' @title Get a vector with metacluster weights for the metric 'compute_fs_mapping_similarity'.
#' @description The outuput is a vector with length == n. of clusters where each element represent the proportions of that metacluster.
#' The proportions take into account both the integrated cells and the validations cell together.
#' @param vector_meta_int a factor containing the metacluster annotations for each (integrated) cell
#' @param vector_meta_val a factor containing the metacluster annotations for each (validation) cell
#'
#' @note Metaclusters which are empty will have a weight == 0
compute_meta_weights <- function(vector_meta_int, vector_meta_val) {
  vector_meta_tot <- c(vector_meta_int, vector_meta_val)
  meta_total <- table(vector_meta_tot, exclude = NULL)
  proportions(meta_total)
}

#' @title Compute flowsom mapping similarity score for a sample pair.
#'
#' @description
#' The function computes a similarity score for 2 FlowSOM tree given in input (along with cell type annotations).
#' It returns a series of object including :
#' - FlowSOM trees used
#' - dissimilarity score: a score from 0 to 100 where 0 means perfect overlap between
#' percentages of cell types within a metacluster for a sample pair (integratied & validation data);
#' and 100 means absence overlap.
#' - similarity score: 100 - dissimilarity score. A score of 100 means that the 2 FlowSOM trees
#' are perfectly overlapping in terms of cell type percentages
#' - vector with differences in percentages of cellt types for each metacluster
#' - vector with weight for each metacluster. The weight is the proportion of cell types in a metacluster out of the total number cells
#' (validation + integrated cells)
#'
#' @param fs_tree_int FlowSOM tree object of integrated data via NewData()
#' @param ct_ann_int vector containing the metacluster annotations for each (integrated) cell
#' @param fs_tree_val FlowSOM tree object created on validation data
#' @param ct_ann_val vector containing the metacluster annotations for each (validation) cell
#'
#' @note When computing the absolute difference matrix, metaclusters that have no differences in cell proportions
#' or metaclusters which are empty in both validation and integrated data(== SOM node empty)
#' will result in all zeroes rows. This has no effect on the final metric
compute_fs_mapping_similarity <- function(
  fs_tree_int,
  ct_ann_int,
  fs_tree_val,
  ct_ann_val
) {
  # Get metacluster-level annotations for each sample
  meta_int <- factor(
    paste0("meta ", GetMetaclusters(fs_tree_int)),
    levels = paste0("meta ", levels(fs_tree_int$metaclustering))
  )
  meta_val <- factor(
    paste0("meta ", GetMetaclusters(fs_tree_val)),
    levels = paste0("meta ", levels(fs_tree_val$metaclustering))
  )

  # Get metacluster-level proportions for the paired sample
  meta_weights <- compute_meta_weights(meta_int, meta_val)

  # Get metacluster x cell type proportions
  table_metaxcell_int <- compute_meta_pct(meta_int, ct_ann_int)
  table_metaxcell_val <- compute_meta_pct(meta_val, ct_ann_val)

  table_metaxcell_absdiff <- abs(table_metaxcell_val - table_metaxcell_int)
  # # Debug:
  # heatmap(table_metaxcell_int, Rowv = NA, Colv = NA)
  # heatmap(table_metaxcell_val, Rowv = NA, Colv = NA)
  # heatmap(table_metaxcell_absdiff, Rowv = NA, Colv = NA)

  #Differences in cell type proportions of metaclusters
  meta_differences <- rowSums(table_metaxcell_absdiff)

  #Final metric
  FSOM_dissimilarity <- sum(meta_weights * meta_differences)
  FSOM_similarity <- 100 - FSOM_dissimilarity

  list(
    similarity = FSOM_similarity,
    dissimilarity = FSOM_dissimilarity,
    tree_integrated = fs_tree_int,
    tree_validation = fs_tree_val,
    differences = meta_differences,
    weights = meta_weights
  )
}
