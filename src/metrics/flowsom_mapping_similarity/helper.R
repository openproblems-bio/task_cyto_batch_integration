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

#' @title Create a dataframe with proportions of cell types for each cluster
#'
#' @param vector_clust  factor containing the cluster annotations for each cell
#' @param vector_celltype factor (or any vector) containing the cluster annotations for each (integrated) cell
#'
#' @note it is important that vector_clust is a factor

compute_clust_pct <- function(vector_clust, vector_celltype){
  tab <- table(vector_clust, vector_celltype)
  pct_tab <- prop.table(tab, margin = 1) * 100
  pct_tab <- round(pct_tab, 2)
  pct_tab[is.na(pct_tab)] <- 0 #Nans as zeroes
  return(pct_tab)
}

#' @title Get a vector with cluster weights for the metric 'compute_fs_mapping_similarity'.
#' @description The outuput is a vector with length == n. of clusters where each element represent the proportions of that cluster.
#' The proportions take into account both the integrated cells and the validations cell together.
#' @param vector_clust_int a factor containing the cluster annotations for each (integrated) cell
#' @param vector_clust_val a factor containing the cluster annotations for each (validation) cell
#'
#' @note clusters which are empty will have a weight == 0
compute_clust_weights <- function(vector_clust_int, vector_clust_val){
  vector_clust_tot <- c(vector_clust_int,vector_clust_val)
  clust_total <- table(vector_clust_tot, exclude=NULL)
  clust_prop <- proportions(clust_total)
  return(clust_prop)
}

#' @title Compute flowsom mapping similarity score for a sample pair.
#'
#' @description
#' The function computes a similarity score for 2 FlowSOM tree given in input (along with cell type annotations).
#' It returns a series of object including :
#' - FlowSOM trees used
#' - dissimilarity score: a score from 0 to 100 where 0 means perfect overlap between
#' percentages of cell types within a cluster for a sample pair (integratied & validation data);
#' and 100 means absence overlap.
#' - similarity score: 100 - dissimilarity score. A score of 100 means that the 2 FlowSOM trees
#' are perfectly overlapping in terms of cell type percentages
#' - vector with differences in percentages of cellt types for each cluster
#' - vector with weight for each cluster. The weight is the proportion of cell types in a cluster out of the total number cells
#' (validation + integrated cells)
#'
#' @param fs_tree_int FlowSOM tree object of integrated data via NewData()
#' @param ct_ann_int vector containing the cluster annotations for each (integrated) cell
#' @param fs_tree_val FlowSOM tree object created on validation data
#' @param ct_ann_val vector containing the cluster annotations for each (validation) cell
#'
#' @note When computing the absolute difference matrix, clusters that have no differences in cell proportions
#' or clusters which are empty in both validation and integrated data(== SOM node empty)
#' will result in all zeroes rows. This has no effect on the final metric
compute_fs_mapping_similarity <- function(
  fs_tree_int,
  ct_ann_int,
  fs_tree_val,
  ct_ann_val
) {

  #Get cluster-level annotations for each sample
  clust_levels <- levels(as.factor(c(GetClusters(fs_tree_val),GetClusters(fs_tree_int))))
  clust_int <- factor(
    paste0("clust ", GetClusters(fs_tree_int)), levels = paste0("clust ",clust_levels)
  )
  clust_val <- factor(
    paste0("clust ", GetClusters(fs_tree_val)), levels = paste0("clust ",clust_levels)
  )

  #Get cluster-level proportions for the paired sample
  clust_weights <- compute_clust_weights(clust_int,clust_val)
  
  #Get cluster x cell type proportions
  table_clustxcell_int <- compute_clust_pct(clust_int, ct_ann_int)
  table_clustxcell_val <- compute_clust_pct(clust_val, ct_ann_val)
  
  table_clustxcell_absdiff <- abs(table_clustxcell_val - table_clustxcell_int)
  # # Debug:
  # heatmap(table_clustxcell_int, Rowv = NA, Colv = NA)
  # heatmap(table_clustxcell_val, Rowv = NA, Colv = NA)
  # heatmap(table_clustxcell_absdiff, Rowv = NA, Colv = NA)

  #Differences in cell type proportions of clusters
  clust_differences <- rowSums(table_clustxcell_absdiff)

  #Final metric
  FSOM_dissimilarity <- sum(clust_weights*clust_differences)
  FSOM_similarity <- 100 - FSOM_dissimilarity

  list(
    similarity = FSOM_similarity,
    dissimilarity = FSOM_dissimilarity,
    tree_integrated = fs_tree_int,
    tree_validation = fs_tree_val,
    differences = clust_differences,
    weights = clust_weights
  )
}
