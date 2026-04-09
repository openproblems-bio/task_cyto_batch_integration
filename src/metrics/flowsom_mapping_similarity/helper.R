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
#' The proportions refer to both splits.
#' @param vector_clust_s1 a factor containing the cluster annotations for each cell belonging to split 1.
#' @param vector_clust_s2 a factor containing the cluster annotations for each cell belonging to split 2.
#'
#' @note clusters which are empty will have a weight == 0
compute_clust_weights <- function(vector_clust_s1, vector_clust_s2){
  vector_clust_tot <- c(vector_clust_s1,vector_clust_s2)
  clust_total <- table(vector_clust_tot, exclude=NULL)
  clust_prop <- proportions(clust_total)
  return(clust_prop)
}

#' @title Compute flowsom mapping similarity score for a sample pair.
#'
#' @description
#' The function computes a similarity score for 2 FlowSOM trees built on technical replicates (from 2 splits) given in input (along with cell type annotations).
#' It returns a series of object including :
#' - FlowSOM trees used
#' - dissimilarity score: a score from 0 to 100 where 0 means perfect overlap between
#' percentages of cell types within a cluster for a sample pair (split 1 & split 2 data);
#' and 100 means absence overlap.
#' - similarity score: 100 - dissimilarity score. A score of 100 means that the 2 FlowSOM trees
#' are perfectly overlapping in terms of cell type percentages
#' - vector with differences in percentages of cellt types for each cluster
#' - vector with weight for each cluster. The weight is the proportion of cell types in a cluster out of the total number cells
#' (split 1 + split 2 cells)
#'
#' @param fs_tree_s1 FlowSOM tree object from one split
#' @param ct_ann_s1 vector containing the cluster annotations for each integrated cell belonging to `fs_tree_s1`
#' @param fs_tree_s2 FlowSOM tree object created on the other split
#' @param ct_ann_s2 vector containing the cluster annotations for each integrated cell belonging to `fs_tree_s2`
#'
#' @note When computing the absolute difference matrix, clusters that have no differences in cell proportions
#' or clusters which are empty in both data splits(== SOM node empty)
#' will result in all zeroes rows. This has no effect on the final metric
compute_fs_mapping_similarity <- function(
  fs_tree_s1,
  ct_ann_s1,
  fs_tree_s2,
  ct_ann_s2
) {

  #Get cluster-level annotations for each sample
  clust_levels <- levels(as.factor(c(GetClusters(fs_tree_s1),GetClusters(fs_tree_s2))))
  clust_s1 <- factor(
    paste0("clust ", GetClusters(fs_tree_s1)), levels = paste0("clust ",clust_levels)
  )
  clust_s2 <- factor(
    paste0("clust ", GetClusters(fs_tree_s2)), levels = paste0("clust ",clust_levels)
  )

  #Get cluster-level proportions for the paired sample
  clust_weights <- compute_clust_weights(clust_s1,clust_s2)

  #Get cluster x cell type proportions
  table_clustxcell_s1 <- compute_clust_pct(clust_s1, ct_ann_s1)
  table_clustxcell_s2 <- compute_clust_pct(clust_s2, ct_ann_s2)

  table_clustxcell_absdiff <- abs(table_clustxcell_s2 - table_clustxcell_s1)
  # # Debug:
  # heatmap(table_clustxcell_s1, Rowv = NA, Colv = NA)
  # heatmap(table_clustxcell_s2, Rowv = NA, Colv = NA)
  # heatmap(table_clustxcell_absdiff, Rowv = NA, Colv = NA)

  #Differences in cell type proportions of clusters
  clust_differences <- rowSums(table_clustxcell_absdiff)

  #Final metric
  FSOM_dissimilarity <- sum(clust_weights*clust_differences)
  FSOM_similarity <- 100 - FSOM_dissimilarity

  list(
    similarity = FSOM_similarity,
    dissimilarity = FSOM_dissimilarity,
    tree_split1 = fs_tree_s1,
    tree_split2 = fs_tree_s2,
    differences = clust_differences,
    weights = clust_weights
  )
}
