library(FlowSOM)
library(flowCore)

#' Read the dataset-specific parameters used by the metric:
#' lineage markers, groups, and batches.
#'
#' @param adata AnnData object storing unintegrated data.
#' @return list with `lineage_markers`, `group_a`, `group_b`, `batch_1` and `batch_2`
#' 
get_dataset_parameters <- function(adata) {
  groups <- sort(as.character(unique(adata$obs$group)))
  batches <- sort(unique(adata$obs$batch))

  list(
    lineage_markers = get_lineage_markers(adata),
    group_a = groups[1],
    group_b = groups[2],
    batch_1 = batches[1],
    batch_2 = batches[2]
  )
}


#' Compute MEM scores for a set of clusters.
#'
#' Reimplementation of MEM_function.R in https://github.com/cytolab/cytoMEM.
#' The original implementation is not modular enough to extract what just
#' what we need. 
#' See the test scripts in `tests` folder to verify that this implementation 
#' produces the same results as the original.
#' The implementation here assumes that the "reference" used to compute MEM
#' for a given cluster is EVERY CELLS that are not in the cluster.
#' For example, for cluster 1, the reference is every cell that is NOT in cluster 1.
#' 
#' There are different ways of specifying the reference. It can either be cells
#' that are not in the cluster, it can be a specific set of cells, or 0 
#' (zero.ref in original mem implementation).
#' Using the default of "not in the cluster" probably makes the most sense for us
#' because the score then reflects the score of a marker relative to others
#' not in the cluster (i.e., how much higher or lower it is compared to others).
#' If using zero.ref, the score is then the absolute median expression of the cluster
#' which may look good on the surface, but is more prone to be impacted by
#' batch effects if differences in batches are mainly shifting in marker expression.
#' However, this is not saying the default is not impacted by batch effects, especially
#' the complex one, but we at least eliminate one simpler case.
#' 
#' @param expr_matrix matrix of cell by marker expression
#' @param cluster_id_per_cell character vector indicating the cluster id of 
#'   each row of `expr_matrix`
#' @param iqr_thresh numeric floor applied to per-cluster IQRs to avoid
#'   divide-by-zero (default 0.5, matching the original MEM implementation).
#'   The way this is used is, any cluster (or reference) with IQR < iqr_thresh 
#'   is floored to iqr_thresh before being used in the MEM calculation.
#' 
#' @return mem score matrix of clusters x markers where each row is a cluster
#' and each column is a marker.
#' 
compute_mem_matrix <- function(
  expr_matrix, 
  cluster_id_per_cell, 
  iqr_thresh = 0.5
) {
  cluster_id_per_cell <- as.character(cluster_id_per_cell)
  clusters <- sort(unique(cluster_id_per_cell))
  markers <- colnames(expr_matrix)

  # initialise the matrix for storing intermediate calculations we need.
  # mag stands for "magnitude", which in this case is the absolute median 
  # expression of a marker for a cluster or reference (all other cells not in the cluster).
  mag_cluster <- matrix(
    nrow = length(clusters), ncol = length(markers),
    dimnames = list(clusters, markers)
  )
  iqr_cluster <- mag_cluster
  
  # here, each row stores mag of a marker computed using all cells that are 
  # NOT in the cluster (the reference).
  mag_ref <- mag_cluster
  iqr_ref <- mag_cluster

  for (cluster in clusters) {
    cell_in_cluster <- cluster_id_per_cell == cluster

    # IQR here is the difference the third and the first quartile.
    mag_cluster[cluster, ] <- abs(apply(expr_matrix[cell_in_cluster, ], 2, median, na.rm = TRUE))
    iqr_cluster[cluster, ] <- apply(expr_matrix[cell_in_cluster, ], 2, IQR, na.rm = TRUE)

    mag_ref[cluster, ] <- abs(apply(expr_matrix[!cell_in_cluster, ], 2, median, na.rm = TRUE))
    iqr_ref[cluster, ] <- apply(expr_matrix[!cell_in_cluster, ], 2, IQR, na.rm = TRUE)
  }

  # Floor the IQRs so a cluster with zero IQR does not result in error when we
  # when we use it as a denominator in the MEM calculation. 
  # Note, the flooring is in mem implementation, see line 157 of 
  # https://github.com/cytolab/MEM-examples/blob/master/R/MEM_function.R.

  iqr_cluster[iqr_cluster < iqr_thresh] <- iqr_thresh
  iqr_ref[iqr_ref < iqr_thresh] <- iqr_thresh

  # MEM original implementation has this iqr_ratio term 
  # (stored as IQRcomp in line 155-157 of https://github.com/cytolab/MEM-examples/blob/master/R/MEM_function.R)
  # which i think act like a penalty for how spread out the IQR of the cluster is.
  # I think, if the cluster is tighter than the reference, 
  # the term will be positive because the ratio will be > 1 and subtracting by 1
  # will make the ratio remain positive. 
  # If the cluster is more spread out than the reference, the ratio will be
  # negative because the ratio will be < 1 and subtracting by 1 will make it negative.

  iqr_ratio <- (iqr_ref / iqr_cluster) - 1
  iqr_ratio[iqr_ratio < 0] <- 0

  mag_diff <- mag_cluster - mag_ref
  mem_matrix <- abs(mag_diff) + iqr_ratio

  # Markers that are lower in the cluster than in the reference get a negative score.
  # See line 168 of https://github.com/cytolab/MEM-examples/blob/master/R/MEM_function.R
  mem_matrix[mag_diff < 0] <- -mem_matrix[mag_diff < 0]

  # Following the reference implementation if scale.matrix == linear.
  # See line 176 of https://github.com/cytolab/MEM-examples/blob/master/R/MEM_function.R
  mem_matrix = (mem_matrix / max(abs(mem_matrix))) * 10

  rownames(mem_matrix) <- clusters
  colnames(mem_matrix) <- markers

  mem_matrix
}


#' Cluster and compute MEM scores
#' 
#' Cluster data with FlowSOM and calculate the MEM scores of each cluster.
#'
#' @param adata AnnData object
#' @param layer_name layer to cluster and score on
#' @param lineage_markers markers to cluster and score the clusters on
#' @param fsom_param list of FlowSOM settings, see get_flowsom_parameters()
#' @param label_prefix string prepended to the metacluster number.
#'   Set this to a label such that clusters from different batches or splits
#'   can be distinguished from one another. For example, if the clusters are 
#'   from batch 1, set label_prefix = "batch1_c". 
#' @param return_fsom_object whether to return the FlowSOM object as well (default FALSE)
#' @return list with:
#'  `cluster_id_per_cell` (cluster id of each cell, in the order
#'   the cells appear in `adata`),
#'  `mem` (cluster x marker MEM score matrix),
#'  and optionally `fsom_obj` (the FlowSOM object, if `return_fsom_object` is TRUE)
#' 
#' 
cluster_and_calculate_mem <- function(
  adata,
  layer_name,
  lineage_markers,
  fsom_param,
  label_prefix,
  return_fsom_object = FALSE
) {
  
  expr <- adata$layers[[layer_name]]
  colnames(expr) <- rownames(adata$var)
  expr <- expr[, lineage_markers]

  cat("Clustering with FlowSOM\n")

  fsom <- FlowSOM::FlowSOM(
    expr,
    colsToUse = colnames(expr),
    nClus = fsom_param$n_clusters,
    xdim = fsom_param$xdim,
    ydim = fsom_param$ydim,
    seed = 42,
    silent = FALSE
  )
  cluster_id_per_cell <- paste0(label_prefix, FlowSOM::GetMetaclusters(fsom))

  cat("Calculating MEM scores\n")
  mem_score <- compute_mem_matrix(
    expr_matrix = expr,
    cluster_id_per_cell = cluster_id_per_cell
  )

  res <- list(
    cluster_id_per_cell = cluster_id_per_cell,
    mem = mem_score
  )
  if (return_fsom_object) {
    res$fsom_obj <- fsom
  }
  res
}


#' For each sample, compute the proportion of its cells belonging to each
#' cluster.
#'
#' @param cluster_id_per_cell character vector, cluster id of each cell
#' @param sample_id_per_cell character vector, sample id of each cell
#' @return data frame with columns sample, cluster, proportion.
#'   Every sample x cluster combination is included (0 if absent), so every
#'   sample contributes a value for every cluster. The denominator is all of
#'   that sample's cells, including any that fall in clusters that are never
#'   tested.
#' 
compute_sample_proportions <- function(cluster_id_per_cell, sample_id_per_cell) {
  cluster_id_per_cell <- as.character(cluster_id_per_cell)
  sample_id_per_cell <- as.character(sample_id_per_cell)

  # row = one sample, col = one cluster, value = number of cells in that sample x cluster
  n_cells_per_sample_cluster <- table(sample_id_per_cell, cluster_id_per_cell)
  n_cells_per_sample <- rowSums(n_cells_per_sample_cluster)
  prop_cells_per_sample_cluster <- n_cells_per_sample_cluster / n_cells_per_sample

  # convert 2d matrix to long data frame, so every sample x cluster combination is a row
  prop_cells_per_sample_cluster <- as.data.frame(
    as.table(prop_cells_per_sample_cluster), stringsAsFactors = FALSE
  )
  colnames(prop_cells_per_sample_cluster) <- c("sample", "cluster", "proportion")

  prop_cells_per_sample_cluster
}


#' Test whether each cluster's per-sample proportion differs between the two
#' groups, using a Wilcoxon rank-sum test per cluster.
#'
#' @param cluster_id_per_cell character vector, cluster id of each cell
#' @param sample_id_per_cell character vector, sample id of each cell
#' @param group_id_per_cell character vector, biological group of each cell
#' @param group_a,group_b the two group values being compared
#' @param significance_threshold p-value cutoff
#' @param clusters_to_test which cluster ids to test. Defaults to every cluster
#'   present; pass a subset to keep the other clusters in the proportion
#'   denominator without testing them.
#' @return list with `pvals` (named numeric vector, NA for clusters with no
#'   cells at all) and `significant` (cluster ids whose p-value is <= the
#'   threshold)
#' 
test_differential_abundance <- function(
  cluster_id_per_cell,
  sample_id_per_cell,
  group_id_per_cell,
  group_a,
  group_b,
  significance_threshold,
  clusters_to_test = sort(unique(as.character(cluster_id_per_cell)))
) {
  sample_prop_per_cluster <- compute_sample_proportions(
    cluster_id_per_cell = cluster_id_per_cell,
    sample_id_per_cell = sample_id_per_cell
  )

  # to attach group labels to the sample id in sample_prop_per_cluster.
  # Wilcoxon test needs the group id to do the test.
  sample_group_map <- unique(data.frame(
    sample = as.character(sample_id_per_cell),
    group = as.character(group_id_per_cell),
    stringsAsFactors = FALSE
  ))

  # this is just so we can get the group label for each sample
  sample_prop_per_cluster <- merge(
    sample_prop_per_cluster, 
    sample_group_map, 
    by = "sample"
  )

  cat("Testing differential abundance with Wilcoxon rank-sum test\n")
  pvals <- sapply(clusters_to_test, function(cluster) {
    vals_a <- sample_prop_per_cluster$proportion[
      sample_prop_per_cluster$cluster == cluster & sample_prop_per_cluster$group == group_a
    ]
    vals_b <- sample_prop_per_cluster$proportion[
      sample_prop_per_cluster$cluster == cluster & sample_prop_per_cluster$group == group_b
    ]

    wilcox_test_res <- suppressWarnings(wilcox.test(vals_a, vals_b, alternative = "two.sided"))
    wilcox_test_res$p.value
  })

  list(
    pvals = pvals,
    significant = names(pvals)[!is.na(pvals) & pvals <= significance_threshold],
    clusters_tested = clusters_to_test
  )
}


#' MEM RMSD similarity between every pair of MEM score vectors from two
#' clusters x markers matrices.
#' Reimplementation of the MEM RMSD original implementation available from:
#' in https://github.com/cytolab/cytoMEM (MEM_RMSD.R)
#'
#' @param mem_a matrix, clusters x markers
#' @param mem_b matrix, clusters x markers
#' @return matrix, nrow(mem_a) x nrow(mem_b), similarity over the markers common to both
#' 
compute_mem_rmsd <- function(mem_a, mem_b) {
  common_markers <- intersect(colnames(mem_a), colnames(mem_b))
  mem_a <- mem_a[, common_markers, drop = FALSE]
  mem_b <- mem_b[, common_markers, drop = FALSE]

  rmsd <- matrix(
    nrow = nrow(mem_a), ncol = nrow(mem_b),
    dimnames = list(rownames(mem_a), rownames(mem_b))
  )

  for (i in seq_len(nrow(mem_a))) {
    for (j in seq_len(nrow(mem_b))) {
      rmsd[i, j] <- sqrt(mean((mem_a[i, ] - mem_b[j, ])^2))
    }
  }

  100 - (rmsd / 20 * 100)
}


#' Find a match for each cluster using MEM RMSD similarity.
#'
#' For every cluster in mem_a, find the cluster in mem_b whose MEM vector is the best match.
#' Then reverse it in that for every cluster in mem_b, find the cluster in mem_a whose MEM vector is the best match.
#' 
#' The best match here is defined as the cluster with the highest MEM RMSD similarity score.
#' 
#' @param mem_a cluster x marker MEM matrix for a clustering.
#' @param mem_b cluster x marker MEM matrix for a clustering
#' 
#' @return list with `pairs` (data frame with one row per cluster in mem_a:
#'   cluster_a, best_match_b, best_match_b_matching_a, is_mutual) and
#'   `similarity` (the full batch A x batch B similarity matrix)
#' 
find_cluster_match <- function(mem_a, mem_b) {
  # rows are mem_a clusters, cols are mem_b clusters, values are MEM RMSD similarity
  similarity_matrix <- compute_mem_rmsd(mem_a = mem_a, mem_b = mem_b)
  clusters_a <- rownames(similarity_matrix)
  clusters_b <- colnames(similarity_matrix)

  # for each cluster a, find the best match in cluster B
  clusters_a_best_match <- sapply(clusters_a, function(cluster_a) {
    clusters_b[which.max(similarity_matrix[cluster_a, ])]
  })
  # do the reverse for each cluster b, find the best match in cluster A
  clusters_b_best_match <- sapply(clusters_b, function(cluster_b) {
    clusters_a[which.max(similarity_matrix[, cluster_b])]
  })

  # checking whether for each cluster in A, its best match in B is also best matched back to it in A.
  is_cluster_a_best_match_mutual <- sapply(clusters_a, function(cluster_a) {
    matching_cluster_b <- clusters_a_best_match[cluster_a]

    # what is the best matching cluster A for that cluster B (matching_cluster_b)?
    best_match_a_for_b <- clusters_b_best_match[matching_cluster_b]

    # are they the same
    is_mutual <- best_match_a_for_b == cluster_a
    
    return(is_mutual)
  })

  # what we return is just a dataframe
  cluster_a_matches <- data.frame(
    cluster_a = clusters_a,
    best_match_b = clusters_a_best_match,
    # this is for subsequent sanity check to see if the best match for B
    # is the same as the cluster A column we started with.
    best_match_b_matching_a = clusters_b_best_match[clusters_a_best_match],
    is_mutual = is_cluster_a_best_match_mutual,
    stringsAsFactors = FALSE
  )

  list(cluster_a_matches = cluster_a_matches, similarity = similarity_matrix)
}


#' For every row of `mem_query`, name the row of `mem_reference` whose MEM
#' vector it is closest to.
#'
#' @param mem_query matrix, clusters x markers
#' @param mem_reference matrix, clusters x markers
#' @return list with `nearest` (named character vector, query cluster ->
#'   reference cluster) and `similarity` (the full query x reference MEM RMSD
#'   similarity matrix, kept for diagnostics)
#' 
find_nearest_by_mem <- function(mem_query, mem_reference) {
  similarity <- compute_mem_rmsd(mem_a = mem_query, mem_b = mem_reference)
  ref_clusters <- colnames(similarity)
  query_clusters <- rownames(similarity)
  
  nearest_ref_cluster_for_query <- ref_clusters[apply(similarity, 1, which.max)]

  list(
    nearest = setNames(nearest_ref_cluster_for_query, query_clusters),
    similarity = similarity
  )
}


#' Cluster one integrated split, map its clusters back onto the batch-1
#' clusters that were kept as differentially abundant, and re-test whether
#' those clusters are still differentially abundant after integration.
#'
#' Clusters are matched against the *whole* batch-1 clustering rather than
#' only the retained ones, so a cluster whose nearest counterpart was never
#' retained drops out instead of being forced onto a retained cluster. Its
#' cells still count toward each sample's total, keeping the proportions on
#' the same scale as the unintegrated test.
#'
#' @param adata AnnData object, one integrated split, controls already removed
#' @param dataset_param list of dataset settings, see get_dataset_parameters()
#' @param fsom_param list of FlowSOM settings, see get_flowsom_parameters()
#' @param b1_mem cluster x marker MEM matrix for the unintegrated batch 1
#' @param retained_cluster_b1 character vector of the batch 1 cluster ids that
#'   survived the mutual-match and differential-abundance filters
#' @param significance_threshold p-value cutoff
#' @param label_prefix string prepended to this split's cluster numbers
#' @return list with `da_results` (the test_differential_abundance() output,
#'   keyed by batch 1 cluster), `cluster_id_per_cell` (this split's own cluster
#'   ids), `mapped_b1_cluster_per_cell` (the batch 1 cluster each cell was
#'   pooled under, or "unmatched") and `similarity_to_batch1`
#' 
process_integrated_split <- function(
  adata,
  dataset_param,
  fsom_param,
  b1_mem,
  retained_cluster_b1,
  significance_threshold,
  label_prefix
) {
  clustering <- cluster_and_calculate_mem(
    adata = adata,
    layer_name = "integrated",
    lineage_markers = dataset_param$lineage_markers,
    fsom_param = fsom_param,
    label_prefix = label_prefix
  )

  # for each cluster, find the nearest cluster in batch 1
  # keep it simple, ignore b2 because if cluster 1 b1 nearest is cluster 2 in b2 (when b1 is compared against b2),
  # even if cluster 10 in integrated nearest in b1 is cluster 1, it may not necessarily
  # mean that its nearest cluster in b2 is cluster 2. 
  # so keep it simple, just compare with b1.
  # note, mem_rmsd_for_clusters is a list. 
  
  mem_rmsd_for_clusters <- find_nearest_by_mem(
    mem_query = clustering$mem, 
    mem_reference = b1_mem
  )

  # keep only those which b1 clusters are found to be DA and are mutually matched
  nearest_b1_clusters <- mem_rmsd_for_clusters$nearest
  # note, names are the integrated split clusters, and values are the nearest b1 clusters.
  nearest_b1_clusters <- nearest_b1_clusters[nearest_b1_clusters %in% retained_cluster_b1]
 
  # prepare the dataframe for DA
  cell_metadata_for_da <- data.frame(
    cluster_id_per_cell = clustering$cluster_id_per_cell,
    sample_id_per_cell = adata$obs$sample,
    group_id_per_cell = adata$obs$group,
    stringsAsFactors = FALSE
  )
  # add the mapped b1
  cell_metadata_for_da$mapped_b1_cluster <- nearest_b1_clusters[cell_metadata_for_da$cluster_id_per_cell]
  # Cells in a cluster that did not map to a retained b1 cluster are not tested,
  # but they still belong to their sample, so they need a label. Left as NA they
  # would be dropped by table(). Hence call them "unmatched".
  cell_metadata_for_da$mapped_b1_cluster[is.na(cell_metadata_for_da$mapped_b1_cluster)] <- "unmatched"

  # Cells are pooled by the b1 cluster they mapped to, so that is what gets
  # tested - one test per retained b1 cluster, not one per split cluster.
  abundance <- test_differential_abundance(
    cluster_id_per_cell = cell_metadata_for_da$mapped_b1_cluster,
    sample_id_per_cell = cell_metadata_for_da$sample_id_per_cell,
    group_id_per_cell = cell_metadata_for_da$group_id_per_cell,
    group_a = dataset_param$group_a,
    group_b = dataset_param$group_b,
    significance_threshold = significance_threshold,
    clusters_to_test = setdiff(unique(cell_metadata_for_da$mapped_b1_cluster), "unmatched")
  )

  list(
    da_results = abundance,
    cluster_id_per_cell = clustering$cluster_id_per_cell,
    mapped_b1_cluster = mem_rmsd_for_clusters$nearest,
    similarity_to_batch1 = mem_rmsd_for_clusters$similarity
  )
}
