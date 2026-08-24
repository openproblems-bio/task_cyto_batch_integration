library(FlowSOM)
library(flowCore)

# Functions are declared in the order script.R reaches them. Where a function
# is only ever called by another one in here, it sits directly above its caller.


#' Read the FlowSOM settings the dataset was published with.
#'
#' Every component clusters with the same grid size and metacluster count, so
#' the numbers live in the dataset's uns rather than being hardcoded here.
#'
#' @param adata AnnData object, the unintegrated dataset
#' @return list with `n_clusters` and the SOM grid size `xdim`, `ydim`
#' 
get_flowsom_parameters <- function(adata) {
  list(
    n_clusters = adata$uns$parameter_num_clusters,
    xdim = adata$uns$parameter_som_xdim,
    ydim = adata$uns$parameter_som_ydim
  )
}


#' Read the dataset-specific settings the metric is run with: the markers to
#' cluster on, the two biological groups being compared, and the two batches.
#'
#' Batch 1 is the batch that supplies the reference clusters and batch 2 is the
#' one used to confirm them; which is which is just the sort order of the batch
#' labels, so the choice is at least stable across runs.
#'
#' @param adata AnnData object, the unintegrated dataset with controls removed
#' @return list with `lineage_markers`, `group_a`, `group_b`, `batch_1` and `batch_2`
#' 
get_dataset_parameters <- function(adata) {
  groups <- sort(as.character(unique(adata$obs$group)))
  batches <- sort(unique(adata$obs$batch))

  list(
    lineage_markers = rownames(adata$var)[adata$var$marker_type == "lineage"],
    group_a = groups[1],
    group_b = groups[2],
    batch_1 = batches[1],
    batch_2 = batches[2]
  )
}


#' Compute MEM scores for a set of clusters.
#'
#' Reimplements the MEM algorithm (Diggins et al. 2017, Nature Methods; see
#' https://github.com/cytolab/cytoMEM, MEM_function.R) as the original R
#' code is too rigid in what it produces and how it is called..
#' 
#'
#' @param expr_matrix matrix of cell by marker expression
#' @param cluster_id_per_cell character vector giving the cluster id of each row of
#'   `expr_matrix`, so it has one entry per cell and is in the same order as the
#'   rows. Every cluster it names gets a MEM score row.
#' @param iqr_thresh numeric floor applied to per-cluster IQRs to avoid
#'   divide-by-zero (default 0.5, matching the reference implementation).
#' @return mem score matrix of clusters x markers where each row is a cluster
#' and each column is a marker. Values are the MEM score of that marker in that cluster.
#' 
compute_mem_matrix <- function(expr_matrix, cluster_id_per_cell, iqr_thresh = 0.5) {
  cluster_id_per_cell <- as.character(cluster_id_per_cell)
  clusters <- sort(unique(cluster_id_per_cell))
  markers <- colnames(expr_matrix)

  # note, MAG stands for the marker's magnitude, which is the absolute median expression 
  # of a marker in the cluster or reference.
  # for this implementation, ref is essentially anything else that is not in the cluster.
  # so if you have 0-10 clusters, for cluster 1, ref will be whatever that is
  # NOT in cluster 1.
  # these variables below, are just for storing the absolute median (mag) and IQR values, 
  # not for any calculations.
  mag_cluster <- matrix(
    nrow = length(clusters), ncol = length(markers),
    dimnames = list(clusters, markers)
  )
  mag_ref <- mag_cluster
  iqr_cluster <- mag_cluster
  iqr_ref <- mag_cluster

  # MEM calls this the marker's "magnitude": the absolute median expression,
  # measured once within the cluster and once over every other cell.
  for (cluster in clusters) {
    # to get expr matrix for each cluster.
    in_cluster_mask <- cluster_id_per_cell == cluster

    mag_cluster[cluster, ] <- abs(apply(expr_matrix[in_cluster_mask, ], 2, median, na.rm = TRUE))
    iqr_cluster[cluster, ] <- apply(expr_matrix[in_cluster_mask, ], 2, IQR, na.rm = TRUE)

    mag_ref[cluster, ] <- abs(apply(expr_matrix[!in_cluster_mask, ], 2, median, na.rm = TRUE))
    iqr_ref[cluster, ] <- apply(expr_matrix[!in_cluster_mask, ], 2, IQR, na.rm = TRUE)
  }

  # Floor the IQRs so a cluster with (near-)zero spread doesn't blow up the ratio below.
  iqr_cluster[iqr_cluster < iqr_thresh] <- iqr_thresh
  iqr_ref[iqr_ref < iqr_thresh] <- iqr_thresh

  # MEM original implementation has this iqr_ratio term which i think
  # act like a penalty for how spread out the IQR of the cluster is.
  # This iqr_ratio is added after calculation the difference in MAG beftween
  # cluster and reference.
  # I think, if the cluster is very tight, the term will be large and positive.
  # Why positive? Because the iqr_ref / iqr_cluster will be > 1, so subtracting 1 will remain positive.
  # If the cluster is very spread out, ratio will be small, and if the spread in
  # cluster is larger than the spread in the reference, the ratio will be < 1, and 
  # subtracting 1 will make it negative.
  # Note, the flooring is in mem implementation, see line 157 of https://github.com/cytolab/MEM-examples/blob/master/R/MEM_function.R.

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
#' Cluster data with FlowSOM and find the MEM scores of each cluster.
#'
#' @param adata AnnData object storing just the cells to cluster.
#' @param layer_name layer to cluster and score on ("preprocessed" or "integrated")
#' @param lineage_markers markers to cluster and score the clusters on
#' @param fsom_param list of FlowSOM settings, see get_flowsom_parameters()
#' @param label_prefix string prepended to the metacluster number, so cluster
#'   names stay distinguishable once clusters from different batches/splits
#'   are compared side by side
#' @return list with `cluster_id_per_cell` (cluster id of each cell, in the order
#'   the cells appear in `adata`) and `mem` (cluster x marker MEM score matrix)
#' 
cluster_and_calculate_mem <- function(
  adata,
  layer_name,
  lineage_markers,
  fsom_param,
  label_prefix
) {
  # The layer comes out of AnnData without column names, so the markers have to
  # be named from var before they can be picked out.
  expr <- adata$layers[[layer_name]]
  colnames(expr) <- rownames(adata$var)
  expr <- expr[, lineage_markers]

  fsom <- FlowSOM(
    flowCore::flowFrame(expr),
    colsToUse = colnames(expr),
    nClus = fsom_param$n_clusters,
    xdim = fsom_param$xdim,
    ydim = fsom_param$ydim,
    seed = 42
  )
  cluster_id_per_cell <- paste0(label_prefix, GetMetaclusters(fsom))

  list(
    cluster_id_per_cell = cluster_id_per_cell,
    mem = compute_mem_matrix(
      expr_matrix = expr,
      cluster_id_per_cell = cluster_id_per_cell
    )
  )
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
  sample_prop_per_cluster <- merge(
    sample_prop_per_cluster, 
    sample_group_map, 
    by = "sample"
  )

  pvals <- sapply(clusters_to_test, function(cluster) {
    is_group_a <- sample_prop_per_cluster$group == group_a
    is_group_b <- sample_prop_per_cluster$group == group_b

    vals_a <- sample_prop_per_cluster$proportion[sample_prop_per_cluster$cluster == cluster & is_group_a]
    vals_b <- sample_prop_per_cluster$proportion[sample_prop_per_cluster$cluster == cluster & is_group_b]

    suppressWarnings(wilcox.test(vals_a, vals_b, alternative = "two.sided")$p.value)
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
#' @param mem_a,mem_b cluster x marker MEM matrices for two clusterings.
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

  is_cluster_a_best_match_mutual <- sapply(seq_len(length(clusters_a)), function(cluster_a_idx) {
    cluster_a <- clusters_a[cluster_a_idx]
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
  # note, nearest_cluster_in_b1 is a list. 
  # nearest slot is the one that store the name of the nearest cluster in b1.
  # similarity lot is the similarity matrix between integrated split clusters and b1 clusters.
  nearest_cluster_in_b1 <- find_nearest_by_mem(
    mem_query = clustering$mem, 
    mem_reference = b1_mem
  )

  # keep only those which b1 clusters are found to be DA and are mutually matched
  nearest_clusters_in_b1 <- nearest_cluster_in_b1$nearest
  # note, names are the integrated split clusters, and values are the nearest b1 clusters.
  nearest_clusters_in_b1 <- nearest_clusters_in_b1[nearest_clusters_in_b1 %in% retained_cluster_b1]
 
  # prepare the dataframe for DA
  cell_df_for_da <- data.frame(
    cluster_id_per_cell = clustering$cluster_id_per_cell,
    sample_id_per_cell = adata$obs$sample,
    group_id_per_cell = adata$obs$group,
    stringsAsFactors = FALSE
  )
  # add the mapped b1
  cell_df_for_da$mapped_b1_cluster <- nearest_clusters_in_b1[cell_df_for_da$cluster_id_per_cell]
  # Cells in a cluster that did not map to a retained b1 cluster are not tested,
  # but they still belong to their sample, so they need a label. Left as NA they
  # would be dropped by table(). Hence call them "unmatched".
  cell_df_for_da$mapped_b1_cluster[is.na(cell_df_for_da$mapped_b1_cluster)] <- "unmatched"

  # Cells are pooled by the b1 cluster they mapped to, so that is what gets
  # tested - one test per retained b1 cluster, not one per split cluster.
  abundance <- test_differential_abundance(
    cluster_id_per_cell = cell_df_for_da$mapped_b1_cluster,
    sample_id_per_cell = cell_df_for_da$sample_id_per_cell,
    group_id_per_cell = cell_df_for_da$group_id_per_cell,
    group_a = dataset_param$group_a,
    group_b = dataset_param$group_b,
    significance_threshold = significance_threshold,
    clusters_to_test = setdiff(unique(cell_df_for_da$mapped_b1_cluster), "unmatched")
  )

  list(
    da_results = abundance,
    cluster_id_per_cell = clustering$cluster_id_per_cell,
    mapped_b1_cluster_per_cell = cell_df_for_da$mapped_b1_cluster,
    similarity_to_batch1 = nearest_cluster_in_b1$similarity
  )
}
