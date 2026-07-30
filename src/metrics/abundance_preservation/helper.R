library(FlowSOM)
library(flowCore)


#' Extract a layer's expression matrix from an AnnData object, subset to
#' specific markers and (optionally) specific cells.
#'
#' @param adata AnnData object
#' @param layer_name name of the layer to fetch (e.g. "preprocessed", "integrated")
#' @param markers character vector of var_names to keep
#' @param cell_mask optional logical vector to subset cells before selecting markers
#' @return numeric matrix, cells x markers
get_layer_expression <- function(adata, layer_name, markers, cell_mask = NULL) {
  expr <- adata$layers[[layer_name]]
  colnames(expr) <- rownames(adata$var)

  if (!is.null(cell_mask)) {
    expr <- expr[cell_mask, , drop = FALSE]
  }

  expr[, markers, drop = FALSE]
}


#' Compute MEM (Marker Enrichment Modeling) scores for a set of populations,
#' each compared against all remaining cells.
#'
#' Implements the MEM algorithm (Diggins et al. 2017, Nature Methods; see
#' https://github.com/cytolab/cytoMEM, MEM_function.R): for each population
#' and marker, MEM combines the difference in medians between the population
#' and the reference with a bonus/penalty from their relative IQRs, then
#' rescales every score in the matrix so the largest absolute value is 10.
#'
#' @param expr_matrix numeric matrix, cells x markers
#' @param population_labels character vector, per-cell population label,
#'   length nrow(expr_matrix). May include labels that are not in
#'   `populations_of_interest` (e.g. unlabelled cells) - those cells still
#'   contribute to the reference ("all other cells") pool but are never
#'   scored themselves.
#' @param populations_of_interest character vector, the labels to compute a
#'   MEM score row for.
#' @param iqr_thresh numeric floor applied to per-population IQRs to avoid
#'   divide-by-zero (default 0.5, matching the reference implementation).
#' @return numeric matrix, populations_of_interest x colnames(expr_matrix)
compute_mem_matrix <- function(expr_matrix, population_labels, populations_of_interest, iqr_thresh = 0.5) {
  population_labels <- as.character(population_labels)
  populations_of_interest <- as.character(populations_of_interest)
  markers <- colnames(expr_matrix)

  mag_pop <- matrix(
    nrow = length(populations_of_interest), ncol = length(markers),
    dimnames = list(populations_of_interest, markers)
  )
  mag_ref <- mag_pop
  iqr_pop <- mag_pop
  iqr_ref <- mag_pop

  for (pop in populations_of_interest) {
    in_pop <- population_labels == pop

    mag_pop[pop, ] <- abs(apply(expr_matrix[in_pop, , drop = FALSE], 2, median, na.rm = TRUE))
    iqr_pop[pop, ] <- apply(expr_matrix[in_pop, , drop = FALSE], 2, IQR, na.rm = TRUE)

    mag_ref[pop, ] <- abs(apply(expr_matrix[!in_pop, , drop = FALSE], 2, median, na.rm = TRUE))
    iqr_ref[pop, ] <- apply(expr_matrix[!in_pop, , drop = FALSE], 2, IQR, na.rm = TRUE)
  }

  # Floor the IQRs so a population with (near-)zero spread doesn't blow up the ratio below.
  iqr_pop[iqr_pop < iqr_thresh] <- iqr_thresh
  iqr_ref[iqr_ref < iqr_thresh] <- iqr_thresh

  # When a population's spread on this marker is tighter than the reference's
  # (IQRpop < IQRref), this ratio is positive and gets added to the score
  # below, pushing it further from zero - tight agreement within a
  # population counts as stronger evidence of enrichment/depletion. When the
  # population is more spread out than the reference, the ratio would be
  # negative; it's floored at zero instead, so it's never subtracted.
  iqr_ratio <- (iqr_ref / iqr_pop) - 1
  iqr_ratio[iqr_ratio < 0] <- 0

  mag_diff <- mag_pop - mag_ref
  mem_matrix <- abs(mag_diff) + iqr_ratio

  # Markers that are lower in the population than in the reference get a negative score.
  mem_matrix[mag_diff < 0] <- -mem_matrix[mag_diff < 0]

  scale_max <- max(abs(mem_matrix))
  if (scale_max > 0) {
    mem_matrix <- (mem_matrix / scale_max) * 10
  }

  mem_matrix
}


#' MEM RMSD similarity between every pair of MEM score vectors from two
#' populations x markers matrices, following the "pop files" comparison mode
#' in https://github.com/cytolab/cytoMEM (MEM_RMSD.R): root-mean-square
#' deviation over the markers common to both, rescaled to a 0-100
#' "similarity" score. MEM scores are bounded to [-10, 10], so the largest
#' possible per-marker gap between two populations is 20; RMSD is rescaled
#' against that so 100 means identical MEM vectors and 0 means opposite
#' extremes on every marker.
#'
#' @param mem_a matrix, populations x markers
#' @param mem_b matrix, populations x markers
#' @return matrix, nrow(mem_a) x nrow(mem_b), similarity over the markers common to both
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


#' Map each metacluster (a row of `mem_metacluster`) to the population (a row
#' of `mem_population`) with the highest MEM RMSD similarity. Several
#' metaclusters can map to the same population.
#'
#' @param mem_metacluster matrix, FlowSOM metaclusters x markers
#' @param mem_population matrix, reference cell types x markers
#' @return list with:
#'   `mapping`: named character vector, names are rownames(mem_metacluster),
#'     values are the matched rownames(mem_population)
#'   `similarity`: mem_metacluster x mem_population similarity matrix (see
#'     compute_mem_rmsd()), kept for diagnostics/plotting
map_populations_by_mem <- function(mem_metacluster, mem_population) {
  similarity <- compute_mem_rmsd(mem_metacluster, mem_population)
  matched <- apply(similarity, 1, function(row) colnames(similarity)[which.max(row)])

  list(
    mapping = setNames(matched, rownames(similarity)),
    similarity = similarity
  )
}


#' For each sample, compute the proportion of its cells belonging to each
#' cell type (either the ground-truth cell type, or a metacluster-derived
#' cell type label).
#'
#' @param cell_type character vector, per-cell population label
#' @param sample character vector, per-cell sample id
#' @param group character vector, per-cell biological group label
#' @return data frame with columns sample, group, cell_type, proportion.
#'   Every sample x cell_type combination is included (0 if absent), so every
#'   sample contributes a value for every cell type.
compute_sample_proportions <- function(cell_type, sample, group) {
  cell_type <- as.character(cell_type)
  sample <- as.character(sample)
  group <- as.character(group)

  count_table <- table(sample, cell_type)
  sample_totals <- rowSums(count_table)
  proportion_table <- count_table / sample_totals

  proportions <- as.data.frame(as.table(proportion_table), stringsAsFactors = FALSE)
  colnames(proportions) <- c("sample", "cell_type", "proportion")

  sample_group_map <- unique(data.frame(sample = sample, group = group, stringsAsFactors = FALSE))
  merge(proportions, sample_group_map, by = "sample")
}


#' Wilcoxon rank-sum test comparing per-sample proportions between two
#' groups, run independently for every cell type present in `proportions`,
#' except unlabelled cells: those only ever contribute to the sample totals
#' upstream (in compute_sample_proportions) and are never tested themselves.
#'
#' @return named numeric vector, cell_type -> p-value
run_wilcoxon_per_celltype <- function(proportions, group_a, group_b) {
  cell_types <- unique(proportions$cell_type)
  cell_types <- cell_types[!tolower(cell_types) %in% c("unlabelled", "unlabeled")]

  pvals <- sapply(cell_types, function(ct) {
    vals_a <- proportions$proportion[proportions$cell_type == ct & proportions$group == group_a]
    vals_b <- proportions$proportion[proportions$cell_type == ct & proportions$group == group_b]

    suppressWarnings(wilcox.test(vals_a, vals_b, alternative = "two.sided")$p.value)
  })
  names(pvals) <- cell_types

  pvals
}


#' Cluster one integrated split into FlowSOM metaclusters, map each
#' metacluster to a reference cell type by MEM RMSD, and test whether the
#' resulting (metacluster-derived) cell type proportions differ between
#' groups.
#'
#' @param adata_full AnnData object, a single integrated split, nocontrols,
#'   unlabelled cells included
#' @param lineage_markers character vector of var_names to cluster on
#' @param n_clusters target number of metaclusters
#' @param xdim,ydim SOM grid dimensions
#' @param mem_celltype_reference matrix, reference cell types x lineage
#'   markers, as produced by compute_mem_matrix() on the unintegrated data
#' @param group_a,group_b the two obs$group values being compared
#' @param significance_threshold p-value cutoff used to call a cell type
#'   difference significant
#' @return list with `mapping` (metacluster -> cell type), `similarity`
#'   (metacluster x reference cell type MEM similarity matrix, for
#'   diagnostics/plotting), `pvals` (named by cell type) and `sig` (cell
#'   types with pvals <= significance_threshold)
process_integrated_split <- function(
  adata_full,
  lineage_markers,
  n_clusters,
  xdim,
  ydim,
  mem_celltype_reference,
  group_a,
  group_b,
  significance_threshold
) {
  # FlowSOM clusters every cell in the split, unlabelled ones included: they
  # get a metacluster (and therefore a mapped cell type) just like any other
  # cell.
  expression_matrix <- adata_full$layers[["integrated"]]
  colnames(expression_matrix) <- rownames(adata_full$var)
  flowframe <- flowCore::flowFrame(expression_matrix)

  fsom <- FlowSOM(
    flowframe,
    colsToUse = lineage_markers,
    nClus = n_clusters,
    xdim = xdim,
    ydim = ydim,
    seed = 42
  )
  metaclusters <- paste0("meta_", GetMetaclusters(fsom))

  mem_metacluster <- compute_mem_matrix(
    get_layer_expression(adata_full, "integrated", lineage_markers),
    metaclusters,
    unique(metaclusters)
  )

  # Several metaclusters can map to the same cell type; their cells are
  # simply pooled together under that cell type for the proportion/wilcoxon
  # test below.
  mapping_result <- map_populations_by_mem(mem_metacluster, mem_celltype_reference)
  metacluster_to_celltype <- mapping_result$mapping
  mapped_cell_type <- metacluster_to_celltype[metaclusters]

  proportions <- compute_sample_proportions(
    cell_type = mapped_cell_type,
    sample = adata_full$obs$sample,
    group = adata_full$obs$group
  )
  pvals <- run_wilcoxon_per_celltype(proportions, group_a, group_b)

  list(
    mapping = metacluster_to_celltype,
    similarity = mapping_result$similarity,
    pvals = pvals,
    sig = names(pvals)[pvals <= significance_threshold]
  )
}
