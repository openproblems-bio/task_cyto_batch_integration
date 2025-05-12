library(anndata)
library(flowCore)
library(Biobase)
library(docstring)

#' @title Write FCS from anndata object. Alternatively, read a FlowSet from anndata object
#'
#' @description When `out_dir` is specified, it creates a collection of FCS flies from an anndata object and returns a list with he path of the created FCS files.
#' If no directory is specified, the function reads a FlowSet object from an anndata object.
#' NOTE: the AnnData object must have valid adata$obs$sample and adata$var$channel fields.
#' @param adata Anndata object with valid adata$obs$sample and adata$var$channel fields.
#' @param out_dir A directory where FCS files are written. If out_dir does not exists, it creates a new folder. If NULL, it does not create FCS files, but only reads from anndata to FlowSet. Default NULL.
#' @param layer_name String corresponding to a layer in adata$layers. Default 'preprocessed'.
#' @return A FlowSet 
anndata_to_fcs <- function(adata, out_dir= NULL, layer_name = 'preprocessed') {
 
  channel_names <- adata$var$channel
  sampl_names <- unique(adata$obs$sample)

  par_df <- AnnotatedDataFrame(
    data.frame(
      name = channel_names,
      desc = NA,
      range = NA,
      minRange = NA,
      maxRange = NA
    )
  )

  if (is.null(out_dir)){
    print("WARNING: out_dir not specified, creating FlowSet from Anndata")
    ff_list <- list()
    for (sample in sampl_names) {
      idx <- adata$obs$sample == sample
      events <- as.matrix(adata$layers[[layer_name]][idx, ])
      colnames(events) <- channel_names
      ff <- flowFrame(exprs = events, parameters = par_df)
      ff_list[[sample]] <- ff
    }
    fset <- as(ff_list, "flowSet")
    return(fset)

  } else {

    for (sample in sampl_names) {
      idx <- adata$obs$sample == sample
      events <- as.matrix(adata$layers[[layer_name]][idx, ])
      colnames(events) <- channel_names
      ff <- flowFrame(exprs = events, parameters = par_df)

      if (!dir.exists(out_dir)) {
        dir.create(out_dir)
      }
      filename <- paste0(out_dir,"/",sample, ".fcs")
      write.FCS(ff, filename)
    }
    fcs_files <- list.files(path = out_dir, pattern = "*.fcs", full.names = TRUE)
    return(fcs_files)
  }
}
