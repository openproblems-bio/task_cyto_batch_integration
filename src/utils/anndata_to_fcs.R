library(anndata)
library(flowCore)
library(Biobase)
library(docstring)

anndata_to_fcs <- function(adata, out_dir) {
  #' @title Write FCS from anndata object
  #' 
  #' @description From an AnnData object creates a collection of .FCS files, one for each sample in the Anndata object.
  #' Files are written to the specified directory.
  #' NOTE: the AnnData object must have valid adata$obs$sample and adata$var$channel fields.
  #' @param h5ad_file An h5ad file with  adata$obs$sample and adata$var$channel fields.
  #' @param out_dir A directory where FCS files are written. If out_dir does not exists, it creates a new folder 
  #' @return A list of FCS files created in the specified directory.
  
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
  
  for (sample in sampl_names) {
    idx <- adata$obs$sample == sample
    events <- as.matrix(adata$layers[['preprocessed']][idx, ])
    colnames(events) <- channel_names
    ff <- flowFrame(exprs = events, parameters = par_df)
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    filename <- paste0(out_dir,"/",sample, ".fcs")
    write.FCS(ff, filename)
  }
  cat("FCS files created:\n")
  return(list.files(out_dir, pattern = '*.fcs'))
}