library(anndata)
library(flowCore)
library(Biobase)

# Take a h5ad file as an input and return a collection of fcs for each sample in the Anndata object
# NOTE: the AnnData object must have valid adata$obs$sample and adata$var$channel fields.
# Example usage : `write_fcs_from_h5ad('mydata.h5ad',"path_to_dir")`

write_fcs_from_h5ad <- function(h5ad_file, out_dir) {
  adata <- anndata::read_h5ad(h5ad_file)
  
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
    events <- as.matrix(adata$layers$preprocessed[idx, ])
    colnames(events) <- channel_names
    ff <- flowFrame(exprs = events, parameters = par_df)
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    filename <- paste0(out_dir,"/",sample, ".fcs")
    write.FCS(ff, filename)
  }
}
