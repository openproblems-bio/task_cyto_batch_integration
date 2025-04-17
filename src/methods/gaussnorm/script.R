library(flowCore)
library(anndata)
library(flowStats)

## VIASH START
par <- list(
  input = "resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "gaussNorm",
  temp_dir: '/tmp'
)
## VIASH ENDs

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))
tmp_path <- meta[["temp_dir"]]

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input"]])
markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

cat("Creating FCS files from anndata\n")
anndata_to_fcs(adata, tmp_path)

cat("Creating FlowSet from FCS files\n")
fcs_files <- list.files(path = tmp_path, pattern = '*.fcs', full.names = TRUE)
fset <- read.flowSet(files = fcs_files)
print(fset)

cat("Run gaussNorm\n")
for(marker in markers_to_correct){
  print(marker)
  pass<- tryCatch({
    fset <- gaussNorm(fset, channel.names = marker)$flowset},
    error = function(e) {FALSE})
  
  if(isFALSE(pass)){
    print(paste0("Reducing n. of landmarks for ",marker))
    fset <- gaussNorm(fset, channel.names = marker, max.lms = 1)$flowset
  }
}
# Concatenate all flowFrames in a FlowSet
integrated_matrix <- fsApply(fset,exprs)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = adata$obs[,integer(0)],
  var = adata$var[as.vector(colnames(integrated_matrix)), integer(0)],
  layers = list(integrated = integrated_matrix),
  uns = list(
    dataset_id = adata$uns$dataset_id,
    method_id = meta$name,
    parameters = list()
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
