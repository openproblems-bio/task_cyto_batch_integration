library(anndata)

## VIASH START
par <- list(
  input = "resources_test/.../input.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "rpca_to_goal"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])

cat("Preprocess data\n")
# ... preprocessing ...

cat("Train model\n")
# ... train model ...

cat("Generate predictions\n")
# ... generate predictions ...

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  
)
output$write_h5ad(par[["output"]], compression = "gzip")
