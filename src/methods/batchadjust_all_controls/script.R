library(anndata)
library(flowCore)
## VIASH START
par <- list(
  input = "resources_test/.../input.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "batchadjust_all_controls"
)
## VIASH END

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))
source(paste0(meta$resources_dir, "/BatchAdjust.R"))

# As it only works with FCS files, the method requires substantial I/O
# the startegy used here is the following:
# 1. Read input AnnData file
# 1.1 Store the original order of cells in the Original_ID column
#     (to ensure that the order is restored after I/O operations)
# 2. Divide control cells and non-control cells
# 3. Write FCS file with control cells (from both groups)
#    (using BJ file naming convention)
# 4. Write FCS files with control cells
#    (using BJ file naming convention)
# 5. Run BatchAdjust on the FCS files
# 6. Read batch corrected FCS files via flowCore
# 7. Restore and check cell ordering via Original_ID column


cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])
#use Original_ID column to restore cell order after I/O operations
input$layers["preprocessed"][,"Original_ID"] <- seq(1,dim(input)[1])


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
