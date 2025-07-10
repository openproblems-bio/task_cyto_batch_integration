library(CellMixS)
library(anndataR) #for the conversion anndata -> sce
library(anndata)
library(SingleCellExperiment)
library(BiocParallel)
library(robustbase)

## VIASH START
par <- list(
  input_validation = "resources_test/.../validation.h5ad",
  input_unintegrated = "resources_test/.../unintegrated_censored.h5ad",
  input_integrated = "resources_test/.../integrated.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "cms"
)
## VIASH END

cat("Reading input files\n")
input_unintegrated <- anndataR::read_h5ad(par[["input_unintegrated"]])
input_integrated <- anndataR::read_h5ad(par[["input_integrated"]])

cat("Preprocessing SingleCellExperiment object\n")
#Fetch batch annotations form unintegrated layer
batch_key <- input_unintegrated$obs$batch
#Get markers to correct
markers_to_correct <- input_unintegrated$var$to_correct
#Convert to SingleCellExperiment
input_integrated_sce <- input_integrated$as_SingleCellExperiment()
input_integrated_sce <- input_integrated_sce[markers_to_correct, ]
colData(input_integrated_sce)["batch"] <- batch_key

cat("Compute Cell Mixing Score\n")
cms_res <- cms(
  input_integrated_sce,
  group = "batch",
  assay_name = "integrated",
  k = par[["n_neighbors"]],
  n_dim = par[["n_dim"]],
  BPPARAM = MulticoreParam(workers = 8)
)

cat("Compute Medcouple statistic\n")
cms_distr <- colData(cms_res)[, "cms"]
cms_mc <- robustbase::mc(cms_distr)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = input_integrated$uns$dataset_id,
    method_id = input_integrated$uns$method_id,
    metric_ids = "cms",
    metric_values = cms_mc,
    cms_parameters = list(
      n_neighbors = par[["n_neighbors"]],
      n_dim = par[["n_dim"]]
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
