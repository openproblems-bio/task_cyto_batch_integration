library(anndata)
library(flowCore)
## VIASH START
par <- list(
  input = "resources_test/.../input.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "batchadjust_all_controls",
  temp_dir = "/tmp"
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
# 8. Convert the corrected matrix to AnnData


cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])
#use Original_ID column to restore cell order after I/O operations
input$layers["preprocessed"][, "Original_ID"] <- seq(1, dim(input)[1])

cat("Split cells\n")
input_controls <- input[input$obs$is_control != 0, ]
input_no_controls <- input[input$obs$is_control == 0, ]

print("Cells in control sample:")
print(unique(input_controls$obs$sample))
print(input_controls)
print("Cells in non-control sample:")
print(unique(input_no_controls$obs$sample))
print(input_no_controls)

#avoid NA due to invalid factor level
input_controls$obs$sample <- as.character(input_controls$obs$sample)
# Set sample names for batch-specific control files
input_controls$obs$sample[input_controls$obs$batch == 1] <- "0Batch1_anchor"
input_controls$obs$sample[input_controls$obs$batch == 2] <- "0Batch2_anchor"

cat("Writing FCS files\n")
anndata_to_fcs(input_controls, out_dir = meta[["temp_dir"]])
anndata_to_fcs(input_no_controls, out_dir =  meta[["temp_dir"]])

cat("Writing channels to correct as a text file\n")
markers_to_correct <- as.character(
  input$var$channel[input$var$to_correct == TRUE]
)
writeLines(markers_to_correct,
           con = paste0(meta[["temp_dir"]], "/to_correct_list.txt"))

cat("Running BatchAdjust\n")
perc <- paste0(as.character(par[["percentile"]]), "p")
output_dir <- paste0(meta[["temp_dir"]], "/corrected")

BatchAdjust(
  basedir = meta[["temp_dir"]],
  outdir = output_dir,
  channelsFile = paste0(meta[["temp_dir"]], "/to_correct_list.txt"),
  anchorKeyword = "anchor",
  batchKeyword = "atch", #skip 'b' to make it robust to upper/lowercase
  method = perc,
  transformation = FALSE,
  addExt = NULL,
  plotDiagnostics = FALSE
)

print(list.files(output_dir))
cat("Reading FCS files\n")
fs <- read.flowSet(path = output_dir, pattern = "*.fcs")
fs_all <- as(fs, "flowFrame")
corrected_matrix <- as.data.frame(exprs(fs_all))
#Remove column created when converting to flowframe
corrected_matrix$Original <- NULL
colnames(corrected_matrix) <- input$var_names

cat("Converting to Anndata + ordering check\n")
#Sort cells and check order
corrected_matrix <- corrected_matrix[order(corrected_matrix$Original_ID), ]
order_check <- corrected_matrix$Original_ID == input$layers[['preprocessed']][, 'Original_ID']
if (FALSE %in% order_check) {
  stop("Failed in restoring indexing")
}

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = input$obs[, integer(0)],
  var = input$var[colnames(corrected_matrix), integer(0)],
  layers = list(integrated = corrected_matrix),
  uns = list(
    dataset_id = input$uns$dataset_id,
    method_id = meta$name,
    parameters = list(
      percentile = par[["percentile"]]
    )
  )
)

print(output)
output$write_h5ad(par[["output"]], compression = "gzip")
