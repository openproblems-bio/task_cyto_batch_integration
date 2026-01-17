library(anndata)
library(flowCore)
## VIASH START
par <- list(
  input = "resources_test/debug/batchadjust/_viash_par/input_1/censored_split1.h5ad",
  output = "resources_test/debug/batchadjust/output.h5ad",
  percentile = as.integer('80')
)
meta <- list(
  name = "batchadjust_all_controls",
  temp_dir = "resources_test/tmp",
  resources_dir = "src/methods/batchadjust_one_control"
)
## VIASH END

source(paste0(meta$resources_dir, "/utils.R"))
source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))
source(paste0(meta$resources_dir, "/BatchAdjust.R"))

# only for HPC, the idea is if running on HPC, use a temp dir set in the env variable
tmp_dir <- Sys.getenv("HPC_VIASH_META_TEMP_DIR")
if (tmp_dir != "") {
  # Environment variable is set, use it
  print(paste0("Using HPC temp dir from env: ", tmp_dir))
} else {
  # Environment variable not set, use meta
  tmp_dir <- meta[["temp_dir"]]
  print(paste0("Using meta temp dir: ", tmp_dir))
}

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
# input$layers["preprocessed"][, "Original_ID"] <- seq(1, dim(input)[1])

original_id_in_var <- "Original_ID" %in% input$var_names

input <- add_original_id(input)

cat("Split cells\n")
input_controls <- input[input$obs$is_control == 1, ]
input_no_controls <- input[input$obs$is_control != 1, ]

print("Cells in control sample:")
print(unique(input_controls$obs$sample))
print(input_controls)
print("Cells in non-control sample:")
print(unique(input_no_controls$obs$sample))
print(input_no_controls)

# avoid NA due to invalid factor level
input_controls$obs$sample <- as.character(input_controls$obs$sample)

# Set sample names for batch-specific control files
input_controls$obs$sample[input_controls$obs$batch == 1] <- "Batch1_anchor"
input_controls$obs$sample[input_controls$obs$batch == 2] <- "Batch2_anchor"

# process non control samples' sample names
input_no_controls$obs$sample <- as.character(input_no_controls$obs$sample)

# the non control sampels also need to have batch in the sample name.. Doh..
non_control_has_batch <- any(grepl("Batch", input_no_controls$obs$sample))

if (!non_control_has_batch) {
  cat(
    "Non control samples do not have batch info in sample names!\n",
    "Sample names before modification: ", 
    paste0(unique(input_no_controls$obs$sample), collapse = ", "), "\n",
    "Modifying sample names to include batch info!\n"
  )

  input_no_controls$obs$sample <- paste0(input_no_controls$obs$sample, "_Batch", input_no_controls$obs$batch, "_")

  cat(
    "Sample names after modification:\n", 
    paste0(unique(input_no_controls$obs$sample), collapse = ", "), 
    "\n"
  )
}

# make sure there is _ after the batch1 or batch2, otherwise batchadjust won't find the fcs files.
input_no_controls$obs$sample <- sapply(input_no_controls$obs$sample, fix_batch_underscore_anynum)


cat("Writing FCS files\n")
anndata_to_fcs(input_controls, out_dir = tmp_dir)
anndata_to_fcs(input_no_controls, out_dir =  tmp_dir)

cat("Writing channels to correct as a text file\n")
markers_to_correct <- as.character(
  input$var$channel[input$var$to_correct == TRUE]
)
writeLines(markers_to_correct,
           con = paste0(tmp_dir, "/to_correct_list.txt"))

cat("Running BatchAdjust\n")
perc <- paste0(as.character(par[["percentile"]]), "p")
output_dir <- paste0(tmp_dir, "/corrected")

BatchAdjust(
  basedir = tmp_dir,
  outdir = output_dir,
  channelsFile = paste0(tmp_dir, "/to_correct_list.txt"),
  anchorKeyword = "anchor",
  batchKeyword = "Batch", #skip 'b' to make it robust to upper/lowercase
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

# Remove Original_ID if it was not there in the beginning
if (!original_id_in_var) {
    corrected_matrix$Original_ID <- NULL
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
