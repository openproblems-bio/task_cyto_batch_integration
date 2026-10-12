#!/usr/bin/env Rscript
#
# Run an R metric ad-hoc for all (dataset, method) combinations found in a
# pipeline output directory, without re-running the full pipeline.
#
# The script discovers method outputs under:
#     <input-dir>/<dataset>/method_out/<method>_split1.h5ad
#     <input-dir>/<dataset>/method_out/<method>_split2.h5ad
#
# and writes metric outputs to:
#     <output-dir>/<dataset>/<metric>_<method>.h5ad
#
# It works by replacing the ## VIASH START...END block in the metric script
# with the correct par/meta values and running the result in a subprocess.
# This means any R metric in src/metrics/ can be run ad-hoc with this script.
# (This is the R counterpart to run_metric_adhoc.py, for Python metrics.)
#
# Usage:
#   Rscript scripts/adhoc_runs/run_metric_adhoc.R \
#       --metric abundance_preservation \
#       --input-dir /path/to/run_2026-04-26_21-34-33 \
#       --datasets-dir /path/to/datasets \
#       --output-dir /path/to/analysis_2026-05-13/abundance_preservation
#
#   # Restrict to specific datasets or methods:
#   Rscript scripts/adhoc_runs/run_metric_adhoc.R \
#       --metric abundance_preservation \
#       --input-dir /path/to/run_2026-04-26_21-34-33 \
#       --datasets-dir /path/to/datasets \
#       --output-dir /path/to/output \
#       --datasets mouse_spleen_flow_cytometry \
#       --methods harmonypy combat
#
#   # Skip already-computed outputs:
#   Rscript scripts/adhoc_runs/run_metric_adhoc.R ... --skip-existing

#' Directory containing the currently running script, resolved via Rscript's
#' own --file= argument (there is no R equivalent of Python's __file__).
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Could not determine this script's location (not run via Rscript?).")
  }
  dirname(normalizePath(sub("^--file=", "", file_arg)))
}

#' Walk upward from start_dir until a directory containing _viash.yaml is found.
find_repo_root <- function(start_dir) {
  directory <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(directory, "_viash.yaml"))) {
      return(directory)
    }
    parent <- dirname(directory)
    if (parent == directory) {
      stop(sprintf("Could not locate repo root: no _viash.yaml found above %s.", start_dir))
    }
    directory <- parent
  }
}

REPO_ROOT <- find_repo_root(get_script_dir())


#' Locate the split1 and split2 output files for a (dataset, method) pair.
#'
#' @return list(split1 = ..., split2 = ...), or NULL if either file is missing
find_method_splits <- function(input_dir, dataset, method) {
  method_out_dir <- file.path(input_dir, dataset, "method_out")
  split1 <- file.path(method_out_dir, paste0(method, "_split1.h5ad"))
  split2 <- file.path(method_out_dir, paste0(method, "_split2.h5ad"))

  if (!file.exists(split1)) {
    cat(sprintf("    WARNING: split1 not found: %s\n", split1))
    return(NULL)
  }
  if (!file.exists(split2)) {
    cat(sprintf("    WARNING: split2 not found: %s\n", split2))
    return(NULL)
  }

  list(split1 = split1, split2 = split2)
}


#' Return all method names found as <method>_split1.h5ad in method_out_dir.
discover_methods <- function(method_out_dir) {
  files <- list.files(method_out_dir, pattern = "_split1\\.h5ad$")
  sort(sub("_split1\\.h5ad$", "", files))
}


#' Render a named list of R values as literal `list(...)` source text, so it
#' can be spliced into the metric script as `par <- ...` / `meta <- ...`.
r_list_literal <- function(named_list) {
  entries <- vapply(names(named_list), function(key) {
    sprintf("%s = %s", key, deparse(named_list[[key]]))
  }, character(1))
  paste0("list(\n  ", paste(entries, collapse = ",\n  "), "\n)")
}


#' Inject par/meta into the metric script and execute it in a subprocess.
#'
#' The ## VIASH START...END block in script.R is replaced with the provided
#' par and meta values. meta$resources_dir is set to the metric's own source
#' directory (src/metrics/<metric>/): this is where helper.R already lives,
#' and it's where a real viash build would also place the shared
#' helper_functions.R (viash copies every declared resource into one flat
#' directory). Since this script runs script.R directly rather than through a
#' real build, the shared helper is temporarily copied in next to helper.R
#' for the duration of the run, then removed again.
run_metric <- function(metric_name, par) {
  script_path <- file.path(REPO_ROOT, "src", "metrics", metric_name, "script.R")
  if (!file.exists(script_path)) {
    stop(sprintf(
      "Metric script not found: %s\nCheck that '%s' is a valid metric under src/metrics/.",
      script_path, metric_name
    ))
  }
  metric_dir <- dirname(script_path)

  meta <- list(name = metric_name, resources_dir = metric_dir)
  injected_block <- paste0(
    "## VIASH START\n",
    "par <- ", r_list_literal(par), "\n",
    "meta <- ", r_list_literal(meta), "\n",
    "## VIASH END"
  )

  original_text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  match_pos <- regexpr("(?s)## VIASH START.*?## VIASH END", original_text, perl = TRUE)
  if (match_pos == -1) {
    stop(sprintf("Could not find a '## VIASH START ... ## VIASH END' block in %s", script_path))
  }
  match_len <- attr(match_pos, "match.length")
  patched <- paste0(
    substr(original_text, 1, match_pos - 1),
    injected_block,
    substr(original_text, match_pos + match_len, nchar(original_text))
  )

  shared_helper_src <- file.path(REPO_ROOT, "src", "utils", "helper_functions.R")
  shared_helper_dst <- file.path(metric_dir, "helper_functions.R")
  shared_helper_existed <- file.exists(shared_helper_dst)
  if (!shared_helper_existed && file.exists(shared_helper_src)) {
    file.copy(shared_helper_src, shared_helper_dst)
  }

  temp_script <- tempfile(pattern = "adhoc_", tmpdir = metric_dir, fileext = ".R")
  writeLines(patched, temp_script)

  on.exit({
    unlink(temp_script)
    if (!shared_helper_existed) unlink(shared_helper_dst)
  }, add = TRUE)

  exit_code <- system2("Rscript", shQuote(temp_script))
  if (exit_code != 0) {
    stop(sprintf("Metric run failed (exit code %s): %s", exit_code, metric_name))
  }
}


#' Parse this script's own CLI arguments (base R only, no CLI package needed).
parse_args <- function(raw_args) {
  args <- list(
    metric = NULL, input_dir = NULL, datasets_dir = NULL, output_dir = NULL,
    datasets = NULL, methods = NULL, skip_existing = FALSE
  )

  # --datasets/--methods take one or more values, ending at the next flag or
  # the end of the argument list.
  take_values <- function(start_i) {
    values <- character(0)
    i <- start_i
    while (i <= length(raw_args) && !startsWith(raw_args[i], "--")) {
      values <- c(values, raw_args[i])
      i <- i + 1
    }
    list(values = values, next_i = i)
  }

  i <- 1
  while (i <= length(raw_args)) {
    flag <- raw_args[i]
    if (flag == "--metric") {
      args$metric <- raw_args[i + 1]; i <- i + 2
    } else if (flag == "--input-dir") {
      args$input_dir <- raw_args[i + 1]; i <- i + 2
    } else if (flag == "--datasets-dir") {
      args$datasets_dir <- raw_args[i + 1]; i <- i + 2
    } else if (flag == "--output-dir") {
      args$output_dir <- raw_args[i + 1]; i <- i + 2
    } else if (flag == "--datasets") {
      taken <- take_values(i + 1); args$datasets <- taken$values; i <- taken$next_i
    } else if (flag == "--methods") {
      taken <- take_values(i + 1); args$methods <- taken$values; i <- taken$next_i
    } else if (flag == "--skip-existing") {
      args$skip_existing <- TRUE; i <- i + 1
    } else {
      stop(sprintf("Unknown argument: %s", flag))
    }
  }

  required <- c("metric", "input_dir", "datasets_dir", "output_dir")
  missing <- required[vapply(args[required], is.null, logical(1))]
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required argument(s): %s",
      paste0("--", gsub("_", "-", missing), collapse = ", ")
    ))
  }

  args
}


main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  input_dir <- normalizePath(args$input_dir, mustWork = TRUE)
  datasets_dir <- normalizePath(args$datasets_dir, mustWork = TRUE)
  output_dir <- normalizePath(args$output_dir, mustWork = FALSE)

  # --datasets restricts the run to an explicit list; otherwise fall back to
  # every subdirectory of input_dir, since each dataset is expected to have
  # its own top-level folder there (<input-dir>/<dataset>/method_out/...).
  datasets <- args$datasets
  if (is.null(datasets)) {
    datasets <- sort(basename(list.dirs(input_dir, recursive = FALSE)))
  }

  for (dataset in datasets) {
    cat(sprintf("\nDataset: %s\n", dataset))

    method_out_dir <- file.path(input_dir, dataset, "method_out")
    if (!dir.exists(method_out_dir)) {
      cat("  No method_out directory found, skipping.\n")
      next
    }

    unintegrated <- file.path(datasets_dir, dataset, "unintegrated.h5ad")
    if (!file.exists(unintegrated)) {
      cat(sprintf("  Unintegrated file not found: %s, skipping.\n", unintegrated))
      next
    }

    # --methods restricts to an explicit list; otherwise discover every
    # method that produced output for this dataset.
    methods <- args$methods
    if (is.null(methods)) {
      methods <- discover_methods(method_out_dir)
    }
    if (length(methods) == 0) {
      cat(sprintf("  No method outputs found in %s, skipping.\n", method_out_dir))
      next
    }

    for (method in methods) {
      cat(sprintf("  Method: %s\n", method))

      # Check --skip-existing before doing any of the comparatively expensive
      # file lookups below, so re-runs over a mostly-complete output
      # directory stay cheap.
      output_path <- file.path(output_dir, dataset, sprintf("%s_%s.h5ad", args$metric, method))
      if (args$skip_existing && file.exists(output_path)) {
        cat(sprintf("    Output exists, skipping: %s\n", output_path))
        next
      }

      splits <- find_method_splits(input_dir, dataset, method)
      if (is.null(splits)) next

      dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

      # These four keys are exactly what every metric's VIASH par block
      # expects (see any src/metrics/*/script.R).
      par <- list(
        input_unintegrated = unintegrated,
        input_integrated_split1 = splits$split1,
        input_integrated_split2 = splits$split2,
        output = output_path
      )

      cat(sprintf("    Output: %s\n", output_path))
      run_metric(args$metric, par)
      cat("    Done.\n")
    }
  }

  cat("\nAll done.\n")
}

main()
