library(anndata)
fix_batch_underscore_anynum <- function(x) {
    # Pattern: Batch followed by one or more digits, NOT followed by underscore
    pattern <- "(Batch\\d+)(?!_)"
    
    # Add underscore if missing
    x <- sub(pattern, "\\1_", x, perl = TRUE)
    
    return(x)
}

add_original_id <- function(input) {
    if (!"Original_ID" %in% input$var_names) {
        cat("Adding Original_ID to var and recreating input anndata\n")
        old_var <- input$var
        old_var[] <- lapply(old_var, function(x) if (is.factor(x)) as.character(x) else x)
        old_var <- rbind(old_var,
            data.frame(
                numeric_id=length(input$var_names) + 1,
                channel="Original_ID",
                marker="",
                marker_type="other",
                to_correct=FALSE,
                row.names = "Original_ID"
            )
        )
        old_var[] <- lapply(old_var, function(x) if (is.character(x)) as.factor(x) else x)
        new_mat <- Matrix::as.matrix(input$layers[["preprocessed"]])
        new_mat <- cbind(new_mat, Original_ID = seq_len(nrow(new_mat)))

        # create new anndata
        input <- anndata::AnnData(
            X = new_mat,
            obs = input$obs,
            var = old_var,
            layers = list(preprocessed = new_mat),
            uns = list(
                dataset_description = input$uns$dataset_description,
                dataset_id = input$uns$dataset_id,
                dataset_name = input$uns$dataset_name,
                dataset_organism = input$uns$dataset_organism,
                dataset_reference = input$uns$dataset_reference,
                dataset_summary = input$uns$dataset_summary,
                dataset_url = input$uns$dataset_url
            )
        )
    } else {
        cat("Adding new Original_ID var\n")
        input$layers["preprocessed"][, "Original_ID"] <- seq(1, dim(input)[1])
    }
    return(input)
}
