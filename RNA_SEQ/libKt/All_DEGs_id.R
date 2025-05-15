# DEGs selection function 
All_DEGs_id <- function(folder_path, output_file = NULL) {
  # dir exist?
  if (!dir.exists(folder_path)) {
    stop("The specified folder does not exist.")
  }
  # obtain files path
  file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # csv?
  if (length(file_list) == 0) {
    stop("No CSV files found in the folder.")
  }
  
  columns_list <- lapply(file_list, function(file) {
    tryCatch({
      first_column <- read.csv(file, header = TRUE)[, 1] 
      return(as.character(first_column))
    }, error = function(e) {
      message(paste("Error reading file:", file, "-", e$message))
      return(NULL)
    })
  })
  
  columns_vector <- unlist(columns_list)
  unique_columns <- unique(columns_vector)
  result_df <- data.frame(column = unique_columns, stringsAsFactors = FALSE)
  
  # save result
  if (!is.null(output_file)) {
    write.csv(result_df, output_file, row.names = FALSE)
    message("Result saved to ", output_file)
  }
  
  # get result
  return(result_df)
}
