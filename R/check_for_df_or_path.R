# checking input
check_for_df_or_path = function(inputvar) {
  if ( !is.null(inputvar) ) {
    # check if input is an object or path
    if (is.matrix(inputvar) | is.data.frame(inputvar) | is.list(inputvar)) {
      # object passed
      return(inputvar)
    } else if (typeof(inputvar) == "character") {
      # path passed
      res = tryCatch({
        # try to read as an RDS
        as.data.frame(readRDS(inputvar))
      }, error = function(e) {
        # assume its a matrix
        mat = tryCatch({
          # try to read as an RDS
          as.matrix(read.csv(inputvar, row.names = 1))
        }, error = function(e) {
          cat("Check path - path should be to csv or RDS file \n")
        }, finally = {})
        return(mat)
      }, finally = {})
      return(res)
    }
  } else {
    # didn't work
    cat("input not recognized as RDS, matrix, dataframe or list. \n")
    return(NULL)
  }
}
