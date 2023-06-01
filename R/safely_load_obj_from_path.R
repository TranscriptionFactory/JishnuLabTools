# checking input
safely_load_obj_from_path = function(inputvar) {
  if ( !is.null(inputvar) ) {

    # check if its a yaml file path
    if ( typeof(inputvar) == "character" &&
         stringr::str_detect(stringr::str_to_lower(inputvar), pattern = ".yaml")) {

      # load the yaml
      yaml_input = yaml::yaml.load_file(inputvar)

      # load the x and y and return as dataframe with y as first column
      return(cbind.data.frame(y = check_for_df_or_path(yaml_input$y_path)),
             check_for_df_or_path(yaml_input$x_path))
    } else if (is.list(inputvar) && all(c("x_path", "y_path") %in% names(inputvar))) {
      # is it a loaded yaml file
      # same thing
      return(cbind.data.frame(y = check_for_df_or_path(yaml_input$y_path)),
             check_for_df_or_path(yaml_input$x_path))
    } else if (is.matrix(inputvar) | is.data.frame(inputvar) | is.list(inputvar)) {
      # check if input is an object or path
      # object passed
      return(inputvar)
    } else if (typeof(inputvar) == "character") {
      # path passed
      res = tryCatch({
        # try to read as an RDS
        return(as.data.frame(readRDS(inputvar)))
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
