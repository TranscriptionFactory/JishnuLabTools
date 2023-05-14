#' @export
# load data as a single matrix with first column as Y values
# obj is either a loaded dataframe (matrix, list, or dataframe) or a path to one
# out_path is from yaml file
load_data = function(obj = NULL, yaml, remove_mito_ribo = F, create_dir = T, quantile_filter = 0,
remove_zero_median_cols = F) {

  # did you pass a yaml file
  if ( typeof(yaml) == "character" &&
       stringr::str_detect(stringr::str_to_lower(inputvar), pattern = ".yaml")) {

    # load the yaml
    loaded_yaml = yaml::yaml.load_file(yaml)
  } else {
    # passed a yaml file thats already been loaded (e.g. was an RDS)
    loaded_yaml = yaml
  }

  if (! dir.exists(loaded_yaml$out_path)) {
    if (create_dir) {
      dir.create(loaded_yaml$out_path)
    } else {
      cat("Directory doesn't exist. Check yaml out_path or set create_dir = T \n")
      return()
    }
  }

  if ( !is.null(obj)) {
    # use this parameter if you don't want to deal with editing yaml files
    df = JishnuLabTools:::safely_load_obj_from_path(obj)
  } else {
    # load data from yaml file
    df = JishnuLabTools::safely_load_obj_from_path(yaml)
  }

  if (remove_mito_ribo) {
    df = cbind.data.frame(df[, 1], JishnuLabTools::remove_mitochondrial_ribosomal_genes(df[, -1]))
  }

  if (! is.null(df)) {
    # shuffle rows
    df = df[sample(1:nrow(df)),]

    x = as.matrix(df[, -1])
    y = as.matrix(df[, 1])

    # will just remove columns with zero standard deviation.
    # to remove based on quantile, and set quantile

    cleaned = JishnuLabTools::clean_data(x = x, y = y, edit_data = T,
                                         quantile_filter = quantile_filter, remove_zero_median_cols = remove_zero_median_cols)




    # we are going to save these edited data to our output folder
    xpath = paste0(loaded_yaml$out_path, "x.csv")
    ypath = paste0(loaded_yaml$out_path, "y.csv")
    write.csv(cleaned$x, xpath)
    write.csv(cleaned$y, ypath)

    # update yaml parameters with new paths

    loaded_yaml$x_path = xpath
    loaded_yaml$y_path = ypath

    new_yaml_path = paste0(loaded_yaml$out_path, "yaml_parameters.yaml")
    # save this yaml file for reference
    yaml::write_yaml(loaded_yaml, new_yaml_path)

    return(list(x = cleaned$x, y = cleaned$y, yaml = new_yaml_path))
  } else {
    # error
    return()
  }
}
