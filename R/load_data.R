#' @export
# load data as a single matrix with first column as Y values
# out_path is from yaml file
load_data = function(obj, loaded_yaml, remove_mito_ribo = F, create_dir = T) {

  if (! dir.exists(loaded_yaml$out_path)) {
    if (create_dir) {
      dir.create(loaded_yaml$out_path)
    } else {
      cat("Directory doesn't exist. Check yaml out_path or set create_dir = T \n")
      return()
    }
  }

  df = JishnuLabTools:::check_for_df_or_path(obj)

  if (remove_mito_ribo) {
    df[, -1] = JishnuLabTools::remove_mitochondrial_ribosomal_genes(df[, -1])
  }

  if (! is.null(df)) {
    # shuffle rows
    df = df[sample(1:nrow(df)),]

    x = as.matrix(df[, -1])
    y = as.matrix(df[, 1])

    # will just remove columns with zero standard deviation.
    # to remove based on quantile, set edit_data = T and set quantile
    cleaned = JishnuLabTools::clean_data(x = x, y = y, edit_data = F,
                                         quantile_filter = 0, remove_zero_median_cols = F)


    xpath = paste0(loaded_yaml$out_path, "x.csv")
    ypath = paste0(loaded_yaml$out_path, "y.csv")
    write.csv(cleaned$x, xpath)
    write.csv(cleaned$y, ypath)

    # update yaml parameters with new paths

    loaded_yaml$x_path = xpath
    loaded_path$y_path = ypath

    # save this yaml file for reference
    yaml::write_yaml(loaded_yaml, loaded_yaml$out_path)

    return(list("x" = cleaned$x, "y" = cleaned$y))
  } else {
    # error
    return()
  }
}
