
#' @export
correlation_filter = function(data, correlation_threshold = 0.8,
                              correlation_quantile_filter = 0.9) {

  # scale data
  scaled_df = apply(data, 2, function(x) scale(x, T, F))

  # get unsigned correlation matrix
  df_cor = abs(cor(as.matrix(scaled_df)))

  # make upper triangular
  df_cor[lower.tri(df_cor, diag = T)] = 0

  # for each column, figure out how many variables are correlated with this variable
  # e.g. for correlation_threshold = 0.8, how many variables have a correlation > 0.8
  # with this variable
  df_cor_count = apply(df_cor, 2, function(x) length(which(x > 0.8)))

  # from this, filter out the top quantile (variables with most correlations)
  # from the original dataframe
  df_filtered = data[, -which(df_cor_count > quantile(df_cor_count, correlation_quantile_filter))]

  return(df_filtered)
}

