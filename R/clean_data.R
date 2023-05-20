#' @export
clean_data = function(xdata, ydata, edit_data = T,
                      col_var_quantile_filter = 0,
                      row_var_quantile_filter = 0,
                      remove_zero_median_cols = F,
                      remove_zero_median_rows = F,
                      scale_zeroes_directly = 0,
                      remove_spurious = c("GM42418", "LARS2"),
                      remove_mito_ribo = F,
                      er_input = NULL) {

  mode = ifelse(edit_data == T, 1, 0)

  y_correct_order <- match(rownames(xdata), rownames(ydata))

  # we may not have labeled y values, if not just print a message
  if (any(!is.na(y_correct_order))) {
    ydata <- ydata[y_correct_order]
    # fix ydata rownames
    ydata <- as.matrix(ydata)
    rownames(ydata) <- rownames(xdata)

  } else {
    cat("Y does not have labels. Assuming Y and X have rows in the same order")
  }

  if (mode > 0) {

    if (remove_mito_ribo) {
      xdata = JishnuLabTools::remove_mitochondrial_ribosomal_genes(xdata)
    }


    if (!is.null(remove_spurious) && length(remove_spurious) > 0) {
      spurious_cols = which(stringr::str_to_upper(colnames(xdata)) %in% remove_spurious)
      if (length(spurious_cols) > 0) {
        xdata = xdata[, -spurious_cols]
      }
    }

    if (scale_zeroes_directly > 0) {
      x_zeros = which(x == 0)
      x[x_zeros] = scale_zeroes_directly
    }

    # remove zero SD too
    zero_sd <- which(apply(xdata, 2, sd) == 0)

    # check if we need to remove any columns
    if (length(zero_sd) > 0) {
      xdata <- xdata[, -zero_sd]
    }

    if (remove_zero_median_cols) {

      zero_med_cols = which(apply(xdata, 2, median) == 0)

      if ( length(zero_med_cols) > 0 ) {
        xdata = xdata[, -zero_med_cols]
      }
    }

    if (col_var_quantile_filter > 0) {
      # filter column variance quantile
      colvars <- apply(xdata, 2, var)

      low_var_cols <- which(colvars < quantile(colvars, col_var_quantile_filter))

      if (length(low_var_cols > 0)) {
        xdata <- xdata[, -low_var_cols]
      }
    }

    if (remove_zero_median_rows) {
      # remove empty rows last
      empty_rows <- which(apply(xdata, 1, median) == 0)

      # check if we need to remove any rows
      if (length(empty_rows) > 0) {
        # remove rows
        xdata <- xdata[-empty_rows,]
        ydata <- ydata[-empty_rows]
      }
    }

    if (row_var_quantile_filter) {
      # remove empty rows last
      rowvars <- apply(xdata, 1, var)
      low_var_rows <- which(rowvars < quantile(rowvars, row_var_quantile_filter))

      # check if we need to remove any rows
      if (length(low_var_rows) > 0) {
        # remove rows
        xdata <- xdata[-low_var_rows,]
        ydata <- ydata[-low_var_rows]
      }
    }
  }

  return(list("x" = xdata, "y" = ydata))
}
