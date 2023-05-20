#' @export
clean_data <- function(xdata, ydata, edit_data = F,
                      col_var_quantile_filter = 0, remove_zero_median_cols = F,
                      remove_zero_median_rows = F, scale_zeroes_directly = 0,
                      remove_spurious = c("GM42418", "LARS2"),
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

    if (col_var_quantile_filter >= 0) {
      # filter quantiles and remove empty rows

      empty_cols <- which(apply(xdata, 2, median) == 0)

      colvars <- apply(xdata, 2, var)

      low_var_cols <- which(colvars < quantile(colvars, col_var_quantile_filter))

      if (remove_zero_median_cols) {
        remove_cols <- base::union(empty_cols, low_var_cols)
      } else {
        remove_cols <- low_var_cols
      }

      if (length(remove_cols > 0)) {
        xdata <- xdata[, -remove_cols]
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
    }
  }

  return(list("x" = xdata, "y" = ydata))
}
