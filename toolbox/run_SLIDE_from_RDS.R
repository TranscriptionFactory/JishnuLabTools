#!/usr/bin/env Rscript
library(tidyverse)
library(devtools)
library(EssReg)
library(SLIDEHelper)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

# in the "..." in the function parameters, pass the filtering parameters you want to use
# for clean data
# clean_data = function(xdata, ydata, edit_data = T,
#                       col_var_quantile_filter = 0,
#                       row_var_quantile_filter = 0,
#                       col_sparsity_min_nonzero = 0,
#                       row_sparsity_min_nonzero = 0,
#                       remove_zero_median_cols = F,
#                       remove_zero_median_rows = F,
#                       scale_zeroes_directly = 0,
#                       remove_spurious = c("GM42418", "LARS2"),
#                       remove_mito_ribo = F,
#                       er_input = NULL)

run_SLIDE_from_RDS(data_frame_path = NULL,
                   data_frame = NULL, y_col_name = NULL,
                   y_col_position = 1,
                   output_path = NULL,
                   lambdas = c(0.1, 1.0),
                   deltas = c(0.1, 0.01),
                   k = 10, y_factor = TRUE, y_levels = c(0, 1),
                   eval_type = "cor", rep_cv = 25, n_reps = 25,
                   alpha_level = 0.05, thresh_fdr = 0.2,
                   permute = TRUE,
                   std_cv = FALSE,
                   std_y = FALSE,
                   benchmark = FALSE,
                   lasso = TRUE, plsr = TRUE, pcr = TRUE, ...) {

  # get the arguments we passed to this function
  yaml_args = as.list(match.call())

  # add a "/" to output path
  yaml_args$out_path = paste0(yaml_args$out_path, "/")

  # cleaning data function
  clean_data = function(xdata, ydata, edit_data = T,
                        col_var_quantile_filter = 0,
                        row_var_quantile_filter = 0,
                        col_sparsity_min_nonzero = 0,
                        row_sparsity_min_nonzero = 0,
                        remove_zero_median_cols = F,
                        remove_zero_median_rows = F,
                        scale_zeroes_directly = 0,
                        remove_spurious = c("GM42418", "LARS2"),
                        remove_mito_ribo = F,
                        er_input = NULL) {

    mode = ifelse(edit_data == T, 1, 0)

    # replace all NAs with zeros
    xdata[is.na(xdata)] = 0

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

      if (col_sparsity_min_nonzero == 0 && remove_zero_median_cols == T) {
        col_sparsity_min_nonzero = 0.5
      }

      if (row_sparsity_min_nonzero == 0 && remove_zero_median_rows == T) {
        row_sparsity_min_nonzero = 0.5
      }

      if (remove_mito_ribo) {

        remove_mitochondrial_ribosomal_genes = function(df) {

          # check if genes have cluster names e.g C1.RPS prepended
          cluster_names = which(grepl("^[[:alnum:]].\\d+.RPL|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", stringr::str_to_upper(colnames(df))) == T)

          no_cluster_names = which(grepl("^RPL|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", stringr::str_to_upper(colnames(df))) == T)

          if (length(cluster_names) == 0 && length(no_cluster_names) == 0) {
            return(df)
          } else if (length(cluster_names) > 0) {
            df = df[, -cluster_names]
          } else {
            df = df[, -no_cluster_names]
          }
          return(df)
        }

        xdata = remove_mitochondrial_ribosomal_genes(xdata)
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

      if (col_sparsity_min_nonzero > 0) {

        col_nonzero = which(apply(xdata, 2, function(x) length(which(x != 0)) / length(x)) > col_sparsity_min_nonzero)

        if ( length(col_nonzero) > 0 ) {
          xdata = xdata[, -col_nonzero]
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

      if (row_sparsity_min_nonzero > 0) {
        # remove empty rows last
        row_nonzero = which(apply(xdata, 1, function(x) length(which(x != 0)) / length(x)) > row_sparsity_min_nonzero)

        # check if we need to remove any rows
        if (length(row_nonzero) > 0) {
          # remove rows
          xdata <- xdata[-row_nonzero,]
          ydata <- ydata[-row_nonzero]
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


  # check if we have a dataframe or path to dataframe
  if (!is.null(data_frame_path)) {

    if (stringr::str_detect(str_to_lower(data_frame_path), pattern = ".csv")) {
      data_frame = as.matrix(read.csv(data_frame_path, row.names = 1))
    } else if (stringr::str_detect(str_to_lower(data_frame_path), pattern = ".rds")) {
      data_frame = readRDS()
    }
  }

  # set output path
  out_path = yaml_args$out_path

  cat("Saving output to ", out_path, "\n")

  if (!is.null(y_col_name) && is.character(y_col_name) && y_col_name %in% colnames(data_frame)) {
    x = data_frame[, colnames(data_frame) != y_col_name]
    y = data_frame[, y_col_name]
  } else if (!is.null(y_col_position) && is.numeric(y_col_position) &&
             y_col_position > 1 && y_col_position < ncol(data_frame)) {
    x = data_frame[, -c(y_col_position)]
    y = data_frame[, y_col_position]
  } else {
    cat("Error in selecting Y by position or name in dataframe. Check Y selection \n")
    return()
  }

  # check whether our x is continous or binary
  if (length(unique(y[, 1])) == 2) {
    yaml_args$eval_type = "auc"
  }

  rownames(x) = y[, 1]

  sampling_indices = sample(1:nrow(x))

  # shuffle
  x = x[sampling_indices, ]
  y = y[sampling_indices, 1]

  # Begin filtering data
  cat("Dim df is pre-filter ", dim(x), "\n")

  # don't quantile filter columns if there aren't many
  data = clean_data(xdata = x, ydata = y,
                                    edit_data = T, ...)

  cat("Dim df is post-filter ", dim(x), "\n")

  # save our new x and y
  x = data$x
  y = data$y

  yaml_args$x_path = paste0(yaml_args$out_path, "x.csv")
  yaml_args$y_path = paste0(yaml_args$out_path, "y.csv")
  write.csv(x, yaml_args$x_path)
  write.csv(y, yaml_args$y_path)

  # write yaml file to output
  yaml_path = paste0(yaml_args$out_path, "yaml_input.yaml")
  yaml::write_yaml(yaml_args, yaml_path)

  if (yaml_args$k != 10) {
    # don't do anything if our is k is the default
    if (nrow(x) < 20) {
      yaml_args$k = nrow(x)
    } else if (nrow(x) < 50) {
      yaml_args$k = 5
    } else {
      yaml_args$k = 10
    }
  }

  if (all(lambdas > 0)) {
    lbds = lambdas
  } else {
    cat("Lambdas need to be greater than 0. Setting to 0.1 and 1.0 \n")
    lbds = c(0.1, 1.0)
  }

  if (!all(deltas > 0)) {
    dlts = deltas
  } else {
    cat("Deltas need to be greater than 0. Setting to 0.01 and 0.1 \n")
   dlts = c(0.01, 0.1)
  }

  orig_path = yaml_args$out_path

  run_count = 0
  for (ind1 in 1:length(lbds)) {
    for (ind2 in 1:length(dlts)) {

      run_count = run_count + 1
      er_input =  yaml::yaml.load_file(yaml_path)

      er_input$out_path = orig_path

      er_input$lambda = lbds[ind1]
      er_input$delta = dlts[ind2]

      er_input$out_path = stringr::str_c(er_input$out_path, "run", run_count, "_lambda", lbds[ind1], "_delta", dlts[ind2], "/", sep = "")
      yaml_path = paste0(er_input$out_path, "/yaml_input.yaml")
      yaml::write_yaml(er_input, yaml_path)

      cat("running pipeline 3 with reps \n")

      EssReg::pipelineER3(yaml_path)

      # check for ER results in output folder
      er_results_path = list.files(er_input$out_path, pattern = "final_delta")

      if (length(er_results_path == 0)) {
        cat("ER did not complete successfully, trying next delta/lambda combination \n")
      } else {

        Z_matrix <- SLIDEHelper::CalcZMatrix(x_path, er_path, out_path)

        # run slide on these results
        SLIDE_res = SLIDEHelper::runSLIDE(y_path = yaml_args$y_path,
                              z_path = NULL,
                              z_matrix = Z_matrix,
                              er_path = er_results_path,
                              do_interacts = TRUE,
                              spec = 0.1,
                              niter = 100)
        num_top_feats <- 10
        condition <- yaml_args$eval_type

        SLIDE_res <- GetTopFeatures(yaml_args$x_path, yaml_args$y_path, er_path, out_path,
                                    SLIDE_res, num_top_feats = 10, condition)

        plotSigGenes(SLIDE_res, er_input$out_path, plot_interaction = TRUE)

        #the SLIDE_res has to be the output from GetTopFeatures
        CalcControlPerformance(z_matrix = Z_matrix, y_path, SLIDE_res, niter = 1000, condition, er_input$out_path)
      }
        registerDoParallel(cores)
    }
  }
}
