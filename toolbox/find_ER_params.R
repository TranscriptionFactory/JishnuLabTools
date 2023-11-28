#!/usr/bin/env Rscript
library(tidyverse)
library(ggpubr)
library(caret)

###################################################
# path to yaml file (to load x and y, or just change below
###################################################
yaml_path = ''
yaml_args = yaml::yaml.load_file(yaml_path)

x = as.matrix(read.csv(yaml_args$x_path, row.names = 1, check.names = F))
y = as.matrix(read.csv(yaml_args$y_path, row.names = 1))

###################################################
# deltas and lambdas to test (bad parameters will cause the clustering to fail)
###################################################
deltas = c(0.01, 0.05, 0.1, 0.2)
lambdas = c(1.0, 0.1)



###################################################
# k-folds to test (remove nrow(x) for large n)
###################################################
kfolds = c(5, 10, nrow(x))


# code from test plainER with coeff partitions removed changed to look at delta/lambda
test_plainER_params = function(path_list = NULL, data_list = NULL, delta, lambda, k, thresh_fdr, row_samples = 3) {
  # sampling

  cat("\n ********************************************************************************")
  in_args = as.list(match.call())
  cat("\n using params delta =", delta, "lambda =", lambda, "k =", k, "thresh_fdr =", thresh_fdr, "\n")

  if ( !is.null(path_list) ) {
    x = as.matrix(read.csv(path_list[[1]], row.names = 1))
    y = as.matrix(read.csv(path_list[[2]], row.names = 1))

  } else if ( !is.null(data_list) ) {
    x = as.matrix(data_list[[1]])
    y = as.matrix(data_list[[2]])
  }

  # reorder df in ascending order of column coeff var
  full_x = x#[, order(col_coeff_var)]
  full_y = y

  replicate_results = data.frame()
  cat("\n ******************** Starting Replicates ******************** \n")

  for (row_sample_rep in 1:row_samples) {
    # for now, we just pick random sample of rows to use
    group_inds <- caret::createFolds(factor(y), k = k, list = TRUE, returnTrain = T)

    rowinds = group_inds[[sample(k)[1]]] #sample(1:nrow(sorted_x), floor(nrow(sorted_x)/k))

    x = full_x[rowinds, ]

    # x = apply(x, 2, function(x) scale(x, T, F))
    y = full_y[rowinds, ]

    n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions

    se_est <- apply(x, 2, stats::sd)

    sigma = cor(scale(x))

    control_fdr <- EssReg::threshSigma(x = x,
                                       sigma = sigma,
                                       thresh = thresh_fdr)
    sigma <- control_fdr$thresh_sigma
    kept_entries <- control_fdr$kept_entries

    result_AI <- EssReg::estAI(sigma = sigma, delta = delta, se_est = se_est)
    pure_numb <- sapply(result_AI$pure_list, FUN = function(x) {length(c(x$pos, x$neg))})

    A_hat <- result_AI$AI
    I_hat <- result_AI$pure_vec
    I_hat_list <- result_AI$pure_list

    if ( nrow(result_AI$AI) - length(result_AI$pure_vec) != ncol(sigma[, -result_AI$pure_vec]) ){
      replicate_result = "Failed"
      cat("\n \t\t Replicate", row_sample_rep, "of", row_samples, "Failed  \n")
    } else {
      replicate_result = "Passed"
      cat("\n \t\t Replicate", row_sample_rep, "of", row_samples, "Passed  \n")
    }

    replicate_results = rbind.data.frame(replicate_results,
                                         list(delta = delta, lambda = lambda, k = k,
                                              thresh_fdr = thresh_fdr,
                                              replicate_result = replicate_result,
                                              row_sample_replicate = row_sample_rep))

    next()
  }

  return(replicate_results)
}



find_plainER_params = function(deltas, lambdas, kfolds, ...) {
  res = data.frame()
  total_runs = length(deltas) * length(lambdas) * length(kfolds)
  for (d in deltas) {
    for (l in lambdas) {
      for (k in kfolds) {
        cat("\n runs remaining =", total_runs)
        res = rbind.data.frame(res, test_plainER_params(delta = d, lambda = l, k = k, ...))
        total_runs = total_runs - 1
      }
    }
  }
  res_summary = res %>% group_by(delta, lambda, k) %>% summarise(replicates_passed = length(which(replicate_result == "Passed")))

  res_summary$lambda = factor(res_summary$lambda)
  res_summary$delta = factor(res_summary$delta)
  # res_summary$replicates_passed = factor(res_summary$replicates_passed)

  kfold_label = paste0("k = ", res_summary$k)
  pl = ggpubr::ggdotchart(data = res_summary, x = "delta", y = "lambda", color = "replicates_passed", size = 6, facet.by = "k", nrow = length(kfolds),
                          panel.labs = list(k = unique(kfold_label)),
                          xlab = "Delta", ylab = "Lambda") +
    ggpubr::rotate_x_text(angle = 0) + ggplot2::scale_color_gradientn(colors = c("white", "red", "purple", "blue"), n.breaks = 3, limits = c(0,NA))

  return(list(pl = pl, res = res, res_summary = res_summary))
}




res_list = find_plainER_params(deltas = deltas, lambdas = lambdas, kfolds = kfolds, thresh_fdr = 0.2, data_list = list(x = x, y = y))




ggsave('/ix/djishnu/Aaron/0_for_others/Isha_ER_data/IM_orig_data.png', height = 2*length(kfolds), width = 4, plot = res_list$pl)

