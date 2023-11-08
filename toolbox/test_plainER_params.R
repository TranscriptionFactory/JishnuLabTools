library(tidyverse)
library(ggpubr)
library(caret)

# code from test plainER with coeff partitions removed changed to look at delta/lambda
test_plainER_params = function(path_list = NULL, data_list = NULL, delta, lambda, k, thresh_fdr, row_samples = 3) {
  # sampling

  cat("\n ******************************************************************************** \n")
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

    rowinds = group_inds[[1]] #sample(1:nrow(sorted_x), floor(nrow(sorted_x)/k))

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

    # C_hat <- EssReg::estC(sigma = sigma, AI = A_hat)
    #
    # Gamma_hat <- rep(0, p)
    # Gamma_hat[I_hat] <- diag(sigma[I_hat, I_hat]) - diag(A_hat[I_hat, ] %*% C_hat %*% t(A_hat[I_hat, ]))
    # Gamma_hat[Gamma_hat < 0] <- 1e-2 #### replace negative values with something very close to 0
    #
    #
    # pred_result <- EssReg::prediction(y = y, x = x, sigma = sigma, A_hat = A_hat,
    #                                   Gamma_hat = Gamma_hat, I_hat = I_hat)
    #
    # #### theta_hat (supplement 2.2)
    # theta_hat <- pred_result$theta_hat
    # #### Inference In Latent Factor Regression With Clusterable Features
    # #### Z_tilde = Q*X
    #
    # Q <- try(theta_hat %*% solve(crossprod(x %*% theta_hat) / n, crossprod(theta_hat)), silent = T)
    # if (class(Q)[1] == "try-error") {
    #   Q <- theta_hat %*% MASS::ginv(crossprod(x %*% theta_hat) / n) %*% crossprod(theta_hat)
    # }
    #
    #
    #
    # if (length(result_AI$pure_vec) != nrow(sigma)) { ## check if all vars are pure vars?
    #   sigma_TJ <- EssReg::estSigmaTJ(sigma = sigma, AI = A_hat, pure_vec = result_AI$pure_vec)
    # }
    # cat("\n Completed \n")
  }

  return(replicate_results)
}
