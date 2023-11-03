#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(EssReg)
library(doParallel)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

yaml_path = args[1]


ER_error_helper = function(x = NULL, y = NULL, yaml_args = NULL, yaml_path = NULL,
                           error_file_output_path = NULL) {


  results_log = list()

  if (!is.null(yaml_path) & file.exists(yaml_path)) {
    yaml_args = yaml::yaml.load_file(yaml_path)
    yaml_x = as.matrix(read.csv(yaml_args$x_path, row.names = 1))
    yaml_y = as.matrix(read.csv(yaml_args$y_path, row.names = 1))

    error_file_output_path = yaml_args$out_path
  } else if (is.null(yaml_args)) {
    yaml_args = list()
  }

  if (all(is.null(yaml_args), is.null(x), is.null(y))) {
    cat("\n X,Y data not loaded properly from yaml or passed arguments \n")
    stop()
  }

  if (is.null(x)) {
    # use passed x with preference; only use yaml x if no other option
    x = yaml_x
  }

  if (is.null(y)) {
    # use passed y with preference; only use yaml y if no other option
    y = yaml_y
  }



  check_data = function(x, y) {

    error_messages = c()
    data_messages = c()
    data = as.matrix(apply(x, 2, as.numeric))

    if (any(is.na(data))) {
      error_messages = append(error_messages, paste0("\nNAs in data at", which(is.na(data)), "\n"))
    }

    if (any(is.nan(data))) {
      error_messages = append(error_messages, paste0("\nNaNss in data at", which(is.nan(data)), "\n"))
    }

    if (any(apply(data, 2, sd) == 0)) {
      error_messages = append(error_messages, "\n Columns with sd=0 at", which(apply(data,2,sd)==0), "\n")
    }

    col_variance = apply(data, 2, var)
    col_var_hist = hist(col_variance, xlab = "Column-wise Variance")
    # data_messages = c(data_messages, paste0("Column-wise variance :", col_var_hist$density[1]/sum(col_var_hist$density), "% of variables have a variance below ", col_var_hist$breaks[2], "\n"))
    
    col_cv = apply(data, 2, function(x) sd(x)/mean(x))
    col_cv_hist = hist(col_cv, xlab = "Column-wise Coef. Variation")
    # data_messages = c(data_messages, paste0("Column-wise variance :", col_var_hist$density[1]/sum(col_var_hist$density), "% of variables have a variance below ", col_var_hist$breaks[2], "\n"))

    if (is.null(error_messages)) {
      error_messages = c("\n No Data Errors \n")
    }

    return(list(error_messages = error_messages,
                # data_messages = data_messages,
                col_var_hist = col_var_hist,
                col_cv_hist = col_cv_hist))

  }

  # check data
  # load x,y
  results_log$check_data = check_data(x, y)

  cat(results_log$check_data$error_messages)

  check_initial_ER_steps = function(yaml_args = NULL, yaml_path = NULL, x = NULL, y = NULL) {

    if (is.null(yaml_args)) {
      yaml_args = yaml::yaml.load_file(yaml_path)
    }

    if (is.null(x)) {
      x = as.matrix(read.csv(yaml_args$x_path, row.names = 1))
    }
    if (is.null(y)) {
      y = as.matrix(read.csv(yaml_args$y_path, row.names = 1))
    }

    n <- nrow(x);  p <- ncol(x)
    unscaled_x = x
    x = scale(x, T, T)

    yaml_args$thresh_fdr = ifelse(is.null(yaml_args$thresh_fdr), 0.05, yaml_args$thresh_fdr)
    yaml_args$alpha_level = ifelse(is.null(yaml_args$alpha_level), 0.05, yaml_args$alpha_level)

    yaml_args$delta = ifelse(is.null(yaml_args$delta), 0.1, yaml_args$delta)
    yaml_args$lambda = ifelse(is.null(yaml_args$lambda), 1.0, yaml_args$lambda)

    sigma = cor(x)
    se_est = apply(x, 2, sd)

    control_fdr = EssReg::threshSigma(x, sigma, yaml_args$thresh_fdr)

    delta = ifelse(length(yaml_args$delta) > 1, min(yaml_args$delta), yaml_args$delta)
    lambda = ifelse(length(yaml_args$lambda) > 1, min(yaml_args$lambda), yaml_args$lambda)

    result_AI = EssReg::estAI(sigma, delta, se_est)

    pure_numb = sapply(result_AI$pure_list, FUN = function(x) {length(c(x$pos, x$neg))})

    A_hat <- result_AI$AI
    I_hat <- result_AI$pure_vec
    I_hat_list <- result_AI$pure_list

    if (is.null(I_hat)) {
      cat("Algorithm fails due to non-existence of pure variable.\n")
      stop()
    }

    C_hat <- EssReg::estC(sigma = sigma, AI = A_hat)

    Gamma_hat <- rep(0, p)
    Gamma_hat[I_hat] <- diag(sigma[I_hat, I_hat]) - diag(A_hat[I_hat, ] %*% C_hat %*% t(A_hat[I_hat, ]))
    Gamma_hat[Gamma_hat < 0] <- 1e-2

    pred_result <- EssReg::prediction(y = y, x = x, sigma = sigma, A_hat = A_hat,
                                      Gamma_hat = Gamma_hat, I_hat = I_hat)

    theta_hat <- pred_result$theta_hat

    Q <- try(theta_hat %*% solve(crossprod(x %*% theta_hat) / n, crossprod(theta_hat)), silent = T)
    if (class(Q)[1] == "try-error") {
      Q <- theta_hat %*% MASS::ginv(crossprod(x %*% theta_hat) / n) %*% crossprod(theta_hat)
    }

    if (length(result_AI$pure_vec) != nrow(sigma)) { ## check if all vars are pure vars?
      sigma_TJ <- EssReg::estSigmaTJ(sigma = sigma, AI = A_hat, pure_vec = result_AI$pure_vec)

      AI_hat <- abs(A_hat[I_hat, ]) ## just rows of pure variables
      sigma_bar_sup <- max(solve(crossprod(AI_hat), t(AI_hat)) %*% se_est[I_hat]) ## not sure what this does
      AJ <- EssReg::estAJDant(C_hat = C_hat, sigma_TJ = sigma_TJ,
                              lambda = lambda * delta * sigma_bar_sup,
                              se_est_J = sigma_bar_sup + se_est[-I_hat])

      if (is.null(AJ)) {
        cat("\n Dantzig estimation failed\n")
      } else {
        AJ <- t(solve(C_hat, sigma_TJ))
      }
      A_hat[-result_AI$pure_vec, ] <- AJ
    }


    Gamma_hat[-I_hat] <- diag(sigma[-I_hat, -I_hat]) - diag(A_hat[-I_hat,] %*% C_hat %*% t(A_hat[-I_hat, ]))
    Gamma_hat[Gamma_hat < 0] <- 1e2


    res_beta <- EssReg::estBeta(y = y, x = x, sigma = sigma, A_hat = A_hat,
                                C_hat = C_hat, Gamma_hat = Gamma_hat, I_hat = I_hat,
                                I_hat_list = I_hat_list, alpha_level = yaml_args$alpha_level)
    beta_hat <- res_beta$beta_hat
    beta_conf_int <- res_beta$conf_int
    beta_var <- res_beta$beta_var
    rownames(A_hat) <- colnames(x)
    colnames(A_hat) <- paste0("Z", 1:ncol(A_hat))

    rownames(C_hat) <- colnames(C_hat) <- paste0("Z", 1:ncol(C_hat))

    I_clust <- NULL
    for (i in 1:length(I_hat_list)) {
      clust_name <- paste0("Z", i)
      cluster <- I_hat_list[[i]]
      pos <- cluster$pos
      neg <- cluster$neg
      if (length(pos) > 0) {
        names(pos) <- colnames(x)[pos]
      } else {
        pos <- NULL
      }
      if (length(neg) > 0) {
        names(neg) <- colnames(x)[neg]
      } else {
        neg <- NULL
      }
      I_clust[[clust_name]] <- list("pos" = pos,
                                    "neg" = neg)
    }

    names(I_hat) <- colnames(x)[I_hat]

    cat("\nFinished plainER\n")
    return(list(K = ncol(A_hat),
                A = A_hat,
                C = C_hat,
                I = I_hat,
                I_clust = I_clust, ## original is I_ind
                Gamma = Gamma_hat,
                beta = beta_hat,
                beta_conf_int = beta_conf_int, ## original is beta_CIs
                beta_var = beta_var,
                pred = pred_result,
                opt_lambda = lambda,
                opt_delta = delta / sqrt(log(max(p, n)) / n),
                Q = Q,
                thresh_sigma = sigma))

  }

 cat("\nStarting plainER \n")

  withCallingHandlers(expr = (er_res <<- check_initial_ER_steps(yaml_args)),
                                            error = function(e) {
                                              print(sys.calls())
                                            })

  results_log$plainER = er_res

  cat("\nCompleted plainER (no cross validation). Saving results in ", error_file_output_path, "\n")

  if (is.null(error_file_output_path)) {
    error_file_output_path = getwd()
  } else if (!dir.exists(error_file_output_path)) {
    dir.create(error_file_output_path, recursive = T)
  }

  error_file_output_path = paste0(error_file_output_path, "/ER_results_log.RDS")

  cat("\nSaving results log to ", error_file_output_path, "\n")
  saveRDS(results_log, error_file_output_path)

  cat("\nStarting plainER with cross validation \n")

  pipeline3_output_path = paste0(error_file_output_path, "/pipeline3_error_stacktrace.txt")

  cat("\nStarting pipeline3. Will output pipeline3 error messages to ", pipeline3_output_path, "\n")

  withCallingHandlers(EssReg::pipelineER3(yaml_path),
                      error = function(e) {
                        cat(sys.calls(), file = pipeline3_output_path, append = T)
                      })

  return(results_log)
}

ER_error_helper(yaml_path = yaml_path)
