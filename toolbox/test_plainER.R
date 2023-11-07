

test_plainER = function(path_list = NULL, data_list = NULL, delta, k, thresh_fdr) {
  # sampling

  if ( !is.null(path_list) ) {
    x = as.matrix(read.csv(path_list[[1]], row.names = 1))
    y = as.matrix(read.csv(path_list[[2]], row.names = 1))

  } else if ( !is.null(data_list) ) {
    x = as.matrix(data_list[[1]])
    y = as.matrix(data_list[[2]])
  }

  sinds = sample(1:nrow(x), size = floor(nrow(x)/k))
  x = x[sinds, ]
  y = y[sinds]

  n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions

  col_var = apply(x, 2, var)

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
    cat("\n error with clustering \n")
    return()
  } else {
    cat("\n dims are equal\n")
  }

  C_hat <- EssReg::estC(sigma = sigma, AI = A_hat)

  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(sigma[I_hat, I_hat]) - diag(A_hat[I_hat, ] %*% C_hat %*% t(A_hat[I_hat, ]))
  Gamma_hat[Gamma_hat < 0] <- 1e-2 #### replace negative values with something very close to 0


  pred_result <- EssReg::prediction(y = y, x = x, sigma = sigma, A_hat = A_hat,
                            Gamma_hat = Gamma_hat, I_hat = I_hat)

  #### theta_hat (supplement 2.2)
  theta_hat <- pred_result$theta_hat
  #### Inference In Latent Factor Regression With Clusterable Features
  #### Z_tilde = Q*X

  Q <- try(theta_hat %*% solve(crossprod(x %*% theta_hat) / n, crossprod(theta_hat)), silent = T)
  if (class(Q)[1] == "try-error") {
    Q <- theta_hat %*% MASS::ginv(crossprod(x %*% theta_hat) / n) %*% crossprod(theta_hat)
  }



  if (length(result_AI$pure_vec) != nrow(sigma)) { ## check if all vars are pure vars?
    sigma_TJ <- EssReg::estSigmaTJ(sigma = sigma, AI = A_hat, pure_vec = result_AI$pure_vec)
  }
  cat("\n Completed \n")
}






kfold_plainER_steps = function(k, ...) {
  cpu_times = list()
  for (i in 1:k) {
    cpu_times[[k]] = withCallingHandlers(system.time(test_plainER(k = i, ...)), error = function(e){ print(sys.calls())})
  }
  return(cpu_times)
}



kfold_plainER_steps(k = 2, delta = 0.1, thresh_fdr = 0.2, data_list = list(x = x, y = y))


