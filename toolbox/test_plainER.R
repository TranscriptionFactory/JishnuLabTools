library(ggpubr)

test_plainER = function(path_list = NULL, data_list = NULL, delta, k, thresh_fdr,
                        coeff_var_partitions = 3, row_samples = 3, partition_random_sample_percent = 0.5) {
  # sampling

  cat("\n ******************************************************************************** \n")
  in_args = as.list(match.call())
  cat("\n using params delta =", delta, " k =", k, " thresh_fdr =", thresh_fdr, "\n")

  if ( !is.null(path_list) ) {
    x = as.matrix(read.csv(path_list[[1]], row.names = 1))
    y = as.matrix(read.csv(path_list[[2]], row.names = 1))

  } else if ( !is.null(data_list) ) {
    x = as.matrix(data_list[[1]])
    y = as.matrix(data_list[[2]])
  }

  cat("\n Splitting data into low, medium and high coeff. var groups to sample from \n")

  col_coeff_var = apply(apply(x, 2, function(x) scale(x, T, F)), 2, function(x) sd(x)/mean(x))

  # reorder df in ascending order of column coeff var
  sorted_x = x[, order(col_coeff_var)]
  true_y = y

  partition_results = data.frame()
  partitions = seq(0, 1, 1/coeff_var_partitions)


  for (row_sample_rep in 1:row_samples) {
    # for now, we just pick random sample of rows to use
    rowinds = sample(1:nrow(sorted_x), floor(nrow(sorted_x)/k))

    partition_name = c("Low", "Medium", "High", "From All")
    for (partition_index in 1:(length(partitions))) {

      # select our quantile
      lower_quantile = quantile(col_coeff_var, partitions[partition_index])
      upper_quantile = quantile(col_coeff_var, partitions[partition_index + 1])

      cat("\n ******************** Starting partition ", partition_index, " Replicate ", row_sample_rep, " of ", row_samples, " ******************** \n")

      if (partition_index == 4) {
        coeff_var_partition_to_sample_from = sample(1:ncol(sorted_x), floor(ncol(sorted_x)/3))
      } else {
        coeff_var_partition_to_sample_from = which(col_coeff_var > lower_quantile & col_coeff_var < upper_quantile)

        # take 25% out of this and replace with random from rest of df

        coeff_var_partition_to_sample_from = coeff_var_partition_to_sample_from[-sample(coeff_var_partition_to_sample_from,
                                                                                        floor(length(coeff_var_partition_to_sample_from) * partition_random_sample_percent))]

        coeff_var_partition_to_sample_from = base::union(coeff_var_partition_to_sample_from, sample(1:ncol(sorted_x), floor(length(coeff_var_partition_to_sample_from) * partition_random_sample_percent)))
      }
      sinds = sample(coeff_var_partition_to_sample_from)
      x = sorted_x[rowinds, sinds]
      y = true_y[rowinds]

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
        partition_result = "Failed"
        cat("\n \t\t\t Partition ", partition_index, " (", partition_name[partition_index], "Coeff Columns ) Failed  \n")
      } else {
        partition_result = "Passed"
        cat("\n \t\t\t Partition ", partition_index, " (", partition_name[partition_index], "Coeff Columns ) Passed  \n")
      }

      partition_results = rbind.data.frame(partition_results,
                                           list(delta = delta, k = k,
                                           partition = partition_name[partition_index], partition_result = partition_result,
                                           thresh_fdr = thresh_fdr,
                                           row_sample_replicate = row_sample_rep))

      next()

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
  }
  return(partition_results)
}






# kfold_plainER_steps = function(k, ...) {
#   cpu_times = list()
#   for (i in 1:k) {
#     cpu_times[[k]] = withCallingHandlers(system.time(test_plainER(k = i, ...)), error = function(e){ print(sys.calls())})
#   }
#   return(cpu_times)
# }
kfold_plainER_steps = function(k, delta,...) {
  res = data.frame()
  for (i in 1:k) {
    for (d in delta) {
      res = rbind.data.frame(res, test_plainER(k = i, delta = d, ...))
    }
  }
  return(res)
}

# delta = c(0.01, 0.05, 0.1)
res = kfold_plainER_steps(k = 5, delta = c(0.01, 0.05, 0.1), thresh_fdr = 0.2, data_list = list(x = x, y = y))

ggpubr::ggdotchart(data = res, x = "k", y = "partition", color = "partition_result", size = 5, facet.by = "delta", nrow = 3,
                   xlab = "Fold Size (k); data split into nrow(x)/k", ylab = "Coeff. Var Partition (3 Percentiles)") + ggpubr::rotate_x_text(angle = 0)
# ggpubr::ggdotchart(data = res, x = "k", y = "partition", color = "partition_result", size = 5,
#                    xlab = "Fold Size (k); data split into nrow(x)/k", ylab = "Coeff. Var Partition (3 Percentiles)") + ggpubr::rotate_x_text(angle = 0)


ggsave('/ix/djishnu/Aaron/0_for_others/Isha_ER_data/IM_highvar_data.png', height = 5, width = 5)

