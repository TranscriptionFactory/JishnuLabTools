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
