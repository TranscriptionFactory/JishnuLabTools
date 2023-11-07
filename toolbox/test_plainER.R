



# sampling

sinds = sample(1:nrow(x), size = floor(nrow(x)/2))
x = x[sinds, ]
y = y[sinds]


n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions

col_var = apply(x, 2, var)

se_est <- apply(x, 2, stats::sd)

sigma = cor(scale(x))

result_AI <- EssReg::estAI(sigma = sigma, delta = 0.001, se_est = se_est)
pure_numb <- sapply(result_AI$pure_list, FUN = function(x) {length(c(x$pos, x$neg))})

A_hat <- result_AI$AI
I_hat <- result_AI$pure_vec
I_hat_list <- result_AI$pure_list

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
