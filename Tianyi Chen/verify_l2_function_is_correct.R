 y = df_plot$y[1:30]
 tmax <- length(y)
 x <- 1:tmax
 best_cp <- NULL
 min_loss <- Inf
 best_coefs <- NULL
   cp = 3
   X <- cbind(1, (x - cp), (x > cp) * (x - cp))
   # Fit the model using least squares
   fit <- lm(y ~ X - 1)  # "-1" removes intercept as it's already in design matrix
solve(t(X) %*% (X)) %*% t(X) %*%y 
Px = X%*%solve(t(X) %*% (X)) %*% t(X) 
t(y) %*% (diag(1,tmax) - Px)  %*% y
# Get the coefficients and ensure beta_L != beta_R
coef_fit <- coef(fit)
alpha_hat <- coef_fit[1]
beta_L_hat <- coef_fit[2]
beta_R_hat <- coef_fit[3] + beta_L_hat  # Adjust for slope change

if (beta_L_hat != beta_R_hat) {
  # Calculate sum of squared residuals
  fitted_values <- alpha_hat + beta_L_hat * (x - cp) + (beta_R_hat - beta_L_hat) * (x - cp) * (x > cp)
  loss <- sum((y - fitted_values)^2)
  t(y) %*% (diag(1,tmax) - Px)  %*% y
  
  sum((y - fitted_values)^2)
  # Update best changepoint if this loss is the minimum
  if (loss < min_loss) {
    min_loss <- loss
    best_cp <- cp
    best_coefs <- c(alpha = alpha_hat, beta_L = beta_L_hat, beta_R = beta_R_hat)
  }
}
}