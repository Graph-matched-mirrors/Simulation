find_slope_changepoint_with_plot <- function(y, doplot = TRUE) {
  tmax <- length(y)
  x <- 1:tmax
  best_cp <- NULL
  min_loss <- Inf
  best_coefs <- NULL
  
  for (cp in 2:(tmax - 1)) {
    # Construct design matrix based on the given model
    X <- cbind(1, (x - cp), (x > cp) * (x - cp))
    
    # Fit the model using least squares
    fit <- lm(y ~ X - 1)  # "-1" removes intercept as it's already in design matrix
    
    # Get the coefficients and ensure beta_L != beta_R
    coef_fit <- coef(fit)
    alpha_hat <- coef_fit[1]
    beta_L_hat <- coef_fit[2]
    beta_R_hat <- coef_fit[3] + beta_L_hat  # Adjust for slope change
    
    if (beta_L_hat != beta_R_hat) {
      # Calculate sum of squared residuals
      fitted_values <- alpha_hat + beta_L_hat * (x - cp) + (beta_R_hat - beta_L_hat) * (x - cp) * (x > cp)
      loss <- sum((y - fitted_values)^2)
      
      # Update best changepoint if this loss is the minimum
      if (loss < min_loss) {
        min_loss <- loss
        best_cp <- cp
        best_coefs <- c(alpha = alpha_hat, beta_L = beta_L_hat, beta_R = beta_R_hat)
      }
    }
  }
  
  best_coefs
  # Generate fitted values using the best coefficients
  if (!is.null(best_coefs)) {
    fitted_y <- best_coefs["alpha.X1"] + best_coefs["beta_L.X2"] * (x - best_cp) + 
      (best_coefs["beta_R.X3"] - best_coefs["beta_L.X2"]) * (x - best_cp) * (x > best_cp)
  } else {
    fitted_y <- rep(NA, tmax)  # Return NA values if no changepoint was found
  }
  
  # Plot the results if doplot is TRUE
  if (doplot && !is.null(best_cp)) {
    plot(x, y, pch = 16, col = "black", main = "", xlab = "time", ylab = "mirror")
    lines(x, fitted_y, col = "red", lwd = 2)
    abline(v = best_cp, col = "red", lwd = 2, lty = 2)
    abline(v = tstar, col = "black", lwd = 2, lty = 2)
    legend("topleft", legend = c("Data", "Fitted Line l2", "Estimated_CP l2"),
           col = c("black", "red", "red"), pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = 2)
  }
  
  # Return results
  return(list(changepoint = best_cp, coefficients = best_coefs , error = (best_cp- tstar)/tmax   ))
}

# Example usage with plotting
y <- c(1, 2, 3, 10, 12, 14, 20, 22, 24)
result_with_plot <- find_slope_changepoint_with_plot(y, doplot = TRUE)
print(result_with_plot)

# Example usage without plotting
result_without_plot <- find_slope_changepoint_with_plot(y, doplot = FALSE)
print(result_without_plot)


result_without_plot$changepoint


set.seed(2)
n = 200
p <- 0.3
q <- 0.4
nmc = 10

# Define the values of mm
mm <- c(40)

# Loop over each value of m in mm
for (i in 1:length(mm)) {
  tmax <- m <- mm[i]
  tstar <- tmax / 2
  delta <- (1-0.1)/tmax
  
  for (mc in 1:nmc) {
    df <- doSim_London(n, tmax, delta, p, q, tstar)
    
    # Compute square root of avg_edges
    sqrt_edges <- sqrt(unlist(df$avg_edges))
    
    # Compute distance matrices
    Dhat_W1 <- getD_W1(df$xhat)
    Dhat_W2 <- getD_W2(df$xhat)
    D_dMV  <- getD(df$xhat)
    
    # Perform MDS
    df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
    df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
    df.mds_dMV <- doMDS(D_dMV, doplot = F)
    
    # Calculate errors using l2
    error_avg_degree_l2 <- find_slope_changepoint_with_plot(sqrt_edges, doplot = F)$error
    error_W1_l2 <- find_slope_changepoint_with_plot(df.mds_W1$mds[,1], doplot = F)$error
    error_W2_l2 <- find_slope_changepoint_with_plot(df.mds_W2$mds[,1], doplot = F)$error
    error_dMV_l2 <- find_slope_changepoint_with_plot(df.mds_dMV$mds[,1], doplot = F)$error
    
    error_avg_degree <- linf_error(sqrt_edges)[1]
    error_W1 <- linf_error(df.mds_W1$mds[,1])[1]
    error_W2 <- linf_error(df.mds_W2$mds[,1])[1]
    error_dMV <- linf_error(df.mds_dMV$mds[,1])[1]
    
    if(abs(error_dMV_l2) - abs(error_dMV) < 0 ){
      print(c(error_dMV_l2 , error_dMV))
      find_slope_changepoint_with_plot(df.mds_dMV$mds[,1], doplot = T)
      ecp_l1 = linf_error(df.mds_dMV$mds[,1])[2]
      abline(v = ecp_l1, col='green' ,lwd = 2, lty = 2)
      
      res <- linf_cp(1:tmax, df.mds_dMV$mds[,1] , ecp_l1)
      t = 1:tmax
      y_fit = res[2] + res[3] * (t -  ecp_l1) + 
        (res[4] - res[3]) * (t -  ecp_l1) * (t >  ecp_l1)
      
      lines(1:tmax, y_fit, col = 'green' , lwd = 2)
      legend("bottomright", legend = c("Fitted Line linf", "Estimated_CP linf"),
             col = c( "green", "green"), lty = c( 1, 2), lwd = 2)
      break
    }
  }
}


example_df.mds_dMV_m_12 = df.mds_dMV$mds[,1]
example_df.mds_dMV_m_40 = df.mds_dMV$mds[,1]



par(mfrow = c(1,2))

a = find_slope_changepoint_with_plot(example_df.mds_dMV_m_12, doplot = T)
m = tmax = length(example_df.mds_dMV_m_12)
linf_error(example_df.mds_dMV_m_12)[2]
ecp_l1 <-  linf_error(example_df.mds_dMV_m_12)[2]
print(ecp_l1)
print(linf_error(example_df.mds_dMV_m_12)[2])
abline(v = ecp_l1, col='green' ,lwd = 2, lty = 2)
m = tmax = length(example_df.mds_dMV_m_12)

res <- linf_cp(1:tmax, example_df.mds_dMV_m_12 , ecp_l1)
t = 1:tmax
y_fit = res[2] + res[3] * (t -  ecp_l1) + 
  (res[4] - res[3]) * (t -  ecp_l1) * (t >  ecp_l1)

lines(1:tmax, y_fit, col = 'green' , lwd = 2)
legend("bottomright", legend = c("Fitted Line linf", "Estimated_CP linf"),
       col = c( "green", "green"), lty = c( 1, 2), lwd = 2)


find_slope_changepoint_with_plot(example_df.mds_dMV_m_40, doplot = T)
ecp_l1 = linf_error(example_df.mds_dMV_m_40)[2]
abline(v = ecp_l1, col='green' ,lwd = 2, lty = 2)
m = tmax = length(example_df.mds_dMV_m_40)

res <- linf_cp(1:tmax, example_df.mds_dMV_m_40 , ecp_l1)
t = 1:tmax
y_fit = res[2] + res[3] * (t -  ecp_l1) + 
  (res[4] - res[3]) * (t -  ecp_l1) * (t >  ecp_l1)

lines(1:tmax, y_fit, col = 'green' , lwd = 2)
legend("bottomright", legend = c("Fitted Line linf", "Estimated_CP linf"),
       col = c( "green", "green"), lty = c( 1, 2), lwd = 2)

