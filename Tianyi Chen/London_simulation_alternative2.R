# Define the values of mm
set.seed(2)
n = 200
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
nmc = 500

mm <- c(10, 20, 30, 40)

# Initialize a list to store the results
results_list <- list()

# Define error types
error_types <- c("error_avg_degree", "error_W1", "error_W2", "error_dMV")

# Loop over each value of m in mm
for (i in 1:length(mm)) {
  m <- mm[i]
  
  # Temporary matrix to store results for current m
  error_matrix <- matrix(nrow = 4, ncol = nmc)
  
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
    df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)
    
    # Calculate errors
    error_avg_degree <- linf_error(sqrt_edges)
    error_W1 <- linf_error(df.mds_W1$mds[,1])
    error_W2 <- linf_error(df.mds_W2$mds[,1])
    error_dMV <- linf_error(df.mds_dMV$mds[,1])
    
    print(c(m,error_dMV))
    # Store errors in the temporary matrix
    error_matrix[, mc] <- c(error_avg_degree, error_W1, error_W2, error_dMV)
  }
  
  # Compute mean squared error and confidence interval
  mse_mean <- apply(abs(error_matrix)^2, 1, mean)
  mse_sd <- apply(abs(error_matrix)^2, 1, sd)
  conf_interval <- 1.96 * mse_sd / sqrt(nmc)
  
  # Create a dataframe for the current m value
  temp_df <- data.frame(
    m = rep(m, 4),
    error_type = error_types,
    Lower_CI = mse_mean - conf_interval,
    Mean_MSE = mse_mean,
    Upper_CI = mse_mean + conf_interval
  )
  
  # Store the dataframe in the results list
  results_list[[i]] <- temp_df
}

# Combine all results into a single dataframe
summary_df <- do.call(rbind, results_list)


m
data.frame(
  m = rep(m, 4),
  error_type = error_types,
  Lower_CI = mse_mean - conf_interval,
  Mean_MSE = mse_mean,
  Upper_CI = mse_mean + conf_interval
)

# Display the final summary table
print(summary_df)

summary_df <- do.call(rbind, results_list)

plottt1 <- ggplot(summary_df, aes(x = m, y = Mean_MSE, color = error_type)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm', title =  paste('n =', n, 
                                          'nmc =', nmc, 
                                          'p =', p, 
                                          'q =', q)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

# Print the plot
print(plottt1)
