set.seed(2)
n = 200
p <- 0.4
q <- 0.3
nmc = 500

# Define the values of mm
mm <- c(12, 20, 30, 40)

# Initialize a list to store the results
list_name_l2 <- paste0("results_list_l2", "_n_",n)
assign(list_name_l2, list())

list_name <- paste0("results_list", "_n_",n)
assign(list_name, list())
# Define error types
error_types_l2 <- c("error_avg_degree_l2", "error_W1_l2", "error_W2_l2", "error_dMV_l2")
error_types <- c("error_avg_degree", "error_W1", "error_W2", "error_dMV")

# Loop over each value of m in mm
for (i in 1:length(mm)) {
  tmax <- m <- mm[i]
  tstar <- tmax / 2
  delta <- (1-0.1)/tmax
  # Dynamically create and assign a matrix
  matrix_name_l2 <- paste0("error_matrix_l2_m_", m, "_n_",n)
  matrix_name <- paste0("error_matrix_m_", m, "_n_",n)
  
  assign(matrix_name_l2, matrix(nrow = 4, ncol = nmc))
  assign(matrix_name, matrix(nrow = 4, ncol = nmc))
  
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
    
    # Calculate errors using l2
    error_avg_degree_l2 <- find_slope_changepoint_with_plot(sqrt_edges, doplot = F)$error
    error_W1_l2 <- find_slope_changepoint_with_plot(df.mds_W1$mds[,1], doplot = F)$error
    error_W2_l2 <- find_slope_changepoint_with_plot(df.mds_W2$mds[,1], doplot = F)$error
    error_dMV_l2 <- find_slope_changepoint_with_plot(df.mds_dMV$mds[,1], doplot = F)$error
    
    error_avg_degree <- linf_error(sqrt_edges)[1]
    error_W1 <- linf_error(df.mds_W1$mds[,1])[1]
    error_W2 <- linf_error(df.mds_W2$mds[,1])[1]
    error_dMV <- linf_error(df.mds_dMV$mds[,1])[1]
    #print(error_dMV)
    
    # Retrieve and update the matrix dynamically
    temp_matrix_l2 <- get(matrix_name_l2)
    temp_matrix_l2[, mc] <- c(error_avg_degree_l2, error_W1_l2, error_W2_l2, error_dMV_l2)
    assign(matrix_name_l2, temp_matrix_l2)
    
    temp_matrix <- get(matrix_name)
    temp_matrix[, mc] <- c(error_avg_degree, error_W1, error_W2, error_dMV)
    assign(matrix_name, temp_matrix)
    
    
    if (mc %% 100 == 0) {
      print(c(i, mc))
    }
  }
  
  # Retrieve the updated matrix
  final_matrix <- get(matrix_name)
  
  # Compute mean squared error and confidence interval
  mse_mean <- apply(final_matrix^2, 1, mean)
  mse_sd <- apply(final_matrix^2, 1, sd)
  conf_interval <- 1.96 * mse_sd / sqrt(nmc)
  
  # Create a dataframe for the current m value
  temp_df <- data.frame(
    m = rep(m, 4),
    error_type = error_types,
    Lower_CI = mse_mean - conf_interval,
    Mean_MSE = mse_mean,
    Upper_CI = mse_mean + conf_interval
  )
  
  current_list <- get(list_name)
  current_list[[i]] <- temp_df
  
  assign(list_name, current_list)
  
  final_matrix_l2 <- get(matrix_name_l2)
  
  # Compute mean squared error and confidence interval
  mse_mean_l2 <- apply(final_matrix_l2^2, 1, mean)
  mse_sd_l2 <- apply(final_matrix_l2^2, 1, sd)
  conf_interval_l2 <- 1.96 * mse_sd_l2 / sqrt(nmc)
  
  # Create a dataframe for the current m value
  temp_df_l2 <- data.frame(
    m = rep(m, 4),
    error_type = error_types_l2,
    Lower_CI = mse_mean_l2 - conf_interval_l2,
    Mean_MSE = mse_mean_l2,
    Upper_CI = mse_mean_l2 + conf_interval_l2
  )
  current_list_l2 <- get(list_name_l2)
  current_list_l2[[i]] <- temp_df_l2
  assign(list_name_l2, current_list_l2)
}



summary_df_l2 <- do.call(rbind, c(get(list_name_l2)))

summary_df_both <- do.call(rbind, c(get(list_name_l2), get(list_name)))
summary_df_both$localizer <- ifelse(grepl("_l2", summary_df_both$error_type), "L2", "Linf")
summary_df_both$metric <- gsub("error_", "", gsub("_l2", "", summary_df_both$error_type))


summary_df_both

# Update your ggplot call to include linetype mapping
plottt1 <- ggplot(summary_df_both, aes(x = m, y = Mean_MSE, color = metric, linetype = localizer)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm', title = paste('L2-localizer','n =', n, 
                                         'nmc =', nmc, 
                                         'p =', p, 
                                         'q =', q)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

# Print the plot
print(plottt1)


plottt1 <- ggplot(summary_df_both , aes(x = m, y = Mean_MSE, color = error_type)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm', title =  paste('L2-localizer','n =', n, 
                                          'nmc =', nmc, 
                                          'p =', p, 
                                          'q =', q)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

# Print the plot
print(plottt1)


plottt2 <- ggplot(summary_df_l2 , aes(x = m, y = Mean_MSE, color = error_type)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm', title =  paste('L2-localizer','n =', n, 
                                          'nmc =', nmc, 
                                          'p =', p, 
                                          'q =', q)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

# Print the plot
print(plottt2)

summary_df_l2 = summary_df_l2[-c(3,7,11,15),]

plottt2 <- ggplot(summary_df_l2 , aes(x = m, y = Mean_MSE, color = error_type)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm', title =  paste('L2-localizer','n =', n, 
                                          'nmc =', nmc, 
                                          'p =', p, 
                                          'q =', q)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

# Print the plot
print(plottt2)

#c(error_avg_degree_l2, error_W1_l2, error_W2_l2, error_dMV_l2)
hist(temp_matrix_l2[2,]^2,breaks = 50)
mean(temp_matrix_l2[2,]^2)

table(temp_matrix_l2[2,]^2)

hist(temp_matrix_l2[4,]^2,breaks = 50)
mean(temp_matrix_l2[4,]^2)
table(temp_matrix_l2[4,]^2)
