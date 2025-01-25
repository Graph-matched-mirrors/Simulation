pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)

pacman::p_load("doParallel")
registerDoParallel(detectCores()-1)

## mds_nochange and df_true_nochange generation
set.seed(2)
n = 200
tmax <- m <- 20
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
tstar <- tmax / 2 


df <- doSim_London(n,tmax,delta,p,q,tstar)

#plot(1:tmax, sqrt(unlist(df$avg_edges)))

#True_W1_D = sqrt(true_W1_square_London(tmax,tstar,p,q))
D_True_dMV = true_London_dMV(tmax, tstar , p ,q  )

#summary(as.vector((Dhat_W1 - True_W1_D) / True_W1_D))
Dhat_W1 = getD_W1(df$xhat)
Dhat_W2 = getD_W2(df$xhat)
D_dMV = getD(df$xhat)
#summary(as.vector((Dhat_W2  - D_dMV ) / D_dMV ))

df.mds_W1 <- doMDS(Dhat_W1, doplot = T)
#df.mds_True_W1 <- doMDS(True_W1_D, doplot = T)

df.mds_W2 <- doMDS(Dhat_W2, doplot = T)
df.mds_dMV <- doMDS(D_dMV, doplot = T)
doMDS(D_True_dMV, doplot = T)
## what makes dMV so noisy? 

linf_error(df.mds_W1$mds[,1])*m
linf_error(df.mds_dMV$mds[,1])*m
linf_error(sqrt(unlist(df$avg_edges)))*m


n = 200
tmax <- m <- 40
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
tstar <- tmax / 2 

df <- doSim_London(n,tmax,delta,p,q,tstar)
Dhat_W1 = getD_W1(df$xhat)
#Dhat_W2 = getD_W2(df$xhat)
D_dMV = getD(df$xhat)


df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)

linf_error(sqrt(unlist(df$avg_edges)))*m
linf_error(df.mds_W1$mds[,1])*m
linf_error(df.mds_dMV$mds[,1])*m

plot(1:tmax, df.mds_dMV$mds[,1])


set.seed(2)
n = 200
p <- 0.4
q <- 0.3
nmc = 1000

# Define the values of mm
mm <- c(12, 20, 30, 40)

# Initialize a list to store the results
results_list <- list()

# Define error types
error_types <- c("error_avg_degree", "error_W1", "error_W2", "error_dMV")

# Loop over each value of m in mm
for (i in 1:length(mm)) {
  tmax <- m <- mm[i]
  tstar <- tmax / 2
  delta <- (1-0.1)/tmax
  # Dynamically create and assign a matrix
  matrix_name <- paste0("error_matrix_m_", m)
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
    
    # Calculate errors
    error_avg_degree <- linf_error(sqrt_edges)
    error_W1 <- linf_error(df.mds_W1$mds[,1])
    error_W2 <- linf_error(df.mds_W2$mds[,1])
    error_dMV <- linf_error(df.mds_dMV$mds[,1])

    
    #print(error_dMV)
    
    # Retrieve and update the matrix dynamically
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
  mse_mean <- apply(abs(final_matrix)^2, 1, mean)
  mse_sd <- apply(abs(final_matrix)^2, 1, sd)
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

error_matrix_m_12
hist(error_matrix_m_12[4, ], breaks = 30)
hist(error_matrix_m_20[4, ], breaks = 30)
hist(error_matrix_m_30[4, ], breaks = 30)
hist(error_matrix_m_40[4, ], breaks = 30)


error_matrix_m_40

# Combine all results into a single dataframe
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

# Display the final summary table
print(summary_df)



################# IGNORE
ggplot(summary_df, aes(x = as.factor(m), y = Mean_MSE, color = error_type)) +
  geom_point(size = 3) +  # Points for the mean values
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +  # Error bars
  labs(
    title = "Error Analysis across Different Values of m",
    x = "m",
    y = "Mean Squared Error",
    color = "Error Type"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "top",  # Position the legend at the top
    text = element_text(size = 14)  # Increase text size for readability
  )






