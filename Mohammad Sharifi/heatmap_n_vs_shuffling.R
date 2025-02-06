m <- 50
tstar <- 25
delta <- 0.1
p <- 0.4
q <- 0.25
d <- seq(0.05, 1, by = 0.05)
set_of_n <- c(100,300,800)
nmc <- 100
final_errors_diff_n <- NULL
delta <- 0.1
for(i in 1:length(set_of_n)){
  temp_errors <- paired_error_in_shuffling(nmc = nmc, n = set_of_n[i], p = p, q = q, m = m, delta = delta, tstar = tstar, del = d)
  final_errors_diff_n[[paste0("row_mse_", set_of_n[i])]] <- temp_errors[,1]
  final_errors_diff_n[[paste0("row_sds_", set_of_n[i])]] <- temp_errors[,2]
  print(paste(set_of_n[i]," is done! "))
  print(final_errors_diff_n[[paste0("row_mse_", set_of_n[i])]])
  print(final_errors_diff_n[[paste0("row_sds_", set_of_n[i])]])
}

mse_keys <- grep("^row_mse_", names(final_errors_diff_n), value = TRUE) # Find keys starting with "row_mse_"
mse_values <- final_errors_diff_n[mse_keys] # Extract the corresponding elements
mse_matrix <- do.call(cbind, mse_values) # Combine into a matrix
colnames(mse_matrix) <- set_of_n
rownames(mse_matrix) <- c(0,d)


sd_keys <- grep("^row_sds_", names(final_errors_diff_n), value = TRUE) # Find keys starting with "row_mse_"
sd_values <- final_errors_diff_n[sd_keys] # Extract the corresponding elements
sd_matrix <- do.call(cbind, sd_values) # Combine into a matrix
colnames(sd_matrix) <- set_of_n
rownames(sd_matrix) <- c(0,d)

heatmap(mse_matrix, Rowv = NA, Colv = NA,
        main = paste("q =",q, ", p =", p, ", nmc =", nmc),
        xlab = "n", ylab = "shuffling ratio", scale = "none") 
write.csv(final_errors_diff_n, "~/final_errors_diff_n.cvs", row.names = FALSE)
