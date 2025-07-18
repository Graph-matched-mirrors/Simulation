set.seed(10)
m <- 40
tstar <- 20
Num_states <- m 
p <- 0.4
q <- c(0.1, 0.2, 0.35, 0.4, 0.5)
d <- seq(0.05, 1, by = 0.05)
n <- 300
nmc <- 100
final_errors <- NULL
for(i in 1:length(q)){
  temp_errors <- paired_error_in_shuffling(nmc = nmc, n = n, p = p, q = q[i], m = m, Num_states = Num_states, tstar = tstar, del = d)
  final_errors[[paste0("row_mse_", q[i], "_d_l2")]] <- temp_errors[seq(1, nrow(temp_errors)/2), 1]
  final_errors[[paste0("row_sds_", q[i], "_d_l2")]] <- temp_errors[seq(1, nrow(temp_errors)/2), 2]
  final_errors[[paste0("row_mse_", q[i], "_d_li")]] <- temp_errors[seq(1 * nrow(temp_errors)/2 + 1, 2 * nrow(temp_errors)/2), 1]
  final_errors[[paste0("row_sds_", q[i], "_d_li")]] <- temp_errors[seq(1 * nrow(temp_errors)/2 + 1, 2 * nrow(temp_errors)/2), 2]
  print(paste(q[i]," is done! "))
  print(final_errors[[paste0("row_mse_", q[i], "_d_l2")]])
  print(final_errors[[paste0("row_sds_", q[i],"_d_l2")]])
  print(final_errors[[paste0("row_mse_", q[i], "_d_li")]])
  print(final_errors[[paste0("row_sds_", q[i],"_d_li")]])
}

mse_keys <- grep("^row_mse_", names(final_errors), value = TRUE) # Find keys starting with "row_mse_"
mse_values <- final_errors[mse_keys] # Extract the corresponding elements
mse_matrix <- do.call(cbind, mse_values) # Combine into a matrix
colnames(mse_matrix) <- q 
rownames(mse_matrix) <- c(0,d)


sd_keys <- grep("^row_sds_", names(final_errors), value = TRUE) # Find keys starting with "row_mse_"
sd_values <- final_errors[sd_keys] # Extract the corresponding elements
sd_matrix <- do.call(cbind, sd_values) # Combine into a matrix
colnames(sd_matrix) <- q 
rownames(sd_matrix) <- c(0,d)

heatmap(mse_matrix, Rowv = NA, Colv = NA,
        main = paste("n =",n, ", p =", p, ", nmc =", nmc),
        xlab = "q", ylab = "shuffling ratio", scale = "none") 
write.csv(final_errors, "~/final_errors_London.cvs", row.names = FALSE)
