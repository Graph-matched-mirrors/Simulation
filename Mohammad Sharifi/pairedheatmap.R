jump_Lbd=function(cp,p, delta){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=delta
  }
  return(np)
}


jump_Rbd=function(cp,p, delta){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=1-delta
  }
}

jump_middle=function(cp,p, delta){
  u=runif(1)
  if(u<p){
    
    np=cp+delta
    
  }
  
  if(p<u & u<2*p){
    np=cp-delta
  }
  
  if(u>2*p){
    np=cp
  }
  return(np)
}



update_function=function(current_position,p , delta){
  if( abs(current_position-0)<10^(-2) ){
    next_position=jump_Lbd(current_position,p,delta)
  }
  if(current_position> (0+10^(-2)) & current_position <(1-10^(-2)) ){
    next_position=jump_middle(current_position,p,delta)
  }
  if(abs(current_position-1)<10^(-2)){
    next_position=jump_Rbd(current_position,p,delta)
  }
  return(next_position)
}


shuffle_perc_graph <- function(A, del, n){
  G=as.matrix(A)
  print
  dn = floor(n*del)
  if(dn == 0)
    return(A)
  else{
    permu_vec=sample(1:dn)
    random_perm=diag(dn)[permu_vec,]
    a=bdiag(diag(n-dn), random_perm)
    random_perm = as.matrix(a)
    permu_G=as.matrix(random_perm%*%G%*%t(random_perm))
    G_graph=graph_from_adjacency_matrix(permu_G,mode ="undirected")
    return(G_graph)
  }
}

paired_error_in_shuffling_once <- function(n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = c(0.1,0.2)){
  
  xt=matrix(0,nrow = n, ncol = m+1)
  initial_state_all_nodes=sample(seq(0,1,by=delta),n,replace = TRUE)
  xt[,1]=initial_state_all_nodes
  for (i in 1:n) {
    for (j in 2:(tstar)) {
      xt[i,j]=update_function(xt[i,j-1],p,delta)
    }
    for (j in ((tstar+1):(m+1)) ) {
      xt[i,j]=update_function(xt[i,j-1],q,delta)
    }
  }
  
  
  
  
  df <- tibble(time=1:m) %>%
    mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.)))%>%
    mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

  
  for(perc in del){
    df <- df %>%
      mutate(!!paste0("shuffle_", perc) := map(g,~shuffle_perc_graph(., perc, n))) %>%
      mutate(!!paste0("xhat_", perc) := map(!!sym(paste0("shuffle_", perc)), function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  }
  
  ### make graph

  D2=getD(df$xhat)
  df.mds <- doMDS(D2,doplot = FALSE)
  df.iso <- doIso(df.mds, mdsd=10, AtlantaFlag=T)
  errors <- NULL
  errors[1] <- linf_error(df.iso$iso, m)
  i <- 2
  for(perc in del){
    D2_shuffle=getD(df[[paste0("xhat_", perc)]])
    df.mds_shuffle <- doMDS(D2_shuffle,doplot = FALSE)
    df.iso_shuffle <- doIso(df.mds_shuffle, mdsd=10)
    errors[i] <- linf_error(df.iso_shuffle$iso, m)
    i <- i + 1
  }
  print(paste(q,Sys.time()))
  errors
}
paired_error_in_shuffling <- function(nmc = 50, n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = 0.1){
  mc_errors <- sapply(1:nmc, function(i) paired_error_in_shuffling_once(n, p, q, m, delta, tstar, del))
  row_mse <- apply(mc_errors, 1, function(row) mean(abs(row)^2)) # mean(abs(errors[row,])^2)
  row_sds <- apply(mc_errors, 1, function(row) sd(abs(row)^2) / sqrt(nmc)) # sd(abs(errors[row,])^2) / sqrt(nmc)
  cbind(row_mse, row_sds)
}
m <- 50
tstar <- 25
delta <- 0.1
p <- 0.4
q <- seq(0,0.5, by = 0.05)
d <- seq(0.05, 0.30, by = 0.05)
n <- 300
nmc <- 100
final_errors <- NULL
delta <- 0.1
for(i in 1:length(q)){
  temp_errors <- paired_error_in_shuffling(nmc = nmc, n = n, p = p, q = q[i], m = m, delta = delta, tstar = tstar, del = d)
  final_errors[[paste0("row_mse_", q[i])]] <- temp_errors[,1]
  final_errors[[paste0("row_sds_", q[i])]] <- temp_errors[,2]
  print(paste(q[i]," is done! "))
  print(final_errors[[paste0("row_mse_", q[i])]])
  print(final_errors[[paste0("row_sds_", q[i])]])
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
write.csv(final_errors, "~/final_errors.cvs", row.names = FALSE)
