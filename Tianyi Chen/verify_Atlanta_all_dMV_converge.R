optimized_shuffled_X <- function(X, del) {
  n <- nrow(X)
  dn <- floor(del * n)
  
  perm_indices <- sample(n, dn)  # Randomly pick `dn` row indices
  permuted_indices <- sample(perm_indices)  # Shuffle only those indices
  
  X[perm_indices, ] <- X[permuted_indices, ]
  
  return(X)
}

set.seed(1)
p = 0.05
q = 0.45
tmax = m = 30
tstar = 15

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)


True_D=sqrt(true_Atlanta_dmv(p,q,num_state,m,tstar,delta)) ## This is the analytically dMV result.
MDS_True_D =doMDS(True_dmv, doplot = T)

True_shuffled_D_square =true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffled_D = sqrt(True_shuffled_D_square)


del = c(0.1,0.19,0.2,0.36,0.3,0.4,0.5,1)

for (alpha in del) {
  nm <- paste0("True_shuffle_dmv_", alpha)
  assign(nm,
         sqrt((1 - alpha) * True_D^2 + alpha * True_shuffled_D^2),
         envir = .GlobalEnv)
}


MDS_True_shuffle_dmv_0.5 = doMDS(sqrt(True_shuffle_dmv_0.5), doplot = T)



## verify the scaling constant
sqrt(0.5+0.5*True_shuffle_dmv_square[1,2]/(2*eigen(B)$values[1]))
MDS_True_shuffle_dmv_0.5$mds[,1]/MDS_True_D$mds[,1]


n = 1000
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

##compare theoretical variance with real variance
the_var = 0.8^2/12*( (num_state+1)/(num_state-1) )
print(paste('theretical variacne = ', the_var , 'sample variance = ',var(initila_state_all_nodes) )  )

mean(initila_state_all_nodes)
the_var = c^2/12*( (num_state+1)/(num_state-1) )
print(paste('theretical variacne = ', the_var , 'sample variance = ',var(initila_state_all_nodes) )  )
sqrt(2*the_var) #is the entry of true 100% shuffling dMV
True_shuffle_dmv[1,2]

xt[,1]=initila_state_all_nodes
for (i in 1:n) {
  for (j in 2:(tstar)) {
    xt[i,j]=update_function(xt[i,j-1],p)
  }
  for (j in ((tstar+1):(m+1)) ) {
    xt[i,j]=update_function(xt[i,j-1],q)
  }
}

df <- tibble(time=1:m) %>%
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))  

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle1_Xt_", perc) := map( Xt ,~optimized_shuffled_X(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~ shuffle_X_optimized(., perc)) )
  #%>%
  #mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) )
}

for (alpha in del) {
  nm <- paste0("D_shuffle_", alpha)
  col_name <- paste0("shuffle_Xt_", alpha)
  assign(nm,
         getD(df[[col_name]]),
         envir = .GlobalEnv)
}

for (alpha in del) {
  nm <- paste0("D1_shuffle_", alpha)
  col_name <- paste0("shuffle1_Xt_", alpha)
  assign(nm,
         getD(df[[col_name]]),
         envir = .GlobalEnv)
}

summary(na.omit(abs(as.vector((D_shuffle_0.1 - True_shuffle_dmv_0.1)/True_shuffle_dmv_0.1)))) #the max should decrease as n increase with 


D <- getD(df$Xt)

summary(na.omit(abs(as.vector((D- True_D)/True_D))))
summary(na.omit(abs(as.vector((D_shuffle_0.1 - True_shuffle_dmv_0.1)/True_shuffle_dmv_0.1))))
summary(na.omit(abs(as.vector((D_shuffle_0.19 - True_shuffle_dmv_0.19)/True_shuffle_dmv_0.19))))
summary(na.omit(abs(as.vector((D_shuffle_0.2 - True_shuffle_dmv_0.2)/True_shuffle_dmv_0.2))))
summary(na.omit(abs(as.vector((D_shuffle_1 - True_shuffle_dmv_1)/True_shuffle_dmv_1))))


## This is for random shuffle all not fixed part. 
summary(na.omit(abs(as.vector((D1_shuffle_0.1 - True_shuffle_dmv_0.19)/True_shuffle_dmv_0.19))))
summary(na.omit(abs(as.vector((D1_shuffle_0.2 - True_shuffle_dmv_0.36)/True_shuffle_dmv_0.36))))



doMDS(D2_shuffle_0.1,doplot = T)




## relation between YP's scree plot which is the sd of MDS result col wise and the eigenvalue of B matrix
P=diag(rep(1,tmax))-rep(1,tmax)%*%t(rep(1,tmax))/tmax
B=(-1/(2))*P%*%True_D^2%*%P

eigen(B)$values[1]
sum(MDS_True_D$mds[,1]^2)

apply(MDS_True_D$mds, 2, sd)[1]
sqrt(sum(MDS_True_D$mds[,1]^2)/(m-1))
