
error_avg_edges <- error_dMV <- NULL


n = 200
tmax <- m <- 20
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
tstar <- tmax / 2 

nmc = 1000

for(mc in 1:nmc){
  
  df <- doSim_London(n,tmax,delta,p,q,tstar)
  #Dhat_W1 = getD_W1(df$xhat)
  #Dhat_W2 = getD_W2(df$xhat)
  D_dMV = getD(df$xhat)
  
  #df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
  #df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
  df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)
  
  error_avg_edges[mc] = linf_error(sqrt(unlist(df$avg_edges)))
  #linf_error(df.mds_W1$mds[,1])
  error_dMV[mc] = linf_error(df.mds_dMV$mds[,1])
}

hist(error_avg_edges)
hist(error_dMV, breaks = 50)


mean(error_dMV^2)


n = 200
tmax <- m <- 40
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
tstar <- tmax / 2 

nmc = 500

for(mc in 1:nmc){
  
  df <- doSim_London(n,tmax,delta,p,q,tstar)
  #Dhat_W1 = getD_W1(df$xhat)
  #Dhat_W2 = getD_W2(df$xhat)
  D_dMV = getD(df$xhat)
  
  #df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
  #df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
  df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)
  
  error_avg_edges[mc] = linf_error(sqrt(unlist(df$avg_edges)))
  #linf_error(df.mds_W1$mds[,1])
  error_dMV[mc] = linf_error(df.mds_dMV$mds[,1])
}

hist(
  error_avg_edges,
  main = paste('n =', n, 
               'nmc =', nmc, 
               'm =', m, 
               'p =', p, 
               'q =', q),
  xlab = 'Error using Avg degree',
  breaks = 50         # Adjust number of bins
)

hist(
  error_dMV,
  main = paste('n =', n, 
               'nmc =', nmc, 
               'm =', m, 
               'p =', p, 
               'q =', q),
  xlab = 'Error using dMV ',
  breaks = 70         # Adjust number of bins
)
mean(abs(error_dMV)^2)


df <- doSim_London(n,tmax,delta,p,q,tstar)
#Dhat_W1 = getD_W1(df$xhat)
#Dhat_W2 = getD_W2(df$xhat)
D_dMV = getD(df$xhat)

#df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
#df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)

#error_avg_edges[mc] = linf_error(sqrt(unlist(df$avg_edges)))
#linf_error(df.mds_W1$mds[,1])
linf_error(df.mds_dMV$mds[,1])

