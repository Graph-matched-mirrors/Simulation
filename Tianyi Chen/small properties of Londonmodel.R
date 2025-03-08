Xlist <- NULL
Xt <- matrix(0,n,(tmax+1))



n = 20

p = 0.3
q = 0.3

tstar = 4

tmax = 10 

for (t in 2:(tstar+1)) {
  tmp <- runif(n) < p
  Xt[,t] <- Xt[,t] + Xt[,t-1]
  Xt[tmp,t] <- Xt[tmp,t] + delta
}

for (t in (tstar+2):(tmax+1)) {
  tmp <- runif(n) < q
  Xt[,t] <- Xt[,t] + Xt[,t-1]
  Xt[tmp,t] <- Xt[tmp,t] + delta
}
Xt <- Xt[,-1]
Xt <- Xt+0.1

sum(Xt[,4] - Xt[,2])

sum(sort(Xt[,4]) - sort(Xt[,2]))


df <- doSim_London(n, tmax, delta, p, q, tstar)
df <- df %>% mutate(mean = map(Xt, ~mean(as.vector(unlist(.)))  )) 

mm = unlist(df$mean)
mm = mm - mean(mm)

doMDS(D_W2)


# Compute square root of avg_edges
sqrt_edges <- sqrt(unlist(df$avg_edges))

# Compute distance matrices
D_W1 <- getD_W1(df$Xt)
D_W2 <- getD_W2(df$Xt)
D_dMV  <- getD(df$Xt)


D_dMV
D_W2


D_W1^2


doMDS(D_W1)

D_dMV[1,2]
D_W1[1,2]

Xhati = df$xhat[[1]]
Xhatj = df$xhat[[5]]

diff_hat = as.vector(Xhatj - Xhati)

sort_diff_hat =  sort(Xhatj) - sort(Xhati)


hist(diff_hat)

hist(sort_diff_hat)

#sqrt(mean((Xhati - Xhatj)^2))

sqrt(mean(diff_hat^2))
sqrt(mean( sort_diff_hat^2 ))




plot(1:n, diff)

(mean(  sort_diff_hat  ))
(mean(  diff_hat  ))


Xhati = df$Xt[[1]]
Xhatj = df$Xt[[5]]

diff = as.vector(Xhatj - Xhati)

sort_diff =  sort(Xhatj) - sort(Xhati)


hist(diff)

hist(sort_diff)

#sqrt(mean((Xhati - Xhatj)^2))

sqrt(mean(diff^2))
sqrt(mean( sort_diff^2 ))



plot(1:n, diff)

(mean(  sort_diff  ))
(mean(  diff  ))


