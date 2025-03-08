## London model:

n=100
tmax=400
delta= (1 - 0.1)/ tmax
p=0.4
q=0.6
tstar=200
startpoint=0.1

df <- doSim_London(n,tmax,delta, p ,q, tstar,del = c(1)  )

del = c(1)

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) )
}

D2 <- getD(df$Xt)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)

hist(D2/D2_shuffle_1, main = paste('Quantity for London Model','n =', n, 
                                                                          'm =', tmax, 
                                                                          'p =', p, 
                                                                          'q =', q)  )
##Atlanta Model

n = 500
p = 0.4
q = 0.1
m = 20
num_state = 400
c=(.9-0.1)
delta = c/(num_state-1)
tstar = m/2
del = c(0.01, 0.1, 0.5, 1)

xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

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
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
  mutate(g = map(Xt, ~rdpg.sample(.))) %>%
  mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))


for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~shuffle_X(., perc)) )
}



D2hat <- getD(df$xhat)
D2hat_shuffle_1 <- getD(df$shuffle_Xhat_1)

D2hat_shuffle_0.01 <- getD(df$shuffle_Xhat_0.01)


par(mfrow=c(1,2))
hist(D2hat/D2hat_shuffle_1,main = paste('Quantity for Atlanta Model','n =', n, 
                                  'm =', tmax, 
                                  'p =', p, 
                                  'q =', q),xlim = c(0, 1) )


hist(D2hat/D2hat_shuffle_0.01,main = paste('Quantity for Atlanta Model','n =', n, 
                                  'm =', tmax, 
                                  'p =', p, 
                                  'q =', q),xlim = c(0, 1) )


D2 <- getD(df$Xt)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)
D2_shuffle_0.01 <- getD(df$shuffle_Xt_0.01)

par(mfrow=c(1,2))
hist(D2/D2_shuffle_1,main = paste('Quantity for Atlanta Model','n =', n, 
                                  'm =', tmax, 
                                  'p =', p, 
                                  'q =', q),xlim = c(0, 1) )


hist( D2/D2_shuffle_0.01, main = paste('Quantity for Atlanta Model','n =', n, 
                                    'm =', tmax, 
                                    'p =', p, 
                                    'q =', q),xlim = c(0, 1) )

summary(na.omit(as.vector(D2/D2_shuffle_0.01)))
