

## this document is meant to show that entrywise square is will give lienar for N big, but the ratio shuffling make it not linear
## while using non square+ isomap+ mdsd big is still linear. 
p = 0.4
q = 0.2
m = 40
tstar = 20

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)


True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.

True_shuffle_dmv_0.4_square= 0.6*True_dmv_square + 0.4*True_shuffle_dmv_square
doMDS(True_shuffle_dmv_0.4_square, doplot=TRUE)

A_shuffle_dmv_0.4 = doMDS(sqrt(True_shuffle_dmv_0.4_square), doplot=TRUE)

doIso(A_shuffle_dmv_0.4$mds, mdsd = 10 , doplot = T)

ev_in_mds = apply(A_shuffle_dmv_0.4$mds, 2, sd)
elbow = getElbows( ev_in_mds )
elbow
doIso(A_shuffle_dmv_0.4$mds, mdsd = elbow[1], doplot = T)
doIso(A_shuffle_dmv_0.4$mds, mdsd = elbow[2], doplot = T)
doIso(A_shuffle_dmv_0.4$mds, mdsd = elbow[3], doplot = T)


doMDS(True_shuffle_dmv_square, doplot = T)
doMDS(sqrt(True_shuffle_dmv_square), doplot = T)



True_dmv=sqrt(True_dmv_square)

Atlanta_mds = doMDS(True_dmv, doplot=TRUE)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:(m-1), ev_in_mds, main = 'scree plot') ## elbow actually says you need at least d = 3
elbow = getElbows( ev_in_mds )
elbow
## from below you see how it gets more linear 
doIso(Atlanta_mds$mds, mdsd = elbow[1], doplot = T)
doIso(Atlanta_mds$mds, mdsd = elbow[2], doplot = T)
doIso(Atlanta_mds$mds, mdsd = elbow[3], doplot = T)
doIso(Atlanta_mds$mds, mdsd = 10 , doplot = T)


doMDS(True_dmv_square, doplot=TRUE)
doMDS(sqrt(True_dmv_square), doplot=TRUE)


## finite sample behavior
p = 0.4
q = 0.2
m = 40
tstar = 20
c=(.9-0.1)
num_state = 100
delta = c/(num_state-1)
True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta)
True_dmv=sqrt(True_dmv_square)
n = 300
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
df <- NULL
df <- tibble(time=1:m) %>%
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))%>%
  mutate(g = map(Xt, ~rdpg.sample(.))) %>%
  mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

#sum((getD(df$Xt)-getD(df$xhat))^2)
twotoinf = NULL
dMVhat = NULL
diff_perc = NULL
for (kkk in 1:m) {
  k = kkk
  twotoinf[k] = min(max((df$Xt[[k]]-df$xhat[[k]])^2),max((df$Xt[[k]]+df$xhat[[k]])^2))
  dMVhat[k] = min(sum((df$Xt[[k]]-df$xhat[[k]])^2)/n,sum((df$Xt[[k]]+df$xhat[[k]])^2)/n)
  diff_perc [k] = min(max(abs(df$Xt[[k]]+df$xhat[[k]])/df$Xt[[k]]), max(abs(df$Xt[[k]]-df$xhat[[k]])/df$Xt[[k]]))
}

summary(twotoinf)
summary(dMVhat)
summary(diff_perc)

summary(as.vector(abs((getD(df$xhat)-True_dmv))/True_dmv))


#summary(as.vector(True_dmv_square))

#mean(unlist(df$avg_edges))
doIso(doMDS(getD(df$Xt), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$xhat), doplot = F)$mds, mdsd = 10 , doplot = T)


del = c(1,0.1,0.2,0.25,0.3,0.4,0.5,0.6)
#del = c(1)

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~optimized_shuffled_X(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) )
}

doIso(doMDS(getD(df$shuffle_Xhat_0.1), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$shuffle_Xhat_0.2), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$shuffle_Xhat_0.25), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$shuffle_Xhat_0.3), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$shuffle_Xhat_0.4), doplot = F)$mds, mdsd = 10 , doplot = T)
doIso(doMDS(getD(df$shuffle_Xhat_0.6), doplot = F)$mds, mdsd = 10 , doplot = T)








