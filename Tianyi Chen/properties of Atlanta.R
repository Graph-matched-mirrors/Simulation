pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph,iGraphMatch,ggrepel)
registerDoParallel(detectCores()-1)

# Here we summarize properties we obeserve for Atlanta model:


## when n = infty: xhat = x and dmvhat = dmv 

p = 0.4
q = 0.2
m = 50
tstar = m/2

c=(.9-0.1)
num_state = 10
delta = c/(num_state-1)

del = c(1, 0.4) ## shuffling percentage 


## 1, iso+mds with mdsd large can get you picewise linear(no entrywise square the dMV!)
True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

Atlanta_mds = doMDS(True_dmv, doplot=TRUE)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:(m-1), ev_in_mds, main = 'scree plot') ## elbow actually says you need at least d = 3
elbow = getElbows( ev_in_mds )

## from below you see how it gets more linear 
doIso(Atlanta_mds$mds, mdsd = elbow[1], doplot = T)
doIso(Atlanta_mds$mds, mdsd = elbow[2], doplot = T)
doIso(Atlanta_mds$mds, mdsd = elbow[3], doplot = T)


## 2,  large N with entrywise square make the 1d of MDS linear:
num_state = 100
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

Atlanta_mds = doMDS(True_dmv, doplot=TRUE)
Atlanta_mds = doMDS(True_dmv^2, doplot=TRUE)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:(m-1), ev_in_mds, main = 'scree plot') ## elbow actually says you only need d =1 
elbow = getElbows( ev_in_mds )


## 2.5, interesting finding is when N is small, square is not linear enough, but 4th is...
num_state = 10
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

Atlanta_mds = doMDS(True_dmv, doplot=TRUE)
Atlanta_mds = doMDS(True_dmv^2, doplot=TRUE)
Atlanta_mds = doMDS(True_dmv^4, doplot=TRUE)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:length(ev_in_mds), ev_in_mds, main = 'scree plot') ## elbow actually says you only need d =1 
elbow = getElbows( ev_in_mds )


## 100% shuffling behavior

True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )


##c percent shuffling behavior: 

True_shuffle_dmv_0.4_square= 0.6*True_dmv_square + 0.4*True_shuffle_dmv_square



## ZC quantity behavior, change num_state N and see behavior. 
num_state = 10
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffle_dmv = sqrt(True_shuffle_dmv_square)

hist(True_dmv / True_shuffle_dmv , ,xlim = c(0, 1))




## finite sample behavior:
p = 0.4
q = 0.2
m = 30
tstar = m/2

c=(.9-0.1)
num_state = 100
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffle_dmv_0.4_square= 0.6*True_dmv_square + 0.4*True_shuffle_dmv_square

True_dmv=sqrt(True_dmv_square)
True_shuffle_dmv = sqrt(True_shuffle_dmv_square)
True_shuffle_dmv_0.4=sqrt(True_shuffle_dmv_0.4_square)

n = 800
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

##compare theoretical variance with real variance
the_var = 0.8^2/12*( (num_state+1)/(num_state-1) )
print(paste('theretical variacne = ', the_var , 'sample variance = ',var(initila_state_all_nodes) )  )

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
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~shuffle_X(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) )
}

## average degree plots for Atlanta
plot(1:m, unlist(df$avg_edges))


## controll chart plots for Atlanta
diffA <- NULL
for (i in 2:(m)) {
  Ai = as.matrix(df$g[[i]])
  Ai_1 = as.matrix(df$g[[i-1]])
  diffA[i] <-  sum((Ai - Ai_1)^2)
  
}
plot(1:(m-1), diffA[-1])


D2 <- getD(df$Xt)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)

# finite ZC quantity using true Xt
hist(D2/D2_shuffle_1 )

D2_shuffle_0.4 <- getD(df$shuffle_Xt_0.4)


## compare finite behavior with infinity 
summary(as.vector((D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv))
summary(as.vector((D2-True_dmv)/True_dmv))
summary(as.vector((D2_shuffle_0.4-True_shuffle_dmv_0.4)/True_shuffle_dmv_0.4))

D2_shuffle_0.4[1,4]
True_shuffle_dmv_0.4[1,4]

doMDS(D2, doplot=TRUE)
doMDS(D2^2, doplot=TRUE)
doMDS(True_dmv^2, doplot=TRUE)

doMDS(True_shuffle_dmv, doplot=TRUE)
doMDS(D2_shuffle_1, doplot=TRUE)

doMDS(True_shuffle_dmv_0.4, doplot=TRUE)

doMDS(True_shuffle_dmv_0.4^2, doplot=TRUE)
#note above is not the same as below 
doMDS(0.6*True_dmv_square^2 + 0.4*True_shuffle_dmv_square^2, doplot=TRUE)


doMDS(D2_shuffle_0.4, doplot=TRUE)





