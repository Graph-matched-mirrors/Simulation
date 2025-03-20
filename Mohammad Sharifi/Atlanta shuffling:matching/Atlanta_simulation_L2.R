pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
library(ggrepel)
library(purrr)
library(irlba)
registerDoParallel(detectCores()-1)



set.seed(2)

p = 0.4
q = 0.2
tmax = m = 20
tstar = 10

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)

del = c(1) ## shuffling percentage 


## check piecewise linearity using true dMV^2
True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)
Atlanta_mds = doMDS(True_dmv^2, doplot=TRUE)

n = 500


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



###### here is to verify at least xhat and xhat-shuffled has results that make sense 
df <- tibble(time=1:m) %>%
  mutate(xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
  mutate(g = map(xt, ~rdpg.sample(.))) %>%
  mutate(shuffle_g =map(g, ~shuffle_graph(.)) ) %>%
  mutate(avg_edges = map(g,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
  mutate(shuffle_g =map(g, ~shuffle_graph(.)) ) %>%
  mutate(avg_edges = map(g,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%  
  mutate(xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) 

D2 <- getD(df$xhat)
D2_shuffle <- getD(df$xhat_shuffle)

df.mds <- doMDS(D2^2, doplot = F)
df.mds_shuffle <- doMDS(D2_shuffle^2, doplot = F)

find_slope_changepoint_with_plot(df.mds$mds[,1])
find_slope_changepoint_with_plot(df.mds_shuffle$mds[,1])
###########

max_iter = 10




df <- tibble(time=1:m) %>%
  mutate(xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
  mutate(g = map(xt, ~rdpg.sample(.))) %>%
  mutate(shuffle_g =map(g, ~shuffle_graph(.)) ) %>%
  mutate(avg_edges = map(g,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%  
  mutate(xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
  mutate(shuffle_g_GM_pairwise = purrr::accumulate(shuffle_g[-1],
                                          .f = function(acc, curr) graph_mathing(acc, curr, max_iter ),
                                          .init = shuffle_g[[1]]))%>%  # Sequential matching with correction
  mutate(xhat_shuffle_GM_pairwise = map(shuffle_g_GM_pairwise, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
  mutate(shuffle_g_GM_all_to_one = map(shuffle_g, ~graph_mathing(shuffle_g[[1]],.x,max_iter ))) %>%
  mutate(xhat_shuffle_GM_all_to_one = map(shuffle_g_GM_all_to_one, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D2 <- getD(df$xhat)
D2_shuffle <- getD(df$xhat_shuffle)
D2_shuffle_GM_pairwise <- getD(df$xhat_shuffle_GM_pairwise)
D2_shuffle_GM_all_to_one <- getD(df$xhat_shuffle_GM_all_to_one)

# Perform MDS (with plotting enabled, make sure you entrywise square it!!! For picewise linearity)
df.mds <- doMDS(D2^2, doplot = TRUE)
df.mds_shuffle <- doMDS(D2_shuffle^2, doplot = TRUE)
df.mds_shuffle_GM_pairwise <- doMDS(D2_shuffle_GM_pairwise^2, doplot = TRUE)
df.mds_shuffle_GM_all_to_one <- doMDS(D2_shuffle_GM_all_to_one^2, doplot = TRUE)

# Extract the MDS coordinates
mds <- df.mds$mds
mds_shuffle <- df.mds_shuffle$mds
mds_shuffle_GM_pairwise <- df.mds_shuffle_GM_pairwise$mds
mds_shuffle_GM_all_to_one <- df.mds_shuffle_GM_all_to_one$mds



find_slope_changepoint_with_plot(mds[,1])
find_slope_changepoint_with_plot(mds_shuffle[,1])
find_slope_changepoint_with_plot(mds_shuffle_GM_pairwise[,1])
find_slope_changepoint_with_plot(mds_shuffle_GM_all_to_one[,1])



