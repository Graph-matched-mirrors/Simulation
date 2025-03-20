pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)

registerDoParallel(detectCores()-1)
set.seed(2)


p = 0.4
q = 0.2
tmax = m = 20
tstar = 10

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)

nmc=100
max_iter=30
n = 500

a <- Sys.time()
out_dd <- foreach (mc = 1:nmc) %dopar% {
  
  tmp1 <- tmp2 <- tmp3 <- tmp4 <- rep(0,3)
  
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
    mutate(xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
    mutate(g = map(xt, ~rdpg.sample(.))) %>%
    mutate(shuffle_g =map(g, ~shuffle_graph(.)) )
  
  df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
    mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
    mutate(shuffle_g_GM_alltoone = map(shuffle_g, ~graph_mathing(shuffle_g[[1]],.x, max_iter ))) %>%
    mutate(shuffle_g_GM_pairwise = purrr::accumulate(shuffle_g[-1],
                                                     .f = function(acc, curr) graph_mathing(acc, curr, max_iter),
                                                     .init = shuffle_g[[1]]))%>%  # Sequential matching with correction
    mutate(Xhat_shuffle_GM_alltoone = map(shuffle_g_GM_alltoone, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>%
    mutate(Xhat_shuffle_GM_pairwise = map(shuffle_g_GM_pairwise, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  
  D2 <- getD(df$Xhat) 
  D2_shuffle <- getD(df$Xhat_shuffle)
  D2_shuffle_GM_alltoone <- getD(df$Xhat_shuffle_GM_alltoone)
  D2_shuffle_GM_pairwise <- getD(df$Xhat_shuffle_GM_pairwise) 
  
  df.mds <- doMDS(D2^2,doplot = F)
  df.mds_shuffle <- doMDS(D2_shuffle^2,doplot = F)
  df.mds_shuffle_GM_alltoone <- doMDS(D2_shuffle_GM_alltoone^2,doplot = F)
  df.mds_shuffle_GM_pairwise <- doMDS(D2_shuffle_GM_pairwise^2,doplot = F)
  
  mds <- df.mds$mds
  mds_shuffle <- df.mds_shuffle$mds
  mds_shuffle_GM_pairwise <- df.mds_shuffle_GM_pairwise$mds
  mds_shuffle_GM_alltoone <- df.mds_shuffle_GM_alltoone$mds
  
  
  for (dd in 1:3){
    df.iso <- doIso(df.mds$mds, mdsd=dd)$iso
    df.iso_shuffle <- doIso(df.mds_shuffle$mds, mdsd=dd)$iso
    df.iso_shuffle_GM_alltoone <- doIso(df.mds_shuffle_GM_alltoone$mds, mdsd=dd)$iso
    df.iso_shuffle_GM_pairwise <- doIso(df.mds_shuffle_GM_pairwise$mds, mdsd=dd)$iso
    
    tmp1[dd]=find_slope_changepoint_with_plot(df.iso, doplot = F)$error
    tmp2[dd]=find_slope_changepoint_with_plot(df.iso_shuffle, doplot = F)$error
    tmp3[dd]=find_slope_changepoint_with_plot(df.iso_shuffle_GM_alltoone, doplot = F)$error
    tmp4[dd]=find_slope_changepoint_with_plot(df.iso_shuffle_GM_pairwise, doplot = F)$error
  }
  list(tmp1,tmp2,tmp3,tmp4)
}
print(Sys.time() - a)


msehat1=Reduce('cbind', lapply(out_dd, "[[", 1)) ## this will summarize tmp1 

sm_mse1=as.data.frame(matrix(0,3,4))
sm_mse1[,2]=apply(abs(msehat1)^2, 1, mean)
sm_mse1[,1]=apply(abs(msehat1)^2, 1, mean)-apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse1[,3]=apply(abs(msehat1)^2, 1, mean)+apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse1[,4]=1:3

sm_mse1
msehat2=Reduce('cbind', lapply(out_dd, "[[", 2))

sm_mse2=as.data.frame(matrix(0,3,4))
sm_mse2[,2]=apply(abs(msehat2)^2, 1, mean)
sm_mse2[,1]=apply(abs(msehat2)^2, 1, mean)-apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse2[,3]=apply(abs(msehat2)^2, 1, mean)+apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse2[,4]=1:3

sm_mse2

msehat3=Reduce('cbind', lapply(out_dd, "[[", 3))

sm_mse3=as.data.frame(matrix(0,3,4))
sm_mse3[,2]=apply(abs(msehat3)^2, 1, mean)
sm_mse3[,1]=apply(abs(msehat3)^2, 1, mean)-apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse3[,3]=apply(abs(msehat3)^2, 1, mean)+apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse3[,4]=1:3
sm_mse3

msehat4=Reduce('cbind', lapply(out_dd, "[[", 4))

sm_mse4=as.data.frame(matrix(0,3,4))
sm_mse4[,2]=apply(abs(msehat4)^2, 1, mean)
sm_mse4[,1]=apply(abs(msehat4)^2, 1, mean)-apply(abs(msehat4)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse4[,3]=apply(abs(msehat4)^2, 1, mean)+apply(abs(msehat4)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse4[,4]=1:3
sm_mse4


# Combine the data frames and add a new column 'type' to distinguish them
sm_mse1$type <- "True 1-1"
sm_mse2$type <- "shuffled"
sm_mse3$type <- "shuffled then GM (all to one)"
sm_mse4$type <- "shuffled then GM (pairwise)"
sm_mse_all <- rbind(sm_mse1, sm_mse2, sm_mse3, sm_mse4)

# Plot

library(ggplot2)
plottt <- ggplot(sm_mse_all, aes(x=V4, y=V2, color=type, linetype=type)) + 
  geom_line() +
  geom_errorbar(aes(ymin=V1, ymax=V3)) +
  scale_x_continuous(breaks = 1:3) +
  labs(y='relative MSE', x='MDS embedding dim d for the (d -> 1)-iso-mirror', color='Type', linetype='Type') +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"))


plottt <- plottt + labs(title = paste("p =", p, "q =", q,"m=",m, "n =", n, "nmc =", nmc, 'max_iter=', max_iter))
plottt <- plottt + theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25, face = "bold"),
  legend.text = element_text(size = 25),
  legend.title = element_text(size = 25, face = "bold"),
  plot.title = element_text(size = 20, face = "bold") 
)



print(plottt)
