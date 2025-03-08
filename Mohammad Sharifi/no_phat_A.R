
pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)

registerDoParallel(detectCores()-1)
set.seed(2)


n=300
p <- 0.4
q <- 0.3
m <- tmax <- 10
delta <- 0.1 
tstar <- tmax/2
nmc=100
max_iter=30
a <- Sys.time()
out_dd <- foreach (mc = 1:nmc) %dopar% {
  
  tmp1 <- tmp2 <- tmp3 <- tmp4 <- rep(0,3)
  df <- doSim_shuffle_Atlanta(n,tmax,delta,p,q,tstar)
  df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
    mutate (num_e=map(g, ~sum(as.matrix(.x)))) %>%
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
  
  df.mds <- doMDS(D2,doplot = F)
  df.mds_shuffle <- doMDS(D2_shuffle,doplot = F)
  df.mds_shuffle_GM_alltoone <- doMDS(D2_shuffle_GM_alltoone,doplot = F)
  df.mds_shuffle_GM_pairwise <- doMDS(D2_shuffle_GM_pairwise,doplot = F)


  for (dd in 1:3){
    df.iso <- doIso(df.mds, mdsd=dd)
    df.iso_shuffle <- doIso(df.mds_shuffle, mdsd=dd)
    df.iso_shuffle_GM_alltoone <- doIso(df.mds_shuffle_GM_alltoone, mdsd=dd)
    df.iso_shuffle_GM_pairwise <- doIso(df.mds_shuffle_GM_pairwise, mdsd=dd)
 
    tmp1[dd]=L2_changepoint_error(df.iso)
    tmp2[dd]=L2_changepoint_error(df.iso_shuffle)
    tmp3[dd]=L2_changepoint_error(df.iso_shuffle_GM_alltoone)
    tmp4[dd]=L2_changepoint_error(df.iso_shuffle_GM_pairwise)
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


plottt <- plottt + labs(title = paste("p =", p, "q =", q,"m=",m, "n =", n, "nmc =", nmc, 'max_iter=', 100))
plottt <- plottt + theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25, face = "bold"),
  legend.text = element_text(size = 25),
  legend.title = element_text(size = 25, face = "bold"),
  plot.title = element_text(size = 20, face = "bold") 
)



print(plottt)



