## for appendix

pacman::p_load(segmented, igraph, irlba ,RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)
library(ggplot2)
library(dplyr)

## from below you see how it gets more linear 
p = 0.1
q = 0.4

tmax = m = 40
tstar = m/2

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)


## 1, iso+mds with mdsd large can get you picewise linear(no entrywise square the dMV!)
True_D_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_D = sqrt(True_D_square)

Atlanta_mds = doMDS(True_D, doplot=TRUE)

#ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
#plot(1:(m-1), ev_in_mds, main = 'scree plot') ## elbow actually says you need at least d = 3
#elbow = getElbows( ev_in_mds )
#elbow
#doIso(Atlanta_mds$mds, mdsd = elbow[1], doplot = T)
#doIso(Atlanta_mds$mds, mdsd = elbow[2], doplot = T)
#doIso(Atlanta_mds$mds, mdsd = elbow[3], doplot = T)
#doIso(Atlanta_mds$mds, mdsd = 10, doplot = T)


True_shuffled_D_square =true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffled_D = sqrt(True_shuffled_D_square)


del = seq(0,1,by=0.05)

for (alpha in del) {
  nm <- paste0("True_shuffle_dmv_", alpha)
  assign(nm,
         sqrt((1 - alpha) * True_D^2 + alpha * True_shuffled_D^2),
         envir = .GlobalEnv)
}


#MDS_True_shuffle_dmv_0.1 = doMDS(True_shuffle_dmv_0.1,doplot = T)
#ev_in_mds = apply(MDS_True_shuffle_dmv_0.1$mds, 2, sd)
#plot(1:(m-1), ev_in_mds, main = 'scree plot') 
#elbow = getElbows( ev_in_mds )
#elbow
#doIso(MDS_True_shuffle_dmv_0.1$mds, mdsd = elbow[3], doplot = T)
#ture_iso = doIso(MDS_True_shuffle_dmv_0.1$mds, mdsd = 13, doplot = F)

# Check what graphics device is active
dev.cur()

# If it's not the RStudio device, reset it
if(dev.cur() != 2) {  # 2 is usually RStudio's graphics device
  dev.off()
  dev.new()
}
MDS_True_shuffle_dmv_0.2 = doMDS(True_shuffle_dmv_0.2,doplot = T)
ev_in_mds = apply(MDS_True_shuffle_dmv_0.2$mds, 2, sd)
#plot(1:(m-1), ev_in_mds, main = 'scree plot') 
elbow = getElbows( ev_in_mds )
elbow
iso = doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = elbow[3], doplot = T)
plot(1:tmax, iso$iso )
doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = 8, doplot = T)

MDS_True_shuffle_dmv_0.5 = doMDS(True_shuffle_dmv_0.5,doplot = T)
ev_in_mds = apply(MDS_True_shuffle_dmv_0.2$mds, 2, sd)
elbow = getElbows( ev_in_mds )

doIso(MDS_True_shuffle_dmv_0.5$mds, mdsd = 7, doplot = T)


MDS_True_shuffle_dmv_0.1_square = doMDS(True_shuffle_dmv_0.1^2,doplot = T)
#ev_in_mds = apply(MDS_True_shuffle_dmv_0.1_square$mds, 2, sd)
#plot(1:(m-1), ev_in_mds, main = 'scree plot') 
#elbow = getElbows( ev_in_mds )
#elbow

set.seed(2)
n = 300
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

mean(initila_state_all_nodes)
the_var = c^2/12*( (num_state+1)/(num_state-1) )
print(paste('theretical variacne = ', the_var , 'sample variance = ',var(initila_state_all_nodes) )  )
sqrt(2*the_var) #is the entry of true 100% shuffling dMV
True_shuffled_D[1,2]

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
    #mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~ shuffle_X_optimized(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X_optimized(., perc)) )
}


D2_shuffle_0.05_hat <- getD(df$shuffle_Xhat_0.05)
MDS_shuffle_0.05_hat <- doMDS(D2_shuffle_0.05_hat, doplot = T)
find_slope_changepoint_with_plot(MDS_shuffle_0.05_hat$mds[,1])

MDS_shuffle_square_0.05_hat <- doMDS(D2_shuffle_0.05_hat^2, doplot = T)
find_slope_changepoint_with_plot(MDS_shuffle_square_0.05_hat$mds[,1])

ISOd = 10
ISO = doIso(MDS_shuffle_0.05_hat$mds, mdsd = ISOd, doplot = F)
find_slope_changepoint_with_plot(ISO$iso)


ISOd = 3
ture_iso = doIso(MDS_True_shuffle_dmv_0.05$mds, mdsd = ISOd, doplot = F)
find_slope_changepoint_with_plot(ture_iso$iso)

MDS_True_shuffle_dmv_0.05 = doMDS(True_shuffle_dmv_0.05,doplot = T)
MDS_True_shuffle_square_dmv_0.05 = doMDS(True_shuffle_dmv_0.05^2,doplot = T)
find_slope_changepoint_with_plot(MDS_True_shuffle_square_dmv_0.05$mds[,1])




#MDS_True_shuffle_square_dmv_0.05$mds[,1]/MDS_True_shuffle_dmv_0.05$mds[,1] 
############# Begin to make figures


pdf("n=1000m=30N=50_shuffling_0.05.pdf", width = 9, height = 6)
par(mfrow = c(2,3))

plot(1:tmax, MDS_True_shuffle_square_dmv_0.05$mds[,1] ,col='black' ,ylim = c(-0.003,0.005) ,ylab = 'mirror',xlab = 'time'
     ,main = expression("1 d MDS " * d[MV]^2 * " on 5% shuffling") )
points(1:tmax,MDS_shuffle_square_0.05_hat$mds[,1],col = 'blue')


plot(1:tmax, MDS_shuffle_0.05_hat$mds[,1],col = 'blue',ylim = c(-0.03,0.04) , ylab = 'mirror',xlab = 'time'
     , main = expression("1 d MDS " * d[MV] * " on 5% shuffling"))
points(1:tmax, MDS_True_shuffle_dmv_0.05$mds[,1] ,col='black')



ISOd= 3
ture_iso = doIso(MDS_True_shuffle_dmv_0.05$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.05_hat$mds, mdsd = ISOd, doplot = F)

l = min(c( ture_iso$iso,ISO$iso))
u = max(c( ture_iso$iso,ISO$iso))
plot(1:tmax, ture_iso$iso,col='black',ylim=c(l,u), ylab = 'iso-mirror',
     main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 5 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')


ISOd= 5
ture_iso = doIso(MDS_True_shuffle_dmv_0.05$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.05_hat$mds, mdsd = ISOd, doplot = F)
l = min(c( ture_iso$iso,ISO$iso))
u = max(c( ture_iso$iso,ISO$iso))
plot(1:tmax, ture_iso$iso,col='black' ,ylim=c(l,u),ylab = 'iso-mirror'
     ,main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 5 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')

ISOd= 8
ture_iso = doIso(MDS_True_shuffle_dmv_0.05$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.05_hat$mds, mdsd = ISOd, doplot = F)

plot(1:tmax, ture_iso$iso,col='black' ,ylab = 'iso-mirror',
     main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 5 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')


ISOd= 10
ture_iso = doIso(MDS_True_shuffle_dmv_0.05$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.05_hat$mds, mdsd = ISOd, doplot = F)

plot(1:tmax, ture_iso$iso,col='black' ,ylab = 'iso-mirror',main = 
       bquote(
         ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 5 * "% shuffling"
       ))
points(1:tmax,ISO$iso,col = 'blue')

dev.off()
#n=1000m=30N=50_shuffling_0.2.pdf




MDS_True_shuffle_dmv_0.2 = doMDS(True_shuffle_dmv_0.2,doplot = F)
MDS_True_shuffle_square_dmv_0.2 = doMDS(True_shuffle_dmv_0.2^2,doplot = F)



D2_shuffle_0.2_hat <- getD(df$shuffle_Xhat_0.2)
MDS_shuffle_0.2_hat <- doMDS(D2_shuffle_0.2_hat, doplot = F)
MDS_shuffle_square_0.2_hat <- doMDS(D2_shuffle_0.2_hat^2, doplot = F)



MDS_True_shuffle_dmv_0.2$mds[,1]/MDS_True_shuffle_square_dmv_0.2$mds[,1]

############# Begin to make figures

pdf("n=1000m=30N=50_shuffling_0.2.pdf", width = 9, height = 6)

par(mfrow = c(2,3))

l = min(c( MDS_shuffle_square_0.2_hat$mds[,1],MDS_True_shuffle_square_dmv_0.2$mds[,1]))
u = max(c( MDS_shuffle_square_0.2_hat$mds[,1],MDS_True_shuffle_square_dmv_0.2$mds[,1]))
plot(1:tmax, MDS_True_shuffle_square_dmv_0.2$mds[,1] ,col='black',ylim = c(l,u), ylab = 'mirror',xlab = 'time'
     ,main = expression("1 d MDS " * d[MV]^2 * " on 20% shuffling") )
points(1:tmax,MDS_shuffle_square_0.2_hat$mds[,1],col = 'blue')



l = min(c( MDS_shuffle_0.2_hat$mds[,1],MDS_True_shuffle_dmv_0.2$mds[,1]))
u = max(c( MDS_shuffle_0.2_hat$mds[,1],MDS_True_shuffle_dmv_0.2$mds[,1]))
plot(1:tmax, MDS_shuffle_0.2_hat$mds[,1],col = 'blue',ylim = c(l,u) , ylab = 'mirror',xlab = 'time'
     , main = expression("1 d MDS " * d[MV] * " on 20% shuffling"))
points(1:tmax, MDS_True_shuffle_dmv_0.2$mds[,1] ,col='black')

ISOd= 3
ture_iso = doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.2_hat$mds, mdsd = ISOd, doplot = F)

l = min(c( ture_iso$iso,ISO$iso))
u = max(c( ture_iso$iso,ISO$iso))
plot(1:tmax, ture_iso$iso,col='black',ylim=c(l,u), ylab = 'iso-mirror',
     main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 20 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')


ISOd= 5
ture_iso = doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.2_hat$mds, mdsd = ISOd, doplot = F)

l = min(c( ture_iso$iso,ISO$iso))
u = max(c( ture_iso$iso,ISO$iso))
plot(1:tmax, ture_iso$iso,col='black',ylim=c(l,u), ylab = 'iso-mirror',
     main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 20 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')


ISOd= 8
ture_iso = doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.2_hat$mds, mdsd = ISOd, doplot = F)
l = min(c( ture_iso$iso,ISO$iso))
u = max(c( ture_iso$iso,ISO$iso))
plot(1:tmax, ture_iso$iso,col='black' ,ylim=c(l,u),ylab = 'iso-mirror'
     ,main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 20 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')

ISOd= 10
ture_iso = doIso(MDS_True_shuffle_dmv_0.2$mds, mdsd = ISOd, doplot = F)
ISO = doIso(MDS_shuffle_0.2_hat$mds, mdsd = ISOd, doplot = F)

plot(1:tmax, ture_iso$iso,col='black' ,ylab = 'iso-mirror',
     main = bquote(
       ISO ~ on ~ bold(.(as.character(ISOd))) ~ d ~ MDS ~ d[MV] ~ on ~ 20 * "% shuffling"
     ))
points(1:tmax,ISO$iso,col = 'blue')


dev.off()















MDS_True_shuffle_dmv_0.5 = doMDS(True_shuffle_dmv_0.5,doplot = T)
MDS_True_shuffle_square_dmv_0.5 = doMDS(True_shuffle_dmv_0.5^2,doplot = T)
D2_shuffle_0.5_hat <- getD(df$shuffle_Xhat_0.5)
MDS_shuffle_0.5_hat <- doMDS(D2_shuffle_0.5_hat, doplot = T)
MDS_shuffle_square_0.5_hat <- doMDS(D2_shuffle_0.5_hat^2, doplot = T)


plot(1:tmax, MDS_shuffle_0.5_hat$mds[,1],col = 'blue', main = '1d MDS dMV on 50% shuffling')
points(1:tmax, MDS_True_shuffle_dmv_0.5$mds[,1] ,col='black')

plot(1:tmax, MDS_shuffle_square_0.5_hat$mds[,1],col = 'blue', main = '1d MDS d^2MV on 50% shuffling')
points(1:tmax, MDS_True_shuffle_square_dmv_0.5$mds[,1] ,col='black')



D2_shuffle_0.1_hat <- getD(df$shuffle_Xhat_0.1)
MDS_shuffle_0.1_hat <- doMDS(D2_shuffle_0.1_hat, doplot = T)
MDS_shuffle_square_0.1_hat <- doMDS(D2_shuffle_0.1_hat^2, doplot = T)
plot(1:tmax, MDS_shuffle_0.1_hat$mds[,1], col = 'blue')
points(1:tmax, MDS_True_shuffle_dmv_0.1$mds[,1] ,col='black')


plot(1:tmax, MDS_shuffle_square_0.1_hat$mds[,1], col = 'blue')
points(1:tmax, MDS_True_shuffle_dmv_0.1_square$mds[,1],col='black')


ev_in_mds = apply(MDS_shuffle_0.1_hat$mds, 2, sd)
plot(1:(m-1), ev_in_mds, main = 'scree plot') 
elbow = getElbows( ev_in_mds )
elbow
ISO = doIso(MDS_shuffle_0.1_hat$mds, mdsd = 10, doplot = F)
ture_iso = doIso(MDS_True_shuffle_dmv_0.1$mds, mdsd = 10, doplot = F)

plot(1:tmax, ture_iso$iso,col='black')
points(1:tmax,ISO$iso,col = 'blue')


D2_shuffle_0.2_hat <- getD(df$shuffle_Xhat_0.2)
MDS_shuffle_0.2_hat <- doMDS(D2_shuffle_0.2_hat, doplot = F)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:length(ev_in_mds), ev_in_mds, main = 'scree plot') ## elbow actually says you only need d =1 
elbow = getElbows( ev_in_mds )

