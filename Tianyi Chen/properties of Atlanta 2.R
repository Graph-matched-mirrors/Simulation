pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
library(ggrepel)
registerDoParallel(detectCores()-1)

# Here we summarize properties we observe for Atlanta model:


## when n = infty: xhat = x and dmvhat = dmv 

p = 0.4
q = 0.2
m = 40
tstar = 20

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)

del = c(1, 0.4) ## shuffling percentage 


## 1, iso+mds with mdsd large can get you picewise linear(no entrywise square the dMV!)
True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
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
num_state = 3
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

Atlanta_mds = doMDS(True_dmv, doplot=TRUE)
Atlanta_mds = doMDS(True_dmv^2, doplot=TRUE)
Atlanta_mds = doMDS(True_dmv^8, doplot=TRUE)

ev_in_mds = apply(Atlanta_mds$mds, 2, sd)
plot(1:length(ev_in_mds), ev_in_mds, main = 'scree plot') ## elbow actually says you only need d =1 
elbow = getElbows( ev_in_mds )


## 100% shuffling behavior

True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )


##c percent shuffling behavior: 



True_shuffle_dmv_0.4_square= 0.6*True_dmv_square + 0.4*True_shuffle_dmv_square
doMDS(True_shuffle_dmv_0.4_square, doplot=TRUE)
doMDS(sqrt(True_shuffle_dmv_0.4_square), doplot=TRUE)



doMDS(True_shuffle_dmv_square, doplot = T)
doMDS(sqrt(True_shuffle_dmv_square), doplot = T)

doMDS(True_dmv_square, doplot=TRUE)
doMDS(sqrt(True_dmv_square), doplot=TRUE)





## ZC quantity behavior, change num_state N and see behavior. 
num_state = 10
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_dmv=sqrt(True_dmv_square)

True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffle_dmv = sqrt(True_shuffle_dmv_square)

hist(True_dmv / True_shuffle_dmv, xlim = c(0, 1))




## finite sample behavior:
p = 0.4
q = 0.4
m = 4
tstar = 2

c=(.9-0.1)
num_state = 3
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )

True_shuffle_dmv_0.5_square= 0.5*True_dmv_square + 0.5*True_shuffle_dmv_square

True_shuffle_dmv_0.4_square= 0.6*True_dmv_square + 0.4*True_shuffle_dmv_square

#True_shuffle_dmv_0.11_square= (1/num_state^2)*True_dmv_square + (1-1/num_state^2)*True_shuffle_dmv_square

True_shuffle_dmv_0.09_square= (0.09)*True_dmv_square + (1-0.09)*True_shuffle_dmv_square


True_dmv=sqrt(True_dmv_square)
True_shuffle_dmv = sqrt(True_shuffle_dmv_square)
True_shuffle_dmv_0.4=sqrt(True_shuffle_dmv_0.4_square)
True_shuffle_dmv_0.5 = sqrt(True_shuffle_dmv_0.5_square)
#True_shuffle_dmv_0.11 = sqrt(True_shuffle_dmv_0.11_square)
True_shuffle_dmv_0.09 = sqrt(True_shuffle_dmv_0.09_square)


n = 1000
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
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))  
#%>%
#  mutate(g = map(Xt, ~rdpg.sample(.))) 
#%>%
#  mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
#  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
#del = c(1,0.4)
del = c(0.5,1)

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~optimized_shuffled_X(., perc)) ) 
    #%>%
    #mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X(., perc)) )
}
D2 <- getD(df$Xt)
D2_shuffle_0.5 <- getD(df$shuffle_Xt_0.5)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)

summary(na.omit(as.vector(abs(D2_shuffle_0.5-True_shuffle_dmv_0.5)/True_shuffle_dmv_0.5)))
summary(na.omit(as.vector(abs(D2-True_dmv)/True_dmv)))
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv)))

sqrt(D2_shuffle_1^2*0.5+0.5*D2^2)
D2_shuffle_0.5

summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv)))
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv_0.11)/True_shuffle_dmv_0.11)))
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv_0.09)/True_shuffle_dmv_0.09)))


df$Xt[[1]] - df$Xt[[2]]
i = 1
j = 4
stay1=which( (df$Xt[[i]] - df$shuffle_Xt_1[[i]]) == 0)
stay2=which( (df$Xt[[j]] - df$shuffle_Xt_1[[j]]) == 0)

length(stay1)/n


stayboth=intersect(stay1,stay2)

length(stayboth)/n

sqrt(sum((df$shuffle_Xt_1[[1]] - df$shuffle_Xt_1[[2]])^2)/n)
D2_shuffle_1[1,2]



sqrt(sum((df$shuffle_Xt_1[[1]][stayboth,] - df$shuffle_Xt_1[[2]][stayboth,] )^2)/length(stayboth))

notstay = setdiff(1:n,stayboth)

sum((df$shuffle_Xt_1[[1]][notstay,] - df$shuffle_Xt_1[[2]][notstay ,] )^2)/length(notstay)


sqrt(sum((df$shuffle_Xt_1[[1]][notstay,] - df$shuffle_Xt_1[[2]][notstay ,] )^2)/length(notstay ))

D2[1,2]
True_dmv[1,2]


D2 <- getD(df$Xt)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv)))
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv_0.11)/True_shuffle_dmv_0.11)))

hist(D2/D2_shuffle_1, xlim = c(0,2), main = paste('N = ', num_state, 'n = ', n) )
D2/D2_shuffle_1


(True_shuffle_dmv-True_shuffle_dmv_0.11)/True_shuffle_dmv

D2
D2_shuffle_1



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
hist(D2/D2_shuffle_1, xlim = c(0,2), main = paste('N = ', num_state, 'n = ', n) )

# Load the library
library(pheatmap)

# Define a custom color gradient
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Generate heatmap with color legend and no clustering
pheatmap(D2 / D2_shuffle_1, cluster_rows = FALSE, cluster_cols = FALSE, 
         color = custom_colors,
         main = "Heatmap of D2/D2_shuffle_1 (No Clustering)",
         legend = TRUE)

## compare finite behavior with infinity using Xt!! Not xhat
summary(as.vector((D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv))
summary(as.vector((D2-True_dmv)/True_dmv))
#summary(as.vector((D2_shuffle_0.4-True_shuffle_dmv_0.4)/True_shuffle_dmv_0.4))

D2_shuffle_0.4[1,4]
True_shuffle_dmv_0.4[1,4]

doMDS(D2, doplot=TRUE)
doMDS(D2^2, doplot=TRUE)
doMDS(True_dmv^2, doplot=TRUE)

doMDS(True_shuffle_dmv, doplot=TRUE)
doMDS(D2_shuffle_1, doplot=TRUE)

doMDS(True_shuffle_dmv_0.4, doplot=TRUE)
doMDS(True_shuffle_dmv_0.11, doplot=TRUE)


doMDS(True_shuffle_dmv_0.4^2, doplot=TRUE)
#note above is not the same as below 
doMDS(0.6*True_dmv_square^2 + 0.4*True_shuffle_dmv_square^2, doplot=TRUE)


doMDS(D2_shuffle_0.4, doplot=TRUE)





