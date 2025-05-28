

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
