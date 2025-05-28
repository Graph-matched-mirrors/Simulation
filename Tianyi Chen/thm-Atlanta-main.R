## Theorem of Atlanta main:


p = 0
q = 0.5
m = 40
tstar = 20

c=(.9-0.1)
num_state = 200
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically dMV result.
mds = doMDS(True_dmv_square, doplot=TRUE)
psi = mds$mds[,1]

##with time being i/m here is the slope
diff(num_state*(num_state-1)/(2*c^2*m) * psi)/(1/m)

## close to 0.4 and 0.2





True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )