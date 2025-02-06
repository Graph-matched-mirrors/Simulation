pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)

##generate RDPG from given latent position X
rdpg.sample <- function(X, rdpg_scale=FALSE) {
  P <- X %*% t(X)
  if (rdpg_scale) {
    P <- scales::rescale(P)
  }
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph_from_adjacency_matrix(A,"undirected"))
}

##ASE for a network A with embedding dimension d
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  #    require(irlba)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- svds(A,k=d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}

## Procruste/i.e gain W in the estimated d_MV distance. Note in our case latent positions are 1d so this is trivial with w=1or -1
procrustes2 <- function(X, Y) {
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}


## Get distance matrix
getD <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
  if (k==0) {
    ind <- 1:n
  } else {
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    #cat("i = ", i, ", j = ", j, "\n")
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
    }
    
    D <- norm(Xhati - Xhatj, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  #as.dist(D2)
  D2 <- sqrt(D2)
  D2
}

## Apply CMDS on distance matrix 
doMDS <- function(D, doplot=TRUE)
{
  tmax <- m <- nrow(D)
  mds <- cmdscale(D, m-1)
  df.mds <- tibble(ind=1:tmax, time=sprintf("%2d",1:tmax), x=mds[,1], y=mds[,2], z=mds[,3], w=mds[,4])
  
  
  
  if (doplot) {
    
    library('gridExtra')
    
    # Convert the result of apply(mds, 2, sd) into a data frame for ggplot
    sd_data <- data.frame(
      dimension = 1:ncol(mds),  # Assumes mds is a matrix or data frame
      column_sd = apply(mds, 2, sd)
    )
    
    # First plot: standard deviation by dimension (converted to ggplot)
    p0 <- ggplot(sd_data, aes(x = dimension, y = column_sd)) +
      geom_point() +
      geom_line() +
      labs(x = "dimension", y = "column stdev") +
      theme_minimal()
    
    # Second plot: MDS1 vs time
    p1 <- df.mds %>%
      ggplot(aes(x=ind, y=x, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds1")
    
    # Third plot: MDS2 vs time
    p2 <- df.mds %>%
      ggplot(aes(x=ind, y=y, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds2")
    
    # Fourth plot: MDS3 vs time
    p3 <- df.mds %>%
      ggplot(aes(x=ind, y=z, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds3")
    
    # Fifth plot: MDS1 vs MDS2
    p4 <- df.mds %>%
      ggplot(aes(x=x, y=y, color=time)) +
      geom_point(size=3) +
      geom_label_repel(aes(label=time), size=2) +
      theme(legend.position = "none") +
      labs(x="mds1", y="mds2")
    
    # Arrange all four ggplots in a 1x4 grid
    grid.arrange(p0, p1, p2, p3, ncol=4)
    
    
  }
  
  return(list(mds=mds, df.mds=df.mds))
}

#apply ISOMAP on the CMDS result with chosen dimension mdsd from CMDS step defaultly it always embeds to 1 
doIso <- function(df.mds, mdsd=2, isod=1, doplot=F, AtlantaFlag = F)
{
  mds <- df.mds$mds
  df.iso <- NULL
  dis <- vegdist(mds[,1:mdsd,drop=F], "euclidean")
  if(AtlantaFlag){
    dis <- dis^2
  }
  knn <- 1
  success <- FALSE
  while(!success) {
    tryCatch({
      iso = isomap(dis, k=knn, ndim=isod, path="shortest")$points
      success <- TRUE
    },
    error = function(e) {
      knn <<- knn + 1
    })
  }
  iso2 <- tibble(iso=iso[,1]) %>% mutate(i=1:nrow(mds), ind=df.mds$ind, time=df.mds$time, knn=knn)
  df.iso <- rbind(df.iso, cbind(iso2, mdsd=mdsd))
  df.iso <- df.iso %>% group_by(mdsd) %>% mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  if (doplot) {
    p <- df.iso %>% filter(mdsd==mdsd) %>%
      ggplot(aes(x=ind, y=iso, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      theme(legend.position = "none") +
      labs(x="time", y="isomap embedding") +
      # scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short()) +
      theme(axis.text.x=element_text(hjust=0.7))
    # theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3),
    #       axis.text.y = element_text(size = 12),
    #       axis.title = element_text(size = 14, face="bold"))
    #    p <- p + scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short())
    print(p)
    
    df.isok <- df.iso %>% filter(mdsd==mdsd) #%>% mutate(date2 = format(ymd(paste0(date,"-01")),"%m/%y"))
    row.names(df.isok) <- df.isok$time
    fit <- lm(iso ~ i, data=df.isok)
    # print(tidy(fit))
    # print(glance(fit))
    myfor <- augment(fit)
    myfor2 <- myfor %>% mutate(date=.rownames,
                               ranks = rank(.sigma),
                               mycol=sprintf("%2d",rank(.fitted)))
    p <- myfor2 %>%
      ggplot(aes(.fitted, .resid)) +
      geom_point(aes(color=mycol)) +
      geom_hline(yintercept = 0, linetype="dashed", color="grey") +
      geom_smooth(method="loess", se=FALSE) +
      labs(x="Fitted Values", y="Residuals") +
      theme(legend.position = "none",
            axis.title = element_text(size=14, face="bold"))
    p <- p + geom_label_repel(aes(label=date), data=myfor2 %>% filter(ranks %in% 1:3))
    print(p)
  }
  
  return(df.iso)
}



## Implementation of the 3rd step in Algorithm 2 by recasting it as a linear programming problem.
## This function will return the objective function value Sk in the paper for a given change point t 
linf_cp=function(t,y,cp){
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp
  
  #XL; y[1:nl]
  
  XL2=XL
  XL2[,2:3]=-XL[,2:3]
  
  #rbind(XL,XL2); c(y[1:nl],-y[1:nl])
  
  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp
  
  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]
  
  
  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)
  
  library(lpSolveAPI)
  lprec <- make.lp(0,4)
  set.objfn(lprec,c(1,0,0,0))
  for (i in 1:(nrow(X)) ) {
    add.constraint(lprec, X[i,], ">=", Y[i])
  }
  
  set.bounds(lprec, lower = c(0,-Inf,-Inf,-Inf), columns = c(1,2,3,4))
  ColNames <- c('Z', "alpha", "bl","br")
  dimnames(lprec)[[2]] <- ColNames
  solve(lprec)
  return(get.variables(lprec))
}



shuffle_X <- function(X,del){
  n=nrow(X)
  dn=floor(del*n)
  permu_vec=sample(1:dn)
  random_perm=diag(dn)[permu_vec,]
  a=bdiag(diag(n-dn), random_perm)
  random_perm=as.matrix(a)
  return(random_perm %*% X)
}



## Implementation of the 3rd step in Algorithm 2 by recasting it as a linear programming problem.
## This function will return the objective function value Sk in the paper for a given change point t 
linf_cp=function(t,y,cp){
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp
  
  #XL; y[1:nl]
  
  XL2=XL
  XL2[,2:3]=-XL[,2:3]
  
  #rbind(XL,XL2); c(y[1:nl],-y[1:nl])
  
  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp
  
  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]
  
  
  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)
  
  library(lpSolveAPI)
  lprec <- make.lp(0,4)
  set.objfn(lprec,c(1,0,0,0))
  for (i in 1:(nrow(X)) ) {
    add.constraint(lprec, X[i,], ">=", Y[i])
  }
  
  set.bounds(lprec, lower = c(0,-Inf,-Inf,-Inf), columns = c(1,2,3,4))
  ColNames <- c('Z', "alpha", "bl","br")
  dimnames(lprec)[[2]] <- ColNames
  solve(lprec)
  return(get.variables(lprec))
}


linf_error=function(x, tmax){
  obf=NULL
  for (nk in 2:(tmax-1)) { ## find the point which minimize the obj func Sk, that is the change point 
    obf[nk]=linf_cp(1:tmax,x,nk)[1]
  }
  ecp=min(which(obf==min(obf[-1])))
  return((ecp-tstar)/tmax)
}


jump_Lbd=function(cp,p, delta){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=delta
  }
  return(np)
}


jump_Rbd=function(cp,p, delta){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=1-delta
  }
}

jump_middle=function(cp,p, delta){
  u=runif(1)
  if(u<p){
    
    np=cp+delta
    
  }
  
  if(p<u & u<2*p){
    np=cp-delta
  }
  
  if(u>2*p){
    np=cp
  }
  return(np)
}



update_function=function(current_position,p , delta){
  if( abs(current_position-0)<10^(-2) ){
    next_position=jump_Lbd(current_position,p,delta)
  }
  if(current_position> (0+10^(-2)) & current_position <(1-10^(-2)) ){
    next_position=jump_middle(current_position,p,delta)
  }
  if(abs(current_position-1)<10^(-2)){
    next_position=jump_Rbd(current_position,p,delta)
  }
  return(next_position)
}


shuffle_perc_graph <- function(A, del, n){
  G=as.matrix(A)
  print
  dn = floor(n*del)
  if(dn == 0)
    return(A)
  else{
    permu_vec=sample(1:dn)
    random_perm=diag(dn)[permu_vec,]
    a=bdiag(diag(n-dn), random_perm)
    random_perm = as.matrix(a)
    permu_G=as.matrix(random_perm%*%G%*%t(random_perm))
    G_graph=graph_from_adjacency_matrix(permu_G,mode ="undirected")
    return(G_graph)
  }
}

paired_error_in_shuffling_once <- function(n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = c(0.1,0.2)){
  
  xt=matrix(0,nrow = n, ncol = m+1)
  initial_state_all_nodes=sample(seq(0,1,by=delta),n,replace = TRUE)
  xt[,1]=initial_state_all_nodes
  for (i in 1:n) {
    for (j in 2:(tstar)) {
      xt[i,j]=update_function(xt[i,j-1],p,delta)
    }
    for (j in ((tstar+1):(m+1)) ) {
      xt[i,j]=update_function(xt[i,j-1],q,delta)
    }
  }
  
  
  
  
  df <- tibble(time=1:m) %>%
    mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.)))%>%
    mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  
  
  for(perc in del){
    df <- df %>%
      mutate(!!paste0("shuffle_", perc) := map(g,~shuffle_perc_graph(., perc, n))) %>%
      mutate(!!paste0("xhat_", perc) := map(!!sym(paste0("shuffle_", perc)), function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  }
  
  ### make graph
  
  D2=getD(df$xhat)
  df.mds <- doMDS(D2,doplot = FALSE)
  df.iso <- doIso(df.mds, mdsd=10, AtlantaFlag=T)
  errors <- NULL
  errors[1] <- linf_error(df.iso$iso, m)
  i <- 2
  for(perc in del){
    D2_shuffle=getD(df[[paste0("xhat_", perc)]])
    df.mds_shuffle <- doMDS(D2_shuffle,doplot = FALSE)
    df.iso_shuffle <- doIso(df.mds_shuffle, mdsd=10, AtlantaFlag=T)
    errors[i] <- linf_error(df.iso_shuffle$iso, m)
    i <- i + 1
  }
  print(paste(q,Sys.time()))
  errors
}
paired_error_in_shuffling <- function(nmc = 50, n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = 0.1){
  mc_errors <- sapply(1:nmc, function(i) paired_error_in_shuffling_once(n, p, q, m, delta, tstar, del))
  row_mse <- apply(mc_errors, 1, function(row) mean(abs(row)^2)) # mean(abs(errors[row,])^2)
  row_sds <- apply(mc_errors, 1, function(row) sd(abs(row)^2) / sqrt(nmc)) # sd(abs(errors[row,])^2) / sqrt(nmc)
  cbind(row_mse, row_sds)
}


