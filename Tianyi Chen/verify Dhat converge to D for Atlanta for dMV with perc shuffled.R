

vectorized_update_function <- function(current_position, p, delta) {
  # Initialize the next_position vector with the current positions
  n <- length(current_position)
  next_position <- current_position
  
  # --- Identify the indices for each condition ---
  # Use a small tolerance for floating-point comparisons
  epsilon <- 1e-10
  
  # Indices for nodes at the left boundary
  left_idx <- abs(current_position - 0.1) < epsilon
  
  # Indices for nodes at the right boundary
  right_idx <- abs(current_position - 0.9) < epsilon
  
  # Indices for nodes in the middle
  middle_idx <- (current_position > 0.1 + epsilon) & (current_position < 0.9 - epsilon)
  
  # --- Vectorized Jumps for each group ---
  
  # 1. Update Left Boundary Nodes
  if (any(left_idx)) {
    # Decide which nodes will jump
    num_left <- sum(left_idx)
    jump_decision <- runif(num_left) < p
    
    # Apply the jump only to those nodes
    next_position[left_idx][jump_decision] <- 0.1 + delta
  }
  
  # 2. Update Right Boundary Nodes
  if (any(right_idx)) {
    # Decide which nodes will jump
    num_right <- sum(right_idx)
    jump_decision <- runif(num_right) < p
    
    # Apply the jump
    next_position[right_idx][jump_decision] <- 0.9 - delta
  }
  
  # 3. Update Middle Nodes
  if (any(middle_idx)) {
    # Generate one set of random numbers for all middle nodes
    num_middle <- sum(middle_idx)
    u <- runif(num_middle)
    
    # Create logical vectors for jumping up or down
    jump_up <- u < p
    jump_down <- (u >= p) & (u < 2 * p)
    
    # Apply the updates using the logical vectors
    next_position[middle_idx][jump_up] <- next_position[middle_idx][jump_up] + delta
    next_position[middle_idx][jump_down] <- next_position[middle_idx][jump_down] - delta
  }
  
  return(next_position)
}


set.seed(1)
p = 0.05
q = 0.45
tmax = m = 30
tstar = 15

c=(.9-0.1)
num_state = 50
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

xt[,1]=initila_state_all_nodes
# Optimized simulation loop
for (j in 2:tstar) {
  xt[, j] <- vectorized_update_function(xt[, j - 1], p, delta)
}

for (j in (tstar + 1):(m + 1)) {
  xt[, j] <- vectorized_update_function(xt[, j - 1], q, delta)
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
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~shuffle_X_optimized(., perc)) ) ##this function only shuffle fix dn 
  #%>%
  #mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat,~shuffle_X(., perc)) )
}


D2 <- getD(df$Xt)
D2_shuffle_0.5 <- getD(df$shuffle_Xt_0.5)
D2_shuffle_1 <- getD(df$shuffle_Xt_1)

summary(na.omit(as.vector(abs(D2_shuffle_0.5-True_shuffle_dmv_0.5)/True_shuffle_dmv_0.5)))

#(D2_shuffle_0.5[1,2]-True_shuffle_dmv_0.5[1,2])/True_shuffle_dmv_0.5[1,2]
summary(na.omit(as.vector(abs(D2-True_dmv)/True_dmv)))
summary(na.omit(as.vector(abs(D2_shuffle_1-True_shuffle_dmv)/True_shuffle_dmv)))













##########------------------------

#edMV_0.5 = sqrt(mean((df$shuffle_Xt_0.5[[1]]-df$shuffle_Xt_0.5[[2]])^2))
#edMV_0.5
#edMV_1 = sqrt(mean((df$shuffle_Xt_1[[1]]-df$shuffle_Xt_1[[2]])^2))
#edMV = sqrt(mean((df$Xt[[1]]-df$Xt[[2]])^2))

#sqrt(0.5*edMV_1^2+0.5*edMV^2)
#edMV_0.5


True_dmv[1,2]

True_shuffle_dmv[1,2]


D2_shuffle_0.5[1,2]
True_shuffle_dmv_0.5[1,2]


(D2_shuffle_0.5[1,2]-True_shuffle_dmv_0.5[1,2])/True_shuffle_dmv_0.5[1,2]


