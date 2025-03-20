



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
  labs(y='MSE', x='MDS embedding dim d for the (d -> 1)-iso-mirror', color='Type', linetype='Type') +
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
