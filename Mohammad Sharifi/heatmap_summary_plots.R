m <- 40
tstar <- 20
Num_states <- m 
p <- 0.4
q <- seq(0, 0.5, by = 0.05)
d <- seq(0.05, 1, by = 0.05)
n <- 300
nmc <- 100

library(readr)
final_errors <- read_csv("/Users/mohammad/Documents/Simulation/Mohammad Sharifi/final_errors_London.cvs")
final_errors <- read_csv("/Users/mohammad/Documents/Simulation/Mohammad Sharifi/final_errors_Atlanta3.cvs")
final_errors <- read_csv("~/final_errors_Atlanta3.cvs")
View(final_errors)
final_mse=final_errors[,seq(1,10,by=2)]
final_sd=final_errors[,seq(2,11,by=2)]

final_mse$row_mse_0.4

heatmap(as.matrix(final_mse),scale = 'none')

heatmap(as.matrix(final_mse), Rowv = NA, Colv = NA,
        main = paste("n =",n, ", p =", p, ", nmc =", nmc),
        xlab = "q", ylab = "shuffling ratio", scale = "none") 

support = seq(2/m-0.5, (m-1)/m-0.5,by=1/m)
chance_level = sum(support^2/length(support))
print(chance_level)


library(ggplot2)
library(reshape2)
# Convert matrix to data frame for ggplot2
data_melt <- melt(as.matrix(final_mse))
data_melt_sd <- melt(as.matrix(final_sd))

shuffle_ratio=c(0,d)
data_melt$Var1 <- factor(data_melt$Var1, labels = shuffle_ratio)
data_melt$Var2 <- factor(data_melt$Var2, labels = q)

data_melt_sd$Var1 <- factor(data_melt_sd$Var1, labels = shuffle_ratio)
data_melt_sd$Var2 <- factor(data_melt_sd$Var2, labels = q)


colnames(data_melt)=c('shuffle ratio','q','mean')
colnames(data_melt_sd)=c('shuffle ratio','q','sd')

summary=cbind(data_melt,data_melt_sd$sd)

summary$lower=summary$mean-summary$`data_melt_sd$sd`*1.96
summary$upper=summary$mean+summary$`data_melt_sd$sd`*1.96


library(dplyr)  # For filtering the data

# Filter the data for specific shuffle ratios
filtered_summary <- summary %>%
  filter(`shuffle ratio` %in% c(0, 0.25, 0.5, 0.75,1))

filtered_summary <- summary %>%
  filter(q %in% c(0.1,0.2, 0.35, 0.4, 0.5))

plottt <- ggplot(filtered_summary, aes(x = `shuffle ratio`, y = mean, color = q, linetype=q ,group = q)) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(y = 'MSE', x = 'Shuffle ratio', color = 'q') +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

print(plottt)
plottt <- ggplot(filtered_summary, aes(x = `shuffle ratio`, y = mean, color = q, linetype=q ,group = q)) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  annotate("text", x = Inf, y = chance_level, label = "chance_level", hjust = 1.1, vjust = -0.5, color = "red") +
  labs(y = 'MSE', x = 'Shuffle ratio', color = 'q') +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )
print(plottt)
plottt <- ggplot(filtered_summary, aes(x = `shuffle ratio`, y = mean, color = q, linetype = q, group = q)) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  annotate("text", x = Inf, y = chance_level, label = "chance level", hjust = 1.1, vjust = -0.5, color = "red") +
  labs(y = 'MSE', x = 'Shuffle ratio', color = 'q', linetype = 'q') +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02)) +  # <-- force y-axis range and ticks
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold")
  )

print(plottt)


plottt <- ggplot(filtered_summary, aes(x = q, y = mean, color = `shuffle ratio`, group = `shuffle ratio`)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_jitter()+
  labs(y = 'MSE', x = 'q', color = 'Shuffle Ratio') +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

print(plottt)

plottt <- ggplot(summary, aes(x = q, y = mean, color = `shuffle ratio`, group = `shuffle ratio`)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # Set `width` for error bars
  labs(y = 'MSE', x = 'q', color = 'Shuffle Ratio') +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold")
  )

print(plottt)
print(plottt)

plottt <- plottt + labs(title = paste("p =", p, "q =", q,"m=",m, "n =", n, "nmc =", nmc, 'max_iter=', 100))
plottt <- plottt + theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25, face = "bold"),
  legend.text = element_text(size = 25),
  legend.title = element_text(size = 25, face = "bold"),
  plot.title = element_text(size = 20, face = "bold") 
)


data_melt <- melt(as.matrix(final_mse))

shuffle_ratio=c(0,d)
data_melt$Var1 <- factor(data_melt$Var1, labels = shuffle_ratio)
data_melt$Var2 <- factor(data_melt$Var2, labels = q)

# Plot the heatmap
heatmap_plot <- ggplot(data_melt, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#FFF5F5", "#FFD6D6", "#FF9999", "#FF4D4D", "#FF0000", "#B20000", "#7F0000"))(100),
    name = "MSE"
  ) +
  labs(x = "q", y = "Shuffle Ratio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12), # Ensure y-axis text is visible
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Save the plot as a high-resolution image
# ggsave("heatmap_high_res.png", plot = heatmap_plot, dpi = 300, width = 8, height = 6)

# Display the plot
print(heatmap_plot)







