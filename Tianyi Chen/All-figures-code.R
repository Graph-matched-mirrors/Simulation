summary_df_both <- read.csv("London-both-localizer n = 200 nmc = 500 p = 0.4 q = 0.3.csv")
plottt1 <- ggplot(summary_df_both, aes(x = m, y = Mean_MSE, color = metric, linetype = localizer)) +
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm') +
  theme(
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
print(plottt1)
plot_title <- paste('London-both-localizer', 'n =', n, 'nmc =', nmc, 'p =', p, 'q =', q)
filename <- paste0(plot_title, ".pdf")
ggsave(filename = filename, plot = plottt1, device = "pdf", width = 8, height = 6)



