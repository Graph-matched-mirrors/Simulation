

#load("/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_n500_m20_p0.4_q0.2_num_state50_max_iter100_20250321_1208.RData")

#result_name = "/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_n800_m20_p0.4_q0.2_num_state50_max_iter100_20250325_1508.RData"
#result_name = "/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_n500_m20_p0.4_q0.2_num_state50_max_iter100_20250321_1208.RData"
#setwd('/Users/tianyichen/Desktop/Research /PhDresearch/London model with GM/Github/Simulation/Tianyi Chen')

#result_name = "out_dd_n500_m20_p0.4_q0.2_num_state50_max_iter100_20250321_1208.RData"
#result_name = "out_dd_n800_m20_p0.4_q0.2_num_state50_max_iter100_20250325_1508.RData"
#result_name='out_dd_Londonn500_m20_p0.4_q0.3_max_iter100_20250406_0337.RData'

#result_name = "out_dd_n50_m20_p0.4_q0.2_nmc50_num_state50_max_iter10_20250907_1651.RData"
result_name = "out_dd_n500_m20_p0.4_q0.2_nmc300_num_state50_max_iter100_20250908_1644.RData"
# Extract simulation parameters from result_name
pattern <- "n(\\d+)_m(\\d+)_p([0-9.]+)_q([0-9.]+)_nmc(\\d+)_num_state(\\d+)_max_iter(\\d+)"
matches <- regmatches(result_name, regexec(pattern, result_name))[[1]]
n         <- as.numeric(matches[2])
m         <- as.numeric(matches[3])
p         <- as.numeric(matches[4])
q         <- as.numeric(matches[5])
nmc       <- as.numeric(matches[6])
num_state <- as.numeric(matches[7])
max_iter  <- as.numeric(matches[8])

setwd("/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/")
load(result_name)

out_dd

res_matrix <- do.call(rbind, lapply(seq_along(out_dd), function(i) {
  row <- unlist(out_dd[[i]])
  names(row) <- c(
    "true1", "true_iso_d1", "true_iso_d4", "true_iso_d8",
    "shuffle1", "shuffle_iso_d1", "shuffle_iso_d4", "shuffle_iso_d8",
    "gm_alltoone1", "gm_alltoone_iso_d1", "gm_alltoone_iso_d4", "gm_alltoone_iso_d8",
    "gm_pairwise1", "gm_pairwise_iso_d1", "gm_pairwise_iso_d4", "gm_pairwise_iso_d8",
    "W1","Avg_degree"
  )
  return(row)
}))

res_matrix
unlist(out_dd[[1]])

# Number of simulations used in the analysis
nmc <- length(out_dd)

# Calculate summary statistics for each column in res_matrix
col_means <- apply(abs(res_matrix)^2, 2, mean)
col_sds   <- apply(abs(res_matrix)^2, 2, sd)

# Organize the summary into a data frame
summary_df <- data.frame(
  metric = names(col_means),
  mean   = col_means,
  sd     = col_sds
)

# Compute the lower and upper bounds of the 95% Confidence Interval
summary_df$lower_bound <- summary_df$mean - 1.96 * summary_df$sd / sqrt(nmc)
summary_df$upper_bound <- summary_df$mean + 1.96 * summary_df$sd / sqrt(nmc)


library(dplyr)
library(ggplot2)

# Create two new variables: 'type' and 'iso' from the metric names
summary_df <- summary_df %>%
  mutate(type = case_when(
           grepl("^true", metric)       ~ "dMV true alignment",
           grepl("^shuffle", metric)      ~ "dMV shuffled",
           grepl("^gm_alltoone", metric)  ~ "GM all to one",
           grepl("^gm_pairwise", metric)  ~ "GM consecutive"
         ),
         iso = case_when(
           grepl("iso_d1", metric) ~ "iso_d1+D",
           grepl("iso_d4", metric) ~ "iso_d4+D",
           grepl("iso_d8", metric) ~ "iso_d8+D",
           TRUE                  ~ "MDS1+D^2"   # Default case: plain "1"
         ))


summary_df
summary_df$iso <- factor(summary_df$iso, levels = c("MDS1+D^2", "iso_d1+D", "iso_d4+D", "iso_d8+D"))

support = seq(2/m-0.5, (m-1)/m-0.5,by=1/m)
chance_level = sum(support^2/length(support))
chance_level
# Plot: x-axis is the iso label, and different colors indicate different types

summary_df_filtered = summary_df[summary_df$iso!='iso_d1+D',]

plot_summary <- ggplot(summary_df_filtered, aes(x = iso, y = mean, group = type, color = type)) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0, 0.18))+
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = chance_level, label = "chance_level", hjust = 1.1, vjust = -0.5, color = "red") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
  labs(x = "", 
       y = "MSE", 
       #title = paste("n=",n, 'm=',m, 'p=',p, 'q=',q, 'num_state=',num_state, 'max_iter=',max_iter),
       color = "metric") +
  theme_minimal() +
  scale_color_manual(values = c(
    "dMV true alignment" = "#7CAE00",  # same as previous dMV
    "dMV shuffled"       = "#D84315",  # a redder shade than orange, but not pure red
    "GM all to one"       = "#00BFC4",  # similar blue hue for both alltoone and pairwise
    "GM consecutive"       = "#00B0F6"   # similar blue hue for both alltoone and pairwise
  )) +
  theme(axis.text = element_text(size = 25),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        #legend.position = "none",  # This line removes the legend
        legend.title = element_text(size = 20))

print(plot_summary)
plot_title <- paste('Atlanta-l2-localizer', 'n =', n, 'nmc =', nmc, 'p =', p, 'q =', q)
filename <- paste0(plot_title, ".pdf")
filename
ggsave(filename = filename, plot = plot_summary, device = "pdf", width = 10, height = 6)







######
mean(res_matrix[, 'gm_pairwise_iso_d1' ]^2)
mean(res_matrix[, 'gm_pairwise_iso_d4' ]^2)
mean(res_matrix[, 'gm_pairwise_iso_d8' ]^2)


x_ticks <- seq(-0.5, 0.5, by = 1/m)

hist(abs(res_matrix[, 'gm_pairwise_iso_d1']),
     breaks = 50,
     xlim = c(0, 0.5),
     xaxt = "n",
     main = "Histogram of gm_pairwise_iso_d1",
     xlab = "Value")
axis(1, at = x_ticks, labels = x_ticks)

hist(abs(res_matrix[, 'gm_pairwise_iso_d4']),
     breaks = 50,
     xlim = c(0, 0.5),
     xaxt = "n",
     main = "Histogram of gm_pairwise_iso_d4",
     xlab = "Value")
axis(1, at = x_ticks, labels = x_ticks)


hist((res_matrix[, 'gm_pairwise_iso_d4']),
     breaks = 50,
     xlim = c(-0.5, 0.5),
     xaxt = "n",
     main = "Histogram of gm_pairwise_iso_d4",
     xlab = "Value")
axis(1, at = x_ticks, labels = x_ticks)

hist(abs(res_matrix[, 'shuffle1']),
     breaks = 50,
     xlim = c(0, 0.5),
     xaxt = "n",
     main = "Histogram of gm_pairwise_iso_d8",
     xlab = "Value")
axis(1, at = x_ticks, labels = x_ticks)


