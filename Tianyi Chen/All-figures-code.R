summary_df_both <- read.csv("London-both-localizer n = 200 nmc = 500 p = 0.4 q = 0.3.csv")

plottt1 <- ggplot(summary_df_both, aes(x = m, y = Mean_MSE, color = metric, linetype = localizer)) +
  geom_point(size = 1 )+
  geom_line() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
  scale_x_continuous(breaks = mm) +
  labs(y = 'MSE', x = 'm') +
  scale_color_manual(values = c(
    "avg_degree" = "#F8766D",  # pinkish
    "dMV"        = "#7CAE00",  # green
    "W1"         = "#C77CFF",  # purple
    "W2"         = "#00BFC4"   # teal/blue
  )) +
  theme(
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )
print(plottt1)
plot_title <- paste('London-both-localizer', 'n =', n, 'nmc =', nmc, 'p = 0.4', 'q = 0.3')
filename <- paste0(plot_title, ".pdf")
ggsave(filename = filename, plot = plottt1, device = "pdf", width = 8, height = 6)


#out_dd <- read.csv("/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_20250320_131154.csv", header = TRUE, stringsAsFactors = FALSE)
library(ggplot2)



#result_name ="/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData"
#result_name ="/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.2_max_iter100_20250328_0005.RData"
#result_name = '/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.3_q0.4_max_iter100_20250328_0114.RData'
#result_name = '/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData'


## for local load 
setwd('/Users/tianyichen/Desktop/Research /PhDresearch/London model with GM/Github/Simulation/Tianyi Chen')
#result_name = 'out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData'
result_name = 'out_dd_Londonn200_m20_p0.3_q0.4_max_iter100_20250328_0114.RData'

# Extract simulation parameters from result_name
pattern <- "n(\\d+)_m(\\d+)_p([0-9.]+)_q([0-9.]+)_max_iter(\\d+)"
matches <- regmatches(result_name, regexec(pattern, result_name))[[1]]
matches

n         <- as.numeric(matches[2])
m         <- as.numeric(matches[3])
p         <- as.numeric(matches[4])
q         <- as.numeric(matches[5])
max_iter  <- as.numeric(matches[6])



load(result_name)



#out_dd

# Create a summary data frame from out_dd (ignoring the example_df)
summary_df <- do.call(rbind, lapply(seq_along(out_dd), function(i) {
  out <- out_dd[[i]]
  data.frame(
    mc = i,
    true_mds1_iso1 = out$tmp1[1],
    true_mds2_iso1 = out$tmp1[2],
    true_mds3_iso1 = out$tmp1[3],
    shuffled_mds1_iso1 = out$tmp2[1],
    shuffled_mds2_iso1 = out$tmp2[2],
    shuffled_mds3_iso1 = out$tmp2[3],
    gm_alltoone_mds1_iso1 = out$tmp3[1],
    gm_alltoone_mds2_iso1 = out$tmp3[2],
    gm_alltoone_mds3_iso1 = out$tmp3[3],
    gm_pairwise_mds1_iso1 = out$tmp4[1],
    gm_pairwise_mds2_iso1 = out$tmp4[2],
    gm_pairwise_mds3_iso1 = out$tmp4[3],
    tmp_W1 = out$tmp_W1,
    tmp_avg_edges = out$tmp_avg_edges
  )
}))



num_df <- summary_df[, !names(summary_df) %in% "mc"]

squared_df <- num_df^2

nmc <- nrow(summary_df)

# Compute the mean and standard deviation for each metric (each column)
col_means <- colMeans(squared_df)
col_sds   <- apply(squared_df, 2, sd)

lower_bound <- col_means - 1.96 * col_sds / sqrt(nmc)
upper_bound <- col_means + 1.96 * col_sds / sqrt(nmc)

ci_df <- data.frame(lower_bound = lower_bound, mean = col_means, upper_bound = upper_bound)
ci_df$metric <- rownames(ci_df)


ci_df$iso <- ifelse(grepl("mds1", ci_df$metric), "iso1",
                    ifelse(grepl("mds2", ci_df$metric), "iso2",
                           ifelse(grepl("mds3", ci_df$metric), "iso3", NA)))

ci_df$category <- ifelse(grepl("true", ci_df$metric), "dMV true alignment",
                         ifelse(grepl("shuffled", ci_df$metric), "dMV shuffled",
                                ifelse(grepl("alltoone", ci_df$metric), "GM all to one",
                                       ifelse(grepl("pairwise", ci_df$metric), "GM concecutive", NA))))

ci_df$iso[grepl("tmp_W1|tmp_avg_edges", ci_df$metric)] <- "iso1"
ci_df$category[grepl("tmp_W1", ci_df$metric)] <- "W1"
ci_df$category[grepl("tmp_avg_edges", ci_df$metric)] <- "avg degree"

ci_df$category
support = seq(2/m-0.5, (m-1)/m-0.5,by=1/m)
chance_level = sum(support^2/length(support))


ci_df_iso1 <- ci_df[ci_df$iso == "iso1", ]

# Create the plot for iso1
plottt1 <- ggplot(ci_df_iso1, aes(x = iso, y = mean, color = category)) +
  geom_point(size = 4) +
  scale_y_continuous(limits = c(0, 0.04))+
  #geom_hline(yintercept = chance_level, linetype = "dotted", color = "red", linewidth = 1) +
  #annotate("text", x = Inf, y = chance_level, label = "chance_level", hjust = 1.1, vjust = -0.5, color = "red") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
  labs(
    x = "", 
    y = "MSE", 
    #title = paste("n=", n, 'm=', m, 'p=', p, 'q=', q, 'max_iter=', max_iter, " (iso1 only)"),
    color = "metric"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "dMV true alignment" = "#7CAE00",  # same as previous dMV
    "W1"             = "#C77CFF",  # same as before
    "avg degree"      = "#F8766D",  # same as previous avg_degree
    "dMV shuffled"       = "#D84315",  # a redder shade than orange, but not pure red
    "GM all to one"       = "#00BFC4",  # similar blue hue for both alltoone and pairwise
    "GM concecutive"       = "#00B0F6"   # similar blue hue for both alltoone and pairwise
  )) +
  theme(
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    #legend.position = "none",  # This line removes the legend
    plot.title = element_text(size = 18, face = "bold")
  )
print(plottt1)
#
plot_title <- paste('London-l2-localizer-MDS1', 'n =', n, 'nmc =', nmc, 'p =', p, 'q =', q)
filename <- paste0(plot_title, ".pdf")
filename
ggsave(filename = filename, plot = plottt1, device = "pdf", width = 8, height = 6)






p1 <- ggplot(ci_df, aes(x = iso, y = mean, group = category, color = category)) +
  geom_point(size = 1) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = chance_level, label = "chance level", hjust = 1.1, vjust = -0.5, color = "red") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
  labs(x = "MDS d ISOMAP to 1", 
       y =  "MSE", 
       color = "Category") +
  theme_minimal() +
  scale_color_manual(values = c(
    "dMV true alignment" = "#7CAE00",  # same as previous dMV
    "W1"             = "#C77CFF",  # same as before
    "avg degree"      = "#F8766D",  # same as previous avg_degree
    "dMV shuffled"       = "#D84315",  # a redder shade than orange, but not pure red
    "GM all to one"       = "#00BFC4",  # similar blue hue for both alltoone and pairwise
    "GM consecutive"       = "#00B0F6"   # similar blue hue for both alltoone and pairwise
  )) +
  theme(axis.text = element_text(size = 25),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

print(p1)
#plot_title <- paste('London-l2-localizer', 'n =', n, 'nmc =', nmc, 'p =', p, 'q =', q)
#filename <- paste0(plot_title, ".pdf")
#ggsave(filename = filename, plot = plottt1, device = "pdf", width = 8, height = 6)


# Filter the data for only iso1



