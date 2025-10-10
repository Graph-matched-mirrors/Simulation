## This is the file used to summarized the results generated from MCSim_London.R


#out_dd <- read.csv("/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_20250320_131154.csv", header = TRUE, stringsAsFactors = FALSE)
library(ggplot2)



#result_name ="/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData"
#result_name ="/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.2_max_iter100_20250328_0005.RData"
#result_name = '/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.3_q0.4_max_iter100_20250328_0114.RData'
#result_name = '/cis/home/tchen94/tianyi/Simulation/Tianyi Chen/out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData'


## for local load 
setwd('/Users/tianyichen/Desktop/Research /PhDresearch/London model with GM/Github/Simulation/Tianyi Chen')
result_name = 'out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData'

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


ci_df$category <- ifelse(grepl("true", ci_df$metric), "true alignment",
             ifelse(grepl("shuffled", ci_df$metric), "shuffled",
                ifelse(grepl("alltoone", ci_df$metric), "alltoone",
                     ifelse(grepl("pairwise", ci_df$metric), "pairwise", NA))))

ci_df$iso[grepl("tmp_W1|tmp_avg_edges", ci_df$metric)] <- "iso1"

ci_df$category[grepl("tmp_W1", ci_df$metric)] <- "W1"
ci_df$category[grepl("tmp_avg_edges", ci_df$metric)] <- "avg_edges"


support = seq(2/m-0.5, (m-1)/m-0.5,by=1/m)
chance_level = sum(support^2/length(support))


library(dplyr)

filtered_df <- ci_df %>%
  group_by(category) %>%
  slice_min(mean, with_ties = FALSE) %>%
  ungroup()

p2 <- ggplot(filtered_df, aes(x = iso, y = mean, group = category, color = category)) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 0.003))+
  #geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  #annotate("text", x = Inf, y = chance_level, label = "chance level", hjust = 1.1, vjust = -0.5, color = "red") +
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
  ))+
  theme(axis.text = element_text(size = 25),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

print(p2)

p1 <- ggplot(ci_df, aes(x = iso, y = mean, group = category, color = category)) +
  geom_point(size = 2) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "red") +
  annotate("text", x = Inf, y = chance_level, label = "chance level", hjust = 1.1, vjust = -0.5, color = "red") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
  labs(x = "MDS d ISOMAP to 1", 
       y =  "MSE", 
       color = "Category") +
  theme_minimal() +
  theme(axis.text = element_text(size = 25),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

print(p1)
plot_title <- paste('London-both-localizer', 'n =', n, 'nmc =', nmc, 'p =', p, 'q =', q)
filename <- paste0(plot_title, ".pdf")
ggsave(filename = filename, plot = plottt1, device = "pdf", width = 8, height = 6)


        # Filter the data for only iso1








ci_df_iso1 <- ci_df[ci_df$iso == "iso1", ]
ci_df_iso2 <- ci_df[ci_df$iso == "iso2", ]


        # Create the plot for iso1
        ggplot(ci_df_iso1, aes(x = iso, y = mean, group = category, color = category)) +
          geom_point(size = 4) +
          geom_line(linetype = "dashed", linewidth = 1) +
          #geom_hline(yintercept = chance_level, linetype = "dotted", color = "red", linewidth = 1) +
          #annotate("text", x = Inf, y = chance_level, label = "chance_level", hjust = 1.1, vjust = -0.5, color = "red") +
          geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
          labs(
            x = "Category", 
            y = "MSE", 
            title = paste("n=", n, 'm=', m, 'p=', p, 'q=', q, 'max_iter=', max_iter, " (iso1 only)"),
            color = "Category"
          ) +
          theme_minimal() +
          theme(
            axis.text = element_text(size = 14),
            axis.title.x = element_text(size = 18, face = "bold"),
            axis.title.y = element_text(size = 18, face = "bold"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            plot.title = element_text(size = 18, face = "bold")
          )
#




df = out_dd[[300]]$example_df

tstar = m/2
D2 <- getD(df$xhat)
sqrt_edges <- sqrt(unlist(df$avg_edges))
Dhat_W1 <- getD_W1(df$Xhat_shuffle)
Dhat_W1_2 <- getD_W1(df$xhat) 
Dhat_W1_3 <- getD_W1(df$Xhat_shuffle_GM_alltoone)
Dhat_W1_4 <- getD_W1(df$Xhat_shuffle_GM_pairwise)

#norm(Dhat_W1_3 - Dhat_W1_4,'F')
#(Dhat_W1_3 - Dhat_W1_4)/Dhat_W1_4

tmp_avg_edges <- find_slope_changepoint_with_plot(sqrt_edges, doplot = T)$error
  
D2_shuffle <- getD(df$Xhat_shuffle)
D2_shuffle_GM_alltoone <- getD(df$Xhat_shuffle_GM_alltoone)
D2_shuffle_GM_pairwise <- getD(df$Xhat_shuffle_GM_pairwise)
  

df.mds_W1 <- doMDS(Dhat_W1, doplot = T)
df.mds_no_square <- doMDS(D2, doplot = T)
df.mds_shuffle_no_square <- doMDS(D2_shuffle, doplot = T)
df.mds_shuffle_GM_alltoone_no_square <- doMDS(D2_shuffle_GM_alltoone, doplot = T)
df.mds_shuffle_GM_pairwise_no_square <- doMDS(D2_shuffle_GM_pairwise, doplot = T)

find_slope_changepoint_with_plot(sqrt_edges, doplot = T)$error
find_slope_changepoint_with_plot(df.mds_W1$mds[,1], doplot = T)$error
find_slope_changepoint_with_plot(df.mds_no_square$mds[,1], doplot = T)$error
find_slope_changepoint_with_plot(df.mds_shuffle_no_square$mds[,1], doplot = T)$error
find_slope_changepoint_with_plot(df.mds_shuffle_GM_alltoone_no_square $mds[,1], doplot = T)$error
find_slope_changepoint_with_plot(df.mds_shuffle_GM_pairwise_no_square$mds[,1], doplot = T)$error



d_chose = c(1, 2, 3)


k=3
dd = d_chose[k]
dd= 8
    
df.iso <- doIso(df.mds_no_square$mds, mdsd = dd)$iso
df.iso_shuffle <- doIso(df.mds_shuffle_no_square$mds, mdsd = dd)$iso
df.iso_shuffle_GM_alltoone <- doIso(df.mds_shuffle_GM_alltoone_no_square$mds, mdsd = dd)$iso
df.iso_shuffle_GM_pairwise <- doIso(df.mds_shuffle_GM_pairwise_no_square$mds, mdsd = dd)$iso
    
tmp1[k] <- find_slope_changepoint_with_plot(df.iso, doplot = T)$error
tmp2[k] <- find_slope_changepoint_with_plot(df.iso_shuffle, doplot = T)$error
tmp3[k] <- find_slope_changepoint_with_plot(df.iso_shuffle_GM_alltoone, doplot = T)$error
tmp4[k] <- find_slope_changepoint_with_plot(df.iso_shuffle_GM_pairwise, doplot = T)$error


find_slope_changepoint_with_plot(df.iso, doplot = T)$error


OFV_true = OFV_shuffle = OFV_alltoone = OFV_pairwise = matrix(0,m,m)



for (i in 1:m) {
  for (j in 1:m) {
    OFV_true[i,j] =  norm(df$g[[i]][] -  df$g[[j]][],'F')
    OFV_shuffle[i,j] = norm(df$shuffle_g[[i]][] -  df$shuffle_g[[j]][],'F')
    OFV_alltoone[i,j]  =  norm(df$shuffle_g_GM_alltoone[[i]][] -  df$shuffle_g_GM_alltoone[[j]][],'F')
    OFV_pairwise[i,j]  =  norm(df$shuffle_g_GM_pairwise[[i]][] -  df$shuffle_g_GM_pairwise[[j]][],'F')
  }
}

library(pheatmap)



pheatmap(OFV_true,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,    # Show the value in each cell
         color = heat.colors(256),
         main = "Heatmap of OFV with true alignment")

pheatmap(OFV_shuffle,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,    # Show the value in each cell
         color = heat.colors(256),
         main = "Heatmap of OFV before GM with Values")

MDS_OFV=doMDS(OFV)
plot(1:20,sqrt(MDS_OFV$mds[,1]+30))

pheatmap(OFV_alltoone,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,    # Show the value in each cell
         color = heat.colors(256),
         main = "Heatmap of OFV_alltoone with Values")

doMDS(OFV_alltoone)



pheatmap(OFV_pairwise,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,    # Show the value in each cell
         color = heat.colors(256),
         main = "Heatmap of OFV_pairwise with Values")

doMDS(OFV_pairwise)
