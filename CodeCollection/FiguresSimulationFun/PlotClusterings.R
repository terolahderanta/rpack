
# 21.04.2020

# Load packages:

library(rpack)
library(gridExtra)
library(ggplot2)

source("../plot_clusters.R")

# Load data:

test_dat = read.table("../ResultsSimulationFun/test_data.txt", header=T)

clust_centers = read.table("../ResultsSimulationFun/Cluster_centers.txt", header = T)

k = 10

# Permute clusters to same order to get nicer colors:

#old = c(1:10, 99)

# Change labels of figures:

#new = c(2, 1, 3, 8, 6, 7, 5, 9, 4, 10, 99) # corresponds to "c_col"

#test_dat$SqrE_limits_continuous = as.integer(as.character(
#  factor(test_dat$SqrE_limits_continuous, old, new)))

# Overall, there are four (4) different plots.

c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkred","yellow3","darkgreen","bisque4","magenta",
          "royalblue","tomato4","steelblue1",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")

	
##########################################################################################################

# Squared Euclidean, lower and upper limit plus penalty

p1 = plot_clusters(
  coords = test_dat[,1:2],
  weights = test_dat$w,
  clusters = test_dat$SqrE_limits_continuous,
  centers = cbind(clust_centers$SqrdEuc_continuous_x, clust_centers$SqrdEuc_continuous_y),
  subtitle = "(A)",
  outgroup_legend = "outgrp."
)

##########################################################################################################

# Squared Euclidean, NO lower and upper limit plus penalty

p2 = plot_clusters(
  coords = test_dat[,1:2],
  weights = test_dat$w,
  clusters = test_dat$SqrE_nolimits_continuous,
  centers = cbind(clust_centers$SqrdEuc_continuous_no_limit_x, clust_centers$SqrdEuc_continuous_no_limit_y),
  subtitle = "(B)",
  outgroup_legend = "outgrp."
)

##########################################################################################################

gA = ggplotGrob(p1)
gB = ggplotGrob(p2)

# Set the widths

gA$widths = gB$widths

# Arrange the two charts.
# The legend boxes are centered

png("Artificial_squared_Euc.png", width=1000, height=500,)

grid.arrange(gA, gB, ncol = 2)

dev.off()

##########################################################################################################

# Euclidean, lower and upper limit, NO penalty

p3 = plot_clusters(
  coords = test_dat[,1:2],
  weights = test_dat$w,
  clusters = test_dat$E_no_penalty_continuous,
  centers = cbind(clust_centers$Euc_continuous_no_penalty_x, clust_centers$Euc_continuous_no_penalty_y),
  subtitle = "(B)",
  outgroup_legend = "outgrp."
)

##########################################################################################################

# Euclidean, lower and upper limit plus penalty

p4 = plot_clusters(
  coords = test_dat[,1:2],
  weights = test_dat$w,
  clusters = test_dat$E_limits_continuous,
  centers = cbind(clust_centers$Euc_continuous_x, clust_centers$Euc_continuous_y),
  subtitle = "(A)",
  outgroup_legend = "outgrp."
)

##########################################################################################################

gA = ggplotGrob(p4)
gB = ggplotGrob(p3)

# Set the widths

gA$widths = gB$widths

# Arrange the two charts.
# The legend boxes are centered

png("Artificial_Euc.png", width=1000, height=500,)

grid.arrange(gA, gB, ncol = 2)

dev.off()

##########################################################################################################

# Finally, plot the ground truth clusters:

test_dat$orig_group[is.na(test_dat$orig_group)] = 11

PTrue = plot_clusters(
  coords = test_dat[,1:2],
  centers = cbind(rep(0, 10), rep(0, 10)),
  weights = test_dat$w,
  clusters = test_dat$orig_group,
  subtitle = "Ground truth clusters",
  outgroup_label = 11,
  outgroup_legend = "outgrp."
)

PTrue$layers

PTrue$layers[[2]] = NULL # To remove previously defined cluster centers which are nonsense.

PTrue

ggsave("Artificial_OriginalClusters.png", device = 'png', dpi = 600, width = 5.5, height = 5.5)
