# ==================== SETUP ====================
library(dbscan)
library(ggplot2)
library(RColorBrewer)

# ----------------------- User parameters -----------------------
samples <- c("CFPM1D7KC_points_manual_detection", "CFPM1D0KC_points_manual_detection")
legend_labels <- c("D0", "D8")  # Customize here
minPts <- 2
max_points_per_cluster <- 20

eps_values <- seq(10, 200, by = 10)   # For DBSCAN
t1_values <- seq(10, 200, by = 10)    # For Canopy
k_values <- 2:20                      # For K-Means

setwd("C:/Users/maesr/Desktop/Ugent/2nd master biochemistry biotechnology/Master II/UBow/UBow_results/2025_04_22_RobinM_Axioscan/AnalysisKC_M1")
# ---------------------------------------------------------------


# ==================== FUNCTIONS ====================

# Prism-style plot theme
prism_theme <- function() {
  theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2),
      axis.text = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "right"
    )
}

# Plot and save PNG with average points text labels
plot_and_save <- function(df, xvar, yvar, colorvar, outfile, xlab, ylab) {
  p <- ggplot(df, aes_string(x = xvar, y = yvar, color = colorvar, group = colorvar)) +
    geom_line(size = 1.8) +
    geom_point(size = 3) +
    geom_text(aes(label = round(avg_points, 2)), vjust = -1.2, size = 3, show.legend = FALSE) +  # labels above points
    prism_theme() +
    scale_color_brewer(palette = "Set1") +
    labs(title = NULL, x = xlab, y = ylab) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1))))))
  
  ggsave(paste0(outfile, ".svg"), plot = p, width = 10, height = 6, dpi = 300)
  cat("Saved:", outfile, "\n")
}

# Filter clusters by size and minPts
process_clusters <- function(cluster_labels, coords, max_points_per_cluster, minPts = 1) {
  cluster_labels <- as.factor(cluster_labels)
  clustered <- cluster_labels[cluster_labels != 0]
  cluster_sizes <- table(clustered)
  valid_clusters <- cluster_sizes[cluster_sizes <= max_points_per_cluster & cluster_sizes >= minPts]
  n_clusters <- length(valid_clusters)
  avg_points <- if (n_clusters > 0) mean(valid_clusters) else 0
  return(list(n_clusters = n_clusters, avg_points = avg_points))
}

# Simulated Canopy Clustering
canopy_cluster <- function(coords, t1, t2) {
  points_left <- 1:nrow(coords)
  clusters <- rep(0, nrow(coords))
  cluster_id <- 1
  
  while (length(points_left) > 0) {
    center_idx <- sample(points_left, 1)
    center <- coords[center_idx, , drop = FALSE]
    dists <- sqrt(rowSums((t(t(coords) - as.numeric(center)))^2))
    in_canopy <- which(dists <= t1)
    in_core <- which(dists <= t2)
    
    if (length(in_canopy) > 0) {
      clusters[in_canopy] <- cluster_id
      cluster_id <- cluster_id + 1
    }
    points_left <- setdiff(points_left, in_core)
  }
  
  return(clusters)
}

# ==================== DBSCAN ====================
results_dbscan <- data.frame()
for (sample in samples) {
  df <- read.csv(paste0(sample, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  coords <- df[, c("Centroid.X", "Centroid.Y")]
  
  for (eps in eps_values) {
    db <- dbscan(coords, eps = eps, minPts = minPts, borderPoints = TRUE)
    stats <- process_clusters(db$cluster, coords, max_points_per_cluster, minPts)
    results_dbscan <- rbind(results_dbscan, data.frame(
      sample = sample,
      param = eps,
      num_clusters = stats$n_clusters,
      avg_points = stats$avg_points
    ))
  }
}
results_dbscan$sample_label <- factor(results_dbscan$sample, labels = legend_labels)
plot_and_save(results_dbscan, "param", "num_clusters", "sample_label",
              "dbscan_cluster_count_vs_eps_prism_style", "Epsilon (eps)", "Number of Clusters")

# ==================== CANOPY ====================
results_canopy <- data.frame()
for (sample in samples) {
  df <- read.csv(paste0(sample, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  coords <- df[, c("Centroid.X", "Centroid.Y")]
  
  for (t1 in t1_values) {
    t2 <- t1 / 2
    clusters <- canopy_cluster(coords, t1, t2)
    stats <- process_clusters(clusters, coords, max_points_per_cluster, minPts)
    results_canopy <- rbind(results_canopy, data.frame(
      sample = sample,
      param = t1,
      num_clusters = stats$n_clusters,
      avg_points = stats$avg_points
    ))
  }
}
results_canopy$sample_label <- factor(results_canopy$sample, labels = legend_labels)
plot_and_save(results_canopy, "param", "num_clusters", "sample_label",
              "canopy_cluster_count_vs_t1_prism_style", "Threshold t1", "Number of Clusters")
