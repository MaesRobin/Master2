# Load libraries (skip if already loaded)
library(tidyverse)
library(readr)
library(spatstat.geom)
library(spatstat)
library(spatstat.explore)

# Set working directory
setwd("C:/Users/maesr/Desktop/Ugent/2nd master biochemistry biotechnology/Master II/UBow/UBow_results/2025_04_22_RobinM_Axioscan/AnalysisKC_M1")

# Load both datasets
file1 <- "CFPM1D0KC_points_manual_detection.txt"
file2 <- "CFPM1D7KC_points_manual_detection.txt"
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

# Function to read and convert a dataset to point pattern
load_point_pattern <- function(file_path, coord_cols = c("Centroid X", "Centroid Y")) {
  data <- read_tsv(file_path, col_types = cols(.default = "c")) %>%
    mutate(across(all_of(coord_cols), ~ as.numeric(.x))) %>%
    drop_na(all_of(coord_cols))
  
  x <- data[[coord_cols[1]]]
  y <- data[[coord_cols[2]]]
  win <- owin(range(x, na.rm = TRUE), range(y, na.rm = TRUE))
  ppp(x, y, window = win)
}

# Create point pattern objects
pp1 <- load_point_pattern(file1)
pp2 <- load_point_pattern(file2)

# Compute K envelopes with 999 simulations
env1 <- envelope(pp1, Kest, nsim = 999, correction = "Ripley", savefuns = TRUE, verbose = FALSE)
env2 <- envelope(pp2, Kest, nsim = 999, correction = "Ripley", savefuns = TRUE, verbose = FALSE)

# Save individual Ripley plots with envelopes
output_filename_env1 <- paste0("ripley_env_", tools::file_path_sans_ext(file1), ".png")
png(filename = output_filename_env1, width = 800, height = 600)
plot(env1, main = paste("Ripley's K with CSR Envelope for", basename(file1)))
dev.off()
cat("Ripley envelope plot for", basename(file1), "saved as:", output_filename_env1, "\n")

output_filename_env2 <- paste0("ripley_env_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_env2, width = 800, height = 600)
plot(env2, main = paste("Ripley's K with CSR Envelope for", basename(file2)))
dev.off()
cat("Ripley envelope plot for", basename(file2), "saved as:", output_filename_env2, "\n")

# Compute observed K functions
K1 <- Kest(pp1, correction = "Ripley")
K2 <- Kest(pp2, correction = "Ripley")

# Compute theoretical CSR line: K(r) = πr²
r_vals <- K1$r  # assumes both use same r vector
K_CSR <- pi * r_vals^2

# Compute difference: ΔK(r) = K_obs(r) - πr²
dK1 <- K1$iso - env1$theo
dK2 <- K2$iso - env2$theo

# Plot difference curves
output_filename_diff <- paste0("ripley_difference_", tools::file_path_sans_ext(file1), "_vs_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_diff, width = 800, height = 600)
plot(r_vals, dK1, type = "l", col = "black", ylim = range(c(dK1, dK2), na.rm = TRUE),
     xlab = "r", ylab = expression(Delta*K(r)),
     main = "Difference from CSR: ΔK(r) = Kobs(r) - πr²")
lines(r_vals, dK2, col = "blue")
abline(h = 0, col = "red", lty = 2)  # CSR baseline
legend("topright", legend = c(basename(file1), basename(file2), "CSR (ΔK = 0)"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2))
dev.off()

cat("ΔK(r) difference plot saved as:", output_filename_diff, "\n")

# Compute difference of the differences: ΔΔK(r) = ΔK1(r) - ΔK2(r) = K1 - K2
dK_diff <- dK1 - dK2  # or equivalently: K1$iso - K2$iso

# Plot ΔΔK(r)
output_filename_ddiff <- paste0("ripley_double_difference_", tools::file_path_sans_ext(file1), "_vs_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_ddiff, width = 800, height = 600)
plot(r_vals, dK_diff, type = "l", col = "purple",
     xlab = "r", ylab = expression(Delta*Delta*K(r)),
     main = expression("Difference of ΔK(r): ΔΔK(r) = ΔK[1](r) - ΔK[2](r)"))
abline(h = 0, col = "red", lty = 2)
legend("topright", legend = c("ΔΔK(r)", "Zero Line"),
       col = c("purple", "red"), lty = c(1, 2))
dev.off()

cat("ΔΔK(r) difference-of-differences plot saved as:", output_filename_ddiff, "\n")

# Compute signed log2 fold change: log2(|dK1/dK2|) * sign(dK1 - dK2)
epsilon <- 0  # small number to prevent division by zero
fold_change <- log2((abs(dK2) + epsilon)/(abs(dK1) + epsilon))

# Plot signed log2 fold change
output_filename_fc <- paste0("ripley_fold_change_", tools::file_path_sans_ext(file1), "_vs_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_fc, width = 800, height = 600)
plot(r_vals, fold_change, type = "l", col = "darkgreen",
     xlab = "r", ylab = expression("Signed log"[2]*" Fold Change of ΔK(r)"),
     main = expression("Fold Change: log"[2]*"(ΔK[2]/ΔK[1])"))
abline(h = 0, col = "red", lty = 2)
legend("topright", legend = c("Fold Change", "No Change"),
       col = c("darkgreen", "red"), lty = c(1, 2))
dev.off()

cat("Signed fold change plot saved as:", output_filename_fc, "\n")

# --- Intensity and Neighbors ---
# Calculate intensity (points per µm²)
area1 <- area.owin(pp1$window)
area2 <- area.owin(pp2$window)
lambda1 <- pp1$n / area1
lambda2 <- pp2$n / area2

# Observed neighbors = lambda * Kobs
obs_neighbors1 <- lambda1 * K1$iso
obs_neighbors2 <- lambda2 * K2$iso

# Expected neighbors under CSR = lambda * π * r²
expected_neighbors1 <- lambda1 * (pi * r_vals^2)
expected_neighbors2 <- lambda2 * (pi * r_vals^2)

# --- Plot Observed vs Expected Neighbors ---
output_filename_neighbors <- paste0("ripley_neighbors_", tools::file_path_sans_ext(file1), "_vs_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_neighbors, width = 800, height = 600)
plot(r_vals, expected_neighbors1, type = "l", lty = 2, col = "red",
     xlab = "r (µm)", ylab = "Neighbors per Point",
     main = "Observed vs Expected Neighbor Counts", ylim = range(c(expected_neighbors1, expected_neighbors2, obs_neighbors1, obs_neighbors2)))
lines(r_vals, obs_neighbors1, col = "black", lwd = 2)
lines(r_vals, expected_neighbors2, col = "orange", lty = 2)
lines(r_vals, obs_neighbors2, col = "blue", lwd = 2)
legend("topleft", legend = c("Expected - File1", "Observed - File1", "Expected - File2", "Observed - File2"),
       col = c("red", "black", "orange", "blue"), lty = c(2,1,2,1), lwd = c(1,2,1,2))
dev.off()

cat("Observed vs Expected neighbors plot saved as:", output_filename_neighbors, "\n")

# --- L(r) - r Plot for Clustering/Dispersion ---
L1 <- Lest(pp1, correction = "Ripley")
L2 <- Lest(pp2, correction = "Ripley")

output_filename_Lr <- paste0("ripley_Lr_minus_r_", tools::file_path_sans_ext(file1), "_vs_", tools::file_path_sans_ext(file2), ".png")
png(filename = output_filename_Lr, width = 800, height = 600)
# Compute L(r) - r values
L1_minus_r <- L1$iso - L1$r
L2_minus_r <- L2$iso - L2$r
y_range <- range(c(L1_minus_r, L2_minus_r), na.rm = TRUE)

# Plot with corrected y-axis limits
plot(L1$r, L1_minus_r, type = "l", col = "black", lwd = 2,
     xlab = "r (µm)", ylab = expression(L(r) - r),
     main = expression("L(r) - r Plot (Clustering if > 0)"),
     ylim = y_range)
lines(L2$r, L2_minus_r, col = "blue", lwd = 2)
abline(h = 0, col = "red", lty = 2)
legend("topright", legend = c(basename(file1), basename(file2), "CSR baseline"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = c(2, 2, 1))

dev.off()

cat("L(r) - r clustering plot saved as:", output_filename_Lr, "\n")


# Calculate p-values for Day0 and Day7 samples
p_values_day0 <- sapply(r_vals, function(r) {
  mean(env1$simulations[ ,which.min(abs(env1$r - r))] > K1$iso[which.min(abs(K1$r - r))])
})

p_values_day7 <- sapply(r_vals, function(r) {
  mean(env2$simulations[ ,which.min(abs(env2$r - r))] > K2$iso[which.min(abs(K2$r - r))])
})

# Perform a Kolmogorov-Smirnov test to compare the K-functions between Day0 and Day7
ks_test <- ks.test(K1$iso, K2$iso)
cat("p-value from K-S test:", ks_test$p.value, "\n")
# Use Diggle's test to compare the two point patterns (Day0 vs Day7)




