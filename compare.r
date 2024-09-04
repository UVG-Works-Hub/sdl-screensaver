# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Read the CSV files
parallel_data <- read.csv("./analysis/parallel_test_results.csv")
sequential_data <- read.csv("./analysis/sequential_test_results.csv")

# Calculate speedup and efficiency
combined_data <- parallel_data %>%
  inner_join(sequential_data, by = c("Zoom", "MaxIterations", "ScreenWidth", "ScreenHeight"), suffix = c("_parallel", "_sequential")) %>%
  mutate(
    Speedup = AverageFPS_parallel / AverageFPS_sequential,
    Efficiency = Speedup / Threads,
    RelativePerformance = AverageRenderTime_sequential / AverageRenderTime_parallel
  )


# New function to create a bubble plot
create_bubble_plot <- function(data, x, y, size, fill, title) {
  ggplot(data, aes(x = !!sym(x), y = !!sym(y), size = !!sym(size), fill = !!sym(fill))) +
    geom_point(alpha = 0.7, shape = 21, color = "black") +
    scale_size_continuous(range = c(1, 15)) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = title, x = x, y = y) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_x_log10(breaks = unique(data[[x]]), labels = unique(data[[x]])) +
    scale_y_log10(breaks = unique(data[[y]]), labels = unique(data[[y]]))
}


speedup_bubble_plot <- create_bubble_plot(combined_data, "Zoom", "MaxIterations", "Speedup", "Speedup", "Speedup Bubble Plot")
ggsave("./analysis/speedup_bubble_plot.png", speedup_bubble_plot, width = 12, height = 10, dpi = 300)

efficiency_bubble_plot <- create_bubble_plot(combined_data, "Zoom", "MaxIterations", "Efficiency", "Efficiency", "Efficiency Bubble Plot")
ggsave("./analysis/efficiency_bubble_plot.png", efficiency_bubble_plot, width = 12, height = 10, dpi = 300)


# Create line plots for FPS comparison
fps_comparison <- ggplot(combined_data, aes(x = MaxIterations)) +
  geom_line(aes(y = AverageFPS_parallel, color = "Parallel")) +
  geom_line(aes(y = AverageFPS_sequential, color = "Sequential")) +
  facet_wrap(~ Zoom, scales = "free_y") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "FPS Comparison: Parallel vs Sequential", x = "Max Iterations", y = "Average FPS") +
  scale_color_manual(values = c("Parallel" = "blue", "Sequential" = "red"))

ggsave("./analysis/fps_comparison.png", fps_comparison, width = 12, height = 10)

# Calculate overall statistics
overall_stats <- combined_data %>%
  summarise(
    avg_speedup = mean(Speedup),
    max_speedup = max(Speedup),
    avg_efficiency = mean(Efficiency),
    max_efficiency = max(Efficiency),
    avg_relative_performance = mean(RelativePerformance),
    max_relative_performance = max(RelativePerformance)
  )

# Find best and worst cases
best_speedup <- combined_data %>%
  filter(Speedup == max(Speedup)) %>%
  select(Zoom, MaxIterations, Speedup)

worst_speedup <- combined_data %>%
  filter(Speedup == min(Speedup)) %>%
  select(Zoom, MaxIterations, Speedup)

best_efficiency <- combined_data %>%
  filter(Efficiency == max(Efficiency)) %>%
  select(Zoom, MaxIterations, Efficiency)

worst_efficiency <- combined_data %>%
  filter(Efficiency == min(Efficiency)) %>%
  select(Zoom, MaxIterations, Efficiency)

# Write analysis results to a text file
sink("./analysis/performance_analysis.txt")

cat("Julia-Mandelbrot Screensaver Performance Analysis\n")
cat("================================================\n\n")

cat("Overall Statistics:\n")
cat("------------------\n")
cat(sprintf("Average Speedup: %.2f\n", overall_stats$avg_speedup))
cat(sprintf("Maximum Speedup: %.2f\n", overall_stats$max_speedup))
cat(sprintf("Average Efficiency: %.2f\n", overall_stats$avg_efficiency))
cat(sprintf("Maximum Efficiency: %.2f\n", overall_stats$max_efficiency))

cat("Best Cases:\n")
cat("----------\n")
cat("Best Speedup:\n")
print(best_speedup)
cat("\nBest Efficiency:\n")
print(best_efficiency)
cat("\n")

cat("Worst Cases:\n")
cat("------------\n")
cat("Worst Speedup:\n")
print(worst_speedup)
cat("\nWorst Efficiency:\n")
print(worst_efficiency)
cat("\n")
