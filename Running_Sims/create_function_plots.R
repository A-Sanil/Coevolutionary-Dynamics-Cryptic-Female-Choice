# R script to create visualization images for the three offspring functions
# Shows how each function should look: parabola, linear (monotonic), and asymptotic

# Try to load ggplot2, install if needed
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
  library(ggplot2)
}

# Define the functions matching the Julia code
offspring_parabola <- function(score) {
  peak <- 4.0
  opt <- 0.65
  width <- 0.70
  scale <- 1 - ((score - opt) / width)^2
  base <- peak * pmax(0.0, scale)
  return(pmax(1.0, base))
}

offspring_asymptote <- function(score) {
  max_children <- 4.0
  k <- 4.2
  return(max_children * (1 - exp(-k * score)))
}

offspring_linear <- function(score) {
  min_offspring <- 1.0
  max_offspring <- 4.0
  return(min_offspring + (max_offspring - min_offspring) * score)
}

# Create data for plotting
scores <- seq(0, 1, by = 0.01)
parabola_vals <- offspring_parabola(scores)
asymptote_vals <- offspring_asymptote(scores)
linear_vals <- offspring_linear(scores)

df <- data.frame(
  score = rep(scores, 3),
  offspring = c(parabola_vals, asymptote_vals, linear_vals),
  function_type = rep(c("Parabola", "Asymptotic", "Linear (Monotonic)"), each = length(scores))
)

# Create the plot
p <- ggplot(df, aes(x = score, y = offspring, color = function_type, linetype = function_type)) +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Offspring Functions: Expected Shapes",
    x = "Compatibility Score",
    y = "Expected Offspring Count",
    color = "Function Type",
    linetype = "Function Type"
  ) +
  scale_color_manual(values = c("Parabola" = "#E69F00", "Asymptotic" = "#009E73", "Linear (Monotonic)" = "#0072B2")) +
  scale_linetype_manual(values = c("Parabola" = "solid", "Asymptotic" = "dashed", "Linear (Monotonic)" = "dotdash")) +
  xlim(0, 1) +
  ylim(0, 5)

# Save the plot
ggsave("function_shapes_visualization.png", p, width = 12, height = 8, dpi = 300)
cat("Saved: function_shapes_visualization.png\n")

# Also create individual plots for each function
# Parabola
p1 <- ggplot(data.frame(score = scores, offspring = parabola_vals), aes(x = score, y = offspring)) +
  geom_line(color = "#E69F00", size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 14), plot.title = element_text(size = 16, face = "bold")) +
  labs(title = "Parabola Function", x = "Compatibility Score", y = "Expected Offspring Count") +
  xlim(0, 1) + ylim(0, 5)
ggsave("function_parabola.png", p1, width = 8, height = 6, dpi = 300)

# Asymptotic
p2 <- ggplot(data.frame(score = scores, offspring = asymptote_vals), aes(x = score, y = offspring)) +
  geom_line(color = "#009E73", size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 14), plot.title = element_text(size = 16, face = "bold")) +
  labs(title = "Asymptotic Function", x = "Compatibility Score", y = "Expected Offspring Count") +
  xlim(0, 1) + ylim(0, 5)
ggsave("function_asymptotic.png", p2, width = 8, height = 6, dpi = 300)

# Linear (Monotonic)
p3 <- ggplot(data.frame(score = scores, offspring = linear_vals), aes(x = score, y = offspring)) +
  geom_line(color = "#0072B2", size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 14), plot.title = element_text(size = 16, face = "bold")) +
  labs(title = "Linear (Monotonically Increasing) Function", x = "Compatibility Score", y = "Expected Offspring Count") +
  xlim(0, 1) + ylim(0, 5)
ggsave("function_linear.png", p3, width = 8, height = 6, dpi = 300)

cat("All function visualization images created successfully!\n")

