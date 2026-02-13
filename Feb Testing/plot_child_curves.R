# Plot offspring curves from simulation formulas
# Formulas match Child_RunModel.jl / child_testing.jl:
#   parabola:   peak at opt=0.65, width=0.70, floor 1.0
#   asymptote:  max_children=4, k=4.2
#   linear:     min=1, max=4

if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

if (!dir.exists("graphs")) dir.create("graphs", recursive = TRUE)

# --- Formulas (same as Julia) ---
# Parabola: offspring peaks at intermediate compatibility
offspring_parabola <- function(score) {
  peak <- 4.0
  opt <- 0.65
  width <- 0.70
  scale <- 1 - ((score - opt) / width)^2
  base <- peak * max(0, scale)
  max(1.0, base)
}

# Asymptote: rises then levels off
offspring_asymptote <- function(score) {
  max_children <- 4.0
  k <- 4.2
  max_children * (1 - exp(-k * score))
}

# Linear: monotonically increasing
offspring_linear <- function(score) {
  min_offspring <- 1.0
  max_offspring <- 4.0
  min_offspring + (max_offspring - min_offspring) * score
}

# --- Data: compatibility score in [0, 1] ---
scores <- seq(0, 1, length.out = 501)
df <- data.frame(
  score = rep(scores, 3),
  offspring = c(
    sapply(scores, offspring_parabola),
    sapply(scores, offspring_asymptote),
    sapply(scores, offspring_linear)
  ),
  model = rep(c("Parabola", "Asymptote", "Linear"), each = length(scores))
)
df$model <- factor(df$model, levels = c("Parabola", "Asymptote", "Linear"))

# --- Combined plot: all three curves ---
p_combined <- ggplot(df, aes(x = score, y = offspring, color = model, linetype = model)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Parabola" = "#E69F00", "Asymptote" = "#009E73", "Linear" = "#0072B2")
  ) +
  scale_linetype_manual(
    values = c("Parabola" = "solid", "Asymptote" = "dashed", "Linear" = "dotdash")
  ) +
  labs(
    title = "Offspring curves (from simulation formulas)",
    subtitle = "Expected offspring count vs compatibility score [0,1]",
    x = "Compatibility score",
    y = "Expected offspring count",
    color = "Child model",
    linetype = "Child model"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 5)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("graphs/plot_child_curves_combined.png", p_combined, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_child_curves_combined.png\n")

# --- Individual plots per model ---
for (mod in c("Parabola", "Asymptote", "Linear")) {
  sub <- df[df$model == mod, ]
  col <- c("Parabola" = "#E69F00", "Asymptote" = "#009E73", "Linear" = "#0072B2")[mod]
  p <- ggplot(sub, aes(x = score, y = offspring)) +
    geom_line(color = col, linewidth = 1.2) +
    labs(
      title = paste0(mod, " offspring curve"),
      subtitle = "Formula from Child_RunModel.jl",
      x = "Compatibility score",
      y = "Expected offspring count"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 5)) +
    theme_minimal()
  fname <- paste0("graphs/plot_child_curve_", tolower(mod), ".png")
  ggsave(fname, p, width = 8, height = 6, dpi = 300)
  cat("Saved:", fname, "\n")
}

cat("Done. Curves plotted from formulas in Child_RunModel.jl.\n")
