#!/usr/bin/env Rscript
# offspring_mate_curves.R
# Plot offspring-vs-mates curves for logistic, exponential decay, and Gaussian models.

args<-commandArgs(trailingOnly=TRUE)
# default parameters
max_offspring <- 5
r <- 1.4
c <- 3
mu <- 3
sigma <- 1
max_mates <- 10
outfile <- "offspring_mate_curves.png"

if(length(args) > 0) {
  params <- as.list(args)
  for (p in params) {
    kv <- strsplit(p, "=")[[1]]
    if(length(kv)==2) {
      key <- kv[1]
      val <- as.numeric(kv[2])
      if(!is.na(val)) {
        if(key == "max_offspring") max_offspring <- val
        if(key == "r") r <- val
        if(key == "c") c <- val
        if(key == "mu") mu <- val
        if(key == "sigma") sigma <- val
        if(key == "max_mates") max_mates <- val
      } else {
        if(key == "outfile") outfile <- kv[2]
      }
    }
  }
}

mates <- 0:max_mates

logistic <- max_offspring/(1 + exp(-r*(mates-c)))
expdecay <- max_offspring * exp(-r*mates)
gaussian <- max_offspring * exp(-((mates-mu)^2)/(2*sigma^2))

library(ggplot2)

df <- data.frame(mates, logistic, expdecay, gaussian)

df_long <- stack(df[, c("logistic", "expdecay", "gaussian")])
df_long$mates <- rep(mates, times=3)
colnames(df_long) <- c("offspring", "model", "mates")

p <- ggplot(df_long, aes(x=mates, y=offspring, color=model)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  labs(title="Offspring vs Mates Curves",
       x="Number of mates",
       y="Expected offspring") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max_offspring*1.1)) +
  theme(legend.title = element_blank())

ggsave(outfile, p, width=8, height=5)
cat("Saved plot to", outfile, "\n")
