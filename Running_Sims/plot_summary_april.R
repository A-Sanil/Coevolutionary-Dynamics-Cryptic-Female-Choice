library(ggplot2)
library(readr)
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_summary_april.R <summary_csv>")
}
summary_csv = args[1]
df = read_csv(summary_csv)
# Plot MeanMale and MeanFemale with sd ribbons
p1 = ggplot(df, aes(x=Generation)) +
  geom_ribbon(aes(ymin=MeanMale_mean-MeanMale_sd, ymax=MeanMale_mean+MeanMale_sd), fill='blue', alpha=0.2) +
  geom_line(aes(y=MeanMale_mean), color='blue') +
  geom_ribbon(aes(ymin=MeanFemale_mean-MeanFemale_sd, ymax=MeanFemale_mean+MeanFemale_sd), fill='red', alpha=0.2) +
  geom_line(aes(y=MeanFemale_mean), color='red') +
  labs(title='Mean Male & Female by Generation', y='Trait mean') +
  theme_minimal()

p2 = ggplot(df, aes(x=Generation)) +
  geom_ribbon(aes(ymin=MeanCount_mean-MeanCount_sd, ymax=MeanCount_mean+MeanCount_sd), fill='darkgreen', alpha=0.2) +
  geom_line(aes(y=MeanCount_mean), color='darkgreen') +
  labs(title='Mean Sperm Count by Generation', y='Mean sperm count') +
  theme_minimal()

outdir = dirname(summary_csv)
png(filename=file.path(outdir, 'mean_traits_april.png'), width=900, height=600)
print(p1)
dev.off()

png(filename=file.path(outdir, 'mean_spermcount_april.png'), width=900, height=600)
print(p2)
dev.off()
