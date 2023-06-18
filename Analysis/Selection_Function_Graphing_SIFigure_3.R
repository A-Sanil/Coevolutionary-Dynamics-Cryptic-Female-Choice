# R script to generate SI Figure 3
#for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#please send any questions to mkustra@ucsc.edu
#visualizes how strength of selection, male phenotypes, and female phenotypes influence probability of fertilization sucess.

#Loading up libraries
library(patchwork)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggnewscale)

# function that determines probability of fertilization success, assuming sperm number of two males is equal
#mphen2 is male 2 phenotype
#mpehn1 is male 1 phenotype
#a is strength of selection
#d is the female phenotype
prob_success<-function(mphen2,mphen1,a,d){
  probs_m1 <- exp((-(mphen2 - d)^2) / (2*a)) 
  probs_m2 <- exp((-(mphen1 - d)^2) / (2*a))
  return(probs_m1 / (probs_m1 + probs_m2))
}

# function used to make the plots
#mphen is male phenotype
#opt is optimal trait value
plot_function <- function(mphen, opt) {
  #make dataframe for vertical lines
  vline <-
    data.frame(
      Lines = c("Female Trait Value", "Competitor Trait Value"),
      Xint = c(opt, mphen)
    )
  p <-
    ggplot(data = data.frame(mphen2 = 0),
           mapping = aes(mphen2 = mphen2)) +
    theme_classic() +
    xlab("Male Competitor Phenotype") +
    ylab("Probability of Fertilization") +
    xlim(25,75) +
    theme(text = element_text(size = 15)) +
    labs(color = "Functions") +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.1))
  #now add lines for different strengths of selection
  p + stat_function(
    fun = prob_success,
    args = (list(
      mphen1 = mphen, a = 1, d = opt
    )),
    mapping = aes(color = "Strong"),
    size = 1.5
  ) +
    stat_function(
      fun = prob_success,
      args = (list(
        mphen1 = mphen, a = 12.5, d = opt
      )),
      mapping = aes(color = "Moderate"),
      size = 1.5
    ) +
    stat_function(
      fun = prob_success,
      args = (list(
        mphen1 = mphen, a = 50, d = opt
      )),
      mapping = aes(color = "Weak"),
      size = 1.5
    ) +
    scale_color_manual(
      name = "Strength of Selection:",
      breaks = c("Weak", "Moderate", "Strong"),
      values = c("#E69f00", "#000000", "#009E73")
    ) +
    new_scale_color() +
    scale_color_manual(values = c("#619CFF", "#F8766D"), name = "") +
    geom_vline(
      data = vline,
      aes(
        xintercept = Xint,
        linetype = Lines,
        color = Lines
      ),
      size = 1,
      alpha = 0.6
    ) +
    scale_linetype_manual(values = c("dashed", "dashed"), name = "")
}


# Making SIFigure 2
# Making panel A
low_high <- plot_function(45, 50) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")
# Making Panel B
med_med <- plot_function(50, 50) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")
# Making Panel C
high_low <- plot_function(55, 50) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")
# Making the Y axis label
y.label <- textGrob(
  "Probability of \n Fertilization Success",
  gp = gpar(
    fontface = "bold",
    col = "black",
    fontsize = 15
  ),
  rot = 90
)
# Making the X axis label
x.label <- textGrob("Focal Male Sperm Trait Value",
                    gp = gpar(
                      fontface = "bold",
                      col = "black",
                      fontsize = 15
                    ))
# Making a blank plot to help with spacing for aesthetic reasons
blank_p <- plot_spacer() + theme_void()
# Extracting the legend
legend <- get_legend(high_low + theme(
  legend.position = "right",
  legend.box.margin = margin(0, 0, 0, 12)
))
# Making the first plots of A,B,C
plot <- plot_grid(low_high,
                  med_med,
                  high_low,
                  ncol = 3,
                  labels = c("A)", "B)", "C)"))
# Putting together the final plot with the axis label
plot2 <-
  grid.arrange(arrangeGrob(plot, left = y.label, bottom = x.label))
# Adding the legend to the plot with axis labels
SIFigure_2 <-
  plot_grid(plot2,
            NULL,
            legend,
            rel_widths = c(3,-0.2, 1),
            nrow = 1)

SIFigure_2
# Save the plot for SIfigure 2
#save_plot("SIFigure_2.pdf", SIFigure_2, ncol = 2)

#Reviewer repsonse strong selction
#median standard deviation for strong selection is 1.5 when risk of sperm competition is low.
plot_function2 <- function(mphen, opt) {
  #make dataframe for vertical lines
  vline <-
    data.frame(
      Lines = c("Female Trait Value", "Competitor Trait Value"),
      Xint = c(opt, mphen)
    )
  p <-
    ggplot(data = data.frame(mphen2 = 0),
           mapping = aes(mphen2 = mphen2)) +
    theme_classic() +
    xlab("Male Competitor Phenotype") +
    ylab("Probability of Fertilization") +
    xlim(50-3,50+3) +
    theme(text = element_text(size = 15)) +
    labs(color = "Functions") +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.1))
  #now add lines for different strengths of selection
  p + stat_function(
    fun = prob_success,
    args = (list(
      mphen1 = mphen, a = 1, d = opt
    )),
    mapping = aes(color = "Strong"),
    size = 1.5
  ) +
    stat_function(
      fun = prob_success,
      args = (list(
        mphen1 = mphen, a = 12.5, d = opt
      )),
      mapping = aes(color = "Moderate"),
      size = 1.5
    ) +
    stat_function(
      fun = prob_success,
      args = (list(
        mphen1 = mphen, a = 50, d = opt
      )),
      mapping = aes(color = "Weak"),
      size = 1.5
    ) +
    scale_color_manual(
      name = "Strength of Selection:",
      breaks = c("Weak", "Moderate", "Strong"),
      values = c("#E69f00", "#000000", "#009E73")
    ) +
    new_scale_color() +
    scale_color_manual(values = c("#619CFF", "#F8766D"), name = "") +
    geom_vline(
      data = vline,
      aes(
        xintercept = Xint,
        linetype = Lines,
        color = Lines
      ),
      size = 1,
      alpha = 0.6
    ) +
    scale_linetype_manual(values = c("dashed", "dashed"), name = "")
}

plot_function2(50, 50) +
  theme(legend.position = "none",text=element_text(size=16)) +
  xlab("Focal male fertilization success") +
  ylab("Probability of fertilization success")+
  geom_vline(aes(xintercept=50-1.5),size=1,color="black",alpha=0.5,linetype="dotted")+
  geom_vline(aes(xintercept=50+1.5),size=1,color="black",alpha=0.5,linetype="dotted")
