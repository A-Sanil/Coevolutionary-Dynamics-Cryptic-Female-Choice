#SI web app for the paper: Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#please send any questions to mkustra@ucsc.edu
# 1. Getting everything ready ---------------------------------------------

# * 1.a Loading up libraries ----------------------------------------------
library(shiny)
library(ggh4x)
library(shinythemes)
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(profvis)
library(dplyr)
library(tidyr)
library(viridis)
library(pracma)
# * 1.b Loading up and processing data ------------------------------------
run_data <- readRDS("data/run_data_shiny_shorter.rds") %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )
run_data$Level <- factor(run_data$Level)
summary(run_data)
#predicted values
PredAn<-data.frame(Level=c("0.25","0.5","0.75","1.0"),Tradeoff=c("false","false","false","false","true","true","true","true"),PredN=c(16.10988593067346,23.66053786200191,28.9673742726416,33.18317009243252),PredT=c(16.10988593067346,23.66053786200191,28.9673742726416,33.18317009243252,805.4942965336729,1183.0268931000955,1448.36871363208,1659.1585046216262))

run_data <-merge(run_data,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=ifelse(type=="Stabilizing",MeanMale-50,MeanMale-MeanFemale),investScale=ifelse(Tradeoff=="true",investment/50,investment))%>%
  select(-investment,-diff,-diff2)
#Setting up the labeling
A.labs <- c("Weak", "Moderate", "Strong")
names(A.labs) <- c("50.0","12.5","1.0")
Var.labs <- c("Small Var", "Medium Var", "Large Var")
names(Var.labs) <- c("LV", "MV", "HV")
T.labs <- c("No Tradeoff", "Tradeoff")
names(T.labs) <- c("false", "true")
L.labs <- c("2 Loci", "20 Loci")
names(L.labs) <- c("2", "20")
P.labs <- c("N = 1000", "N = 10000")
names(P.labs) <- c("1000", "10000")
#pallete for ploting
pal <- rev(plasma(8)[c(1, 3, 5, 7)])
#default theme
mytheme <-
  theme_bw() + theme(
    legend.position = "bottom",
    #this puts legend on the bottom
    axis.title = (element_text(face = "bold")),
    #Makes the axis line black and  thicker
    text = element_text(size = 15, face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgrey")
  )#makes all the text larger and bold
theme_set(mytheme)

#summarizing data
SumRuns <- run_data %>%
  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize,type) %>%
  mutate(
    intercor = cor(MeanMale, MeanFemale),
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  ) %>%
  summarise_all(funs(mean = mean, sd = sd), na.rm = T)

#Filtering out data for last generation
runs_rawS_LG <- run_data %>%
  filter(Generation == 30000) %>%
  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize,type) %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )

# * 1.c global functions --------------------------------------------------
#analytical model
prob_mate<-function(sperm,alpha,beta){
  return((1-(1/(1+exp(-alpha*(sperm-beta))))))
}
solution<-function(q,alpha,beta){
  (q/2 + lambertWp(q*exp(alpha*beta - q/2)/2))/alpha
  
}

#function that makes plot
plot_function_an<-function(alpha,beta){
  #make data frame for vertical lines
  #make plot
  p1 <-
    ggplot(data = data.frame(sperm = 0),
           mapping = aes(sperm = sperm)) + theme_classic() + xlab("Ejaculate investment") +
    ylab("Probability of mating success") + xlim(0, beta*2) + theme(text = element_text(size =25)) + labs(color = "Functions") + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),limits = c(0, 1)) + stat_function(
      fun = prob_mate,
      args = (list(
        alpha = alpha,
        beta = beta
      )),
      size = 1.5
    )
  
  p2 <-
    ggplot(data = data.frame(q = 0),
           mapping = aes(q = q)) + theme_classic() + xlab("Risk of sperm competition") +
    ylab("ESS ejaculate allocation")  + theme(text = element_text(size =25)) + labs(color = "Functions") + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),limits = c(0, 1)) + stat_function(
      fun = solution,
      args = (list(
        alpha = alpha,
        beta = beta
      )),
      size = 1.5
    )
  return(p1+p2)
}
# * 1.c.1 Strength of selection -------------------------------------------
#function that determins probability of fertilization

prob_success<-function(mphen2,mphen1,sperm2,sperm1,a,d){
  probs_m1 <- exp((-(mphen2 - d)^2) / (2*a)) * sperm2 
  probs_m2 <- exp((-(mphen1 - d)^2) / (2*a)) * sperm1
  return(probs_m1 / (probs_m1 + probs_m2))
}

#function that makes plot
plot_function<-function(mphen,sperm2,sperm1,opt){
  #make data frame for vertical lines
  vline <-
    data.frame(
      Lines = c("Female Choice Trait Value", "Competitor Sperm Trait Value"),
      Xint = c(opt, mphen)
    )
  #make plot
  p <-
    ggplot(data = data.frame(mphen2 = 0),
           mapping = aes(mphen2 = mphen2)) + theme_classic() + xlab("Focal Sperm Trait Value") +
    ylab("Probability of Fertilization") + xlim(30, 99) + theme(text = element_text(size =25)) + labs(color = "Functions") + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),limits = c(0, 1))
  p + stat_function(
    fun = prob_success,
    args = (list(
      mphen1 = mphen,
      sperm2 = sperm2,
      sperm1 = sperm1,
      a = 1,
      d = opt
    )),
    mapping = aes(color = "Strong"),
    size = 1.5
  ) + stat_function(
    fun = prob_success,
    args = (list(
      mphen1 = mphen,
      sperm2 = sperm2,
      sperm1 = sperm1,
      a = 12.5,
      d = opt
    )),
    mapping = aes(color = "Moderate"),
    size = 1.5
  ) + stat_function(
    fun = prob_success,
    args = (list(
      mphen1 = mphen,
      sperm2 = sperm2,
      sperm1 = sperm1,
      a = 50,
      d = opt
    )),
    mapping = aes(color = "Weak"),
    size = 1.5
  ) + scale_color_manual(
    name = "Strength of Selection:",
    breaks = c("Weak", "Moderate", "Strong"),
    values = c("#E69f00", "#000000", "#009E73")
  ) + new_scale_color() + geom_vline(
    data = vline,
    aes(
      xintercept = Xint,
      linetype = Lines,
      color = Lines
    ),
    alpha = 0.5,
    size = 1.5
  ) + scale_color_manual(values = c("blue", "red"), name = "") + scale_linetype_manual(values =c("dashed", "dotted"), name = "")
}

# * 1.c.2 Run Data graphing function --------------------------------------------------------
#input$alpha_r,input$all_r,input$populations,input$rsc_r,run_data
runs_plot<-function(means,alpha,allpop,pops,rsc,Loc,pop,tr,v,ty,Data,Pred=PredAn){
  rsc2<-as.factor(rsc)
  PredN<-Pred[Pred$Level == rsc2 &
                   Pred$Tradeoff == tr,]$PredN
  if(allpop=="yes"){
    #mean traits
    if(means=="Mean"){
    mt <-ggplot(data = Data[Data$A == alpha &
                           Data$Level == rsc &
                           Data$Loci == Loc &
                           Data$PopSize == pop &
                           Data$Tradeoff == tr &
                           Data$Var == v &
                           Data$type==ty, ], aes(x = Generation, group = Rep)) + theme_classic() +
      theme(
        legend.position = "top",
        axis.title = (element_text(
          face = "bold", margin = margin(r = 30)
        )),
        axis.line = element_line(color = "black", size = 2),
        text = element_text(size = 15, face = "bold"),
        strip.background = element_rect(color = "white"),
        strip.placement = "inside",
        strip.text = element_text(size = 15)
      ) + xlab("") + ylab("Average Trait Value") + geom_line(aes(y = MeanMale, color ="Male"), size = 0.3) +
      geom_line(aes(y = MeanFemale, color = "Female")) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_manual(name = "Sex", values = c("#F8766D", "#619CFF"))
    if(ty=="Stabilizing"){
      mt <-mt+geom_hline(aes(yintercept=50),linetype="dashed",alpha=0.5,color="black",size=2)
    }
    #mean count
    mn<-ggplot(data = Data[Data$A == alpha &
                             Data$Level == rsc &
                             Data$Loci == Loc &
                             Data$PopSize == pop &
                             Data$Tradeoff == tr &
                             Data$Var == v &
                             Data$type==ty, ], aes(x = Generation, group = Rep)) + theme_classic() +
      theme(
        legend.position = "top",
        axis.title = (element_text(
          face = "bold", margin = margin(r = 30)
        )),
        axis.line = element_line(color = "black", size = 2),
        text = element_text(size = 15, face = "bold"),
        strip.background = element_rect(color = "white"),
        strip.placement = "inside",
        strip.text = element_text(size = 15)
      ) + xlab("") + ylab("Average Sperm Number") +
      geom_hline(aes(yintercept=PredN),linetype="dashed",alpha=0.5,color="black",size=2)+
      geom_line(aes(y = MeanCount),color="#619CFF")
    }else{
    #CV traits
    sdt<-ggplot(data = Data[Data$A == alpha &
                              Data$Level == rsc &
                              Data$Loci == Loc &
                              Data$PopSize == pop &
                              Data$Tradeoff == tr &
                              Data$Var == v &
                              Data$type==ty, ], aes(x = Generation, group = Rep)) + theme_classic() +
      theme(
        legend.position = "top",
        axis.title = (element_text(
          face = "bold", margin = margin(r = 30)
        )),
        axis.line = element_line(color = "black", size = 2),
        text = element_text(size = 15, face = "bold"),
        strip.background = element_rect(color = "white"),
        strip.placement = "inside",
        strip.text = element_text(size = 15)
      ) + xlab("") + ylab("CV Trait Value") + geom_line(aes(y = CVMale, color ="Male"), size = 0.3) +
      geom_line(aes(y = CVFemale, color = "Female")) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_manual(name = "Sex", values = c("#F8766D", "#619CFF"))
    
    #CV count
    sdn<-ggplot(data = Data[Data$A == alpha &
                              Data$Level == rsc &
                              Data$Loci == Loc &
                              Data$PopSize == pop &
                              Data$Tradeoff == tr &
                              Data$Var == v &
                              Data$type==ty, ], aes(x = Generation, group = Rep)) + theme_classic() +
      theme(
        legend.position = "top",
        axis.title = (element_text(
          face = "bold", margin = margin(r = 30)
        )),
        axis.line = element_line(color = "black", size = 2),
        text = element_text(size = 15, face = "bold"),
        strip.background = element_rect(color = "white"),
        strip.placement = "inside",
        strip.text = element_text(size = 15)
      ) + xlab("") + ylab("CV Sperm Number") +
      geom_line(aes(y = CVCount),color="#619CFF")
    }
    
  }else{#graphing only individual runs
    #Filter Data and rearrange data to give male and female different colors on same graph
    trD<-Data[Data$A == alpha &
                Data$Level == rsc &
                Data$Loci == Loc &
                Data$PopSize == pop &
                Data$Tradeoff == tr & Data$Var == v & Data$Rep %in% pops &
                Data$type==ty, ] %>%
      select(Generation,
             Rep,
             MeanMale,
             MeanFemale,
             CVMale,
             CVFemale,
             MeanCount,
             CVCount) %>%
      rename(
        Mean_Male = MeanMale,
        Mean_Female = MeanFemale,
        CV_Male = CVMale,
        CV_Female = CVFemale
      ) %>%
      pivot_longer(
        !c(Generation, Rep, MeanCount, CVCount),
        names_to = c(".value", "Sex"),
        names_sep = "_"
      ) %>% mutate(
        MeanCount = ifelse(Sex == "Female", NA, MeanCount),
        CVCount = ifelse(Sex == "Female", NA, CVCount)
      )
    
    #mean traits
    if(means=="Mean"){
    mt<-ggplot(data = trD, aes(
      x = Generation,
      y = Mean,
      color = interaction(Sex, Rep)
    )) + theme_classic() + theme(
      legend.position = "top",
      axis.title = (element_text(
        face = "bold", margin = margin(r = 30)
      )),
      axis.line = element_line(color = "black", size = 2),
      text = element_text(size = 15, face = "bold"),
      strip.background = element_rect(color = "white"),
      strip.placement = "inside",
      strip.text = element_text(size = 15)
    ) + xlab("") + ylab("Average Trait Value") + geom_line(aes(group = interaction(Sex, Rep)), size =0.3) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_discrete(name ="Sex.Population #")
    if(ty=="Stabilizing"){
      mt <-mt+geom_hline(aes(yintercept=50),linetype="dashed",alpha=0.5,color="black",size=2)
    }
    #mean count
    mn<-ggplot(data = trD, aes(
      x = Generation,
      y = MeanCount,
      color = interaction(Sex, Rep)
    )) + theme_classic() + theme(
      legend.position = "top",
      axis.title = (element_text(
        face = "bold", margin = margin(r = 30)
      )),
      axis.line = element_line(color = "black", size = 2),
      text = element_text(size = 15, face = "bold"),
      strip.background = element_rect(color = "white"),
      strip.placement = "inside",
      strip.text = element_text(size = 15)
    ) + xlab("") + ylab("Average Sperm Number") +
      geom_hline(aes(yintercept=PredN),linetype="dashed",alpha=0.5,color="black",size=2)+ geom_line(aes(group = interaction(Sex, Rep)), size =0.3) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_discrete(name = "Sex.Population #")
    
    }else{
    #CV traits
    sdt<-ggplot(data = trD, aes(
      x = Generation,
      y = CV,
      color = interaction(Sex, Rep)
    )) + theme_classic() + theme(
      legend.position = "top",
      axis.title = (element_text(
        face = "bold", margin = margin(r = 30)
      )),
      axis.line = element_line(color = "black", size = 2),
      text = element_text(size = 15, face = "bold"),
      strip.background = element_rect(color = "white"),
      strip.placement = "inside",
      strip.text = element_text(size = 15)
    ) + xlab("") + ylab("CV Trait Value") + geom_line(aes(group = interaction(Sex, Rep)), size =0.3) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_discrete(name ="Sex.Population #")
    #CV count
    sdn<-ggplot(data = trD, aes(
      x = Generation,
      y = CVCount,
      color = interaction(Sex, Rep)
    )) + theme_classic() + theme(
      legend.position = "top",
      axis.title = (element_text(
        face = "bold", margin = margin(r = 30)
      )),
      axis.line = element_line(color = "black", size = 2),
      text = element_text(size = 15, face = "bold"),
      strip.background = element_rect(color = "white"),
      strip.placement = "inside",
      strip.text = element_text(size = 15)
    ) + xlab("") + ylab("CV Sperm Number") + geom_line(aes(group = interaction(Sex, Rep)), size =0.3) + guides(color = guide_legend(override.aes = list(size = 5))) + scale_color_discrete(name ="Sex.Population #")
    }
  }
  #due to memory limits, can unfortunatly only show either mean or sd at once :(
  if(means=="Mean"){#show mean
    p<-((mt+mn)&theme(legend.position = "top"))+plot_layout(guides="collect")
  }else{#show CV
    p<-((sdt+sdn)&theme(legend.position = "top"))+plot_layout(guides="collect")
  }
  return(p)
}
# 2. User interface -------------------------------------------------------
ui <- fluidPage(theme = shinytheme("superhero"),
  HTML('<meta name="viewport" content="width=1024">'),
  navbarPage(
    "",

# * 2.a Homepage ---------------------------------------------------------
    tabPanel(
      "Homepage",
      fluidRow(
      titlePanel(h1(
        "Supporting information for:", align = "center"
      ))),fluidRow(h1("'The coevolutionary dynamics of cryptic female choice'", align = "center")),
        fluidRow(
          column(12,
          h2("Abstract:"),
          br(),
          p(
            "In contrast to sexual selection on traits that affect interactions between the sexes before mating, little theoretical research has focused on the coevolution of post-mating traits via cryptic female choice (when females bias fertilization toward specific males). We use simulation models to ask (1) whether and if so how non-directional cryptic female choice (female-by-male interactions in fertilization success) causes deviations from models that focus exclusively on male-mediated post-mating processes and (2) how the risk of sperm competition, the strength of cryptic female choice, and tradeoffs between sperm number and sperm traits interact to influence the coevolutionary dynamics between cryptic female choice and sperm traits. We found that incorporating cryptic female choice can result in males investing much less in their ejaculates than predicted by models with sperm competition only. We also found that cryptic female choice resulted in the evolution of genetic correlations between cryptic female choice and sperm traits, even when the strength of cryptic female choice was weak, and the risk of sperm competition was low. This suggests that cryptic female choice may be important even in systems with low multiple mating. These genetic correlations increased with the risk of sperm competition and as the strength of cryptic female choice increased. When the strength of cryptic female choice and risk of sperm competition was high, extreme co-divergence of sperm traits and cryptic female choice preference occurred even when the sperm trait traded off with sperm number. We also found that male traits lagged behind the evolution of female traits; this lag decreased with increasing strength of cryptic female choice and risk of sperm competition. Overall, our results suggest that cryptic female choice deserves more attention theoretically and may be driving trait evolution in ways just beginning to be explored."
          ))),
      fluidRow(
      column(12, br(),
             h2("Description:"),
             br(),
             p(
               "This is a supplemental web application for the paper: 'The coevolutionary dynamics of cryptic female choice.' At each tab, one can make supplemental figures using data from this paper. One can get to different pages of the web application by using the tabs above. Description of the figures are given on each page, and short descriptions of each tab are given below. You can save graphs by right clicking on the graph."
             ))),
        fluidRow(
          column(4,
                 br(),
                 h4("Analytical sperm competition model:"),
                 br(),
                 p(
                   "Here one can see how the alpha and beta parameter influence the tradeoff between postmating and premating sucess, and the evolutionary stable strategy for ejaculate investment."
                 )),
          column(4,
          br(),
          h4("Strength of selection:"),
          br(),
          p(
            "Here one can see how the strength of selection, cryptic female choice trait, male sperm trait, and sperm number influence the probability of fertilization in the model."
          )),
          column(4,
                 br(),
          h4("Graphing model runs:"),
          br(),
          p(
            "Here one can visualize the evolution of sperm number, male sperm trait, and cryptic female choice trait under different parameter values at the detail of individual runs (i.e., populations). The data for these graphs were used for data analyses in the paper."
          )
        )),
      fluidRow(
        column(4,
               br(),
               h4("Average trait changes:"),
               br(),
               p(
                 "Here one can visualize the evolution of sperm number, male sperm trait, and cryptic female choice trait averaged across all 50 runs of the same unique parameter combinations. The user can see how these traits evolve overtime or distributions of these traits at specific generations. The data for these graphs were used for data analyses in the paper."
               )),
        column(4,
               br(),
               h4("Average trait deviations:"),
               br(),
               p(
                 "Here one can visualize how sperm traits deviate from the optimum (only considering selection on the sperm trait) or how ejaculate allocation differs from ejaculate allocation predicted by the analytical model. Compare these figures to Figure 1B and 3 in the main text."
               )),
    column(4, br(),
           h4("Correlations:"),
           br(),
           p(
             "Here one can visualize the average genetic within-population correlation between male sperm trait and cryptic female choice trait as well as across-population phenotypic correlations.The user can see how these correlations change over time or as distributions (within-population correlation) or scatter plots of the different populations (across-population correlation) of these traits at specific generations.The data for these graphs were used for data analyses in the paper."
          )
    ))),


# * 2.b Analytical model -------------------------------------------------
tabPanel(
  "Analytical sperm competition model",
  sidebarPanel(h3("Figure legend:"),"here you can explore how the alpha and beta parameter influence the tradeoff between premating and postmating sexual selection (left; Equation 1) and the evolutionary stable strategy (ESS) ejaculate investment (right; Equation 6). Alpha is the inverse (e.g., selected 20 is plotting when alpha = 1/20). Compare to Figure S1. In all simlations without a tradeoff between sperm number and sperm trait alpha =1/20 and beta = 50;for simulations with a tradeoff alpha = 1/1000 beta = 2500 to keep the scaling of overall investment the same.",
               fluidRow(
                 numericInput(
                   "alpha",
                   h3("Alpha shape parameter:"),
                   min = 1,
                   max = 10000,
                   value = 20
                 ),
                 numericInput(
                   "beta",
                   h3("Beta shape parameter:"),
                   min = 1,
                   max = 1000000,
                   value = 50
                 )
               )
  ),
  mainPanel(plotOutput("fun_plot_an",height=600))
),



# * 2.c Strength of Selection ---------------------------------------------
    tabPanel(
      "Strength of selection",
      sidebarPanel(h3("Figure legend:"),"the probability of fertilization for the focal male is highest when his trait value matches the optimal trait value or female choice trait value in the case of matching selection (red dashed line; Eq. 9). This probability increases the farther the competitor sperm trait value (blue dashed line) deviates from the optimal trait value relative to the focal male’s trait value regardless of direction. The fertilization advantage of being closer to the optimal trait value increases as strength of selection increases. One can use the sliders to change the focal male sperm number, the competitor male sperm number, the competitor sperm trait value, and the focal female trait value. Optimal trait value can take any positive number. Compare to Figure S3.",
        fluidRow(
          sliderInput(
            "fm",
            h3("Competitor sperm trait value:"),
            min = 30,
            max = 99,
            value = 40,
            step = 0.01
          ),
          sliderInput(
            "om",
            h3("Focal female trait value (optimal sperm trait value):"),
            min = 30,
            max = 99,
            value = 40,
            step = 0.01
          ),
          sliderInput(
            "fs",
            h3("Focal male sperm number:"),
            min = 1,
            max = 700,
            value = 50,
            step = 10
          ),
          sliderInput(
            "cs",
            h3("Competitor male sperm number:"),
            min = 1,
            max = 700,
            value = 50,
            step = 10
          )
        )
      ),
      mainPanel(plotOutput("fun_plot",height=600))
    ),

# * 2.d Graphing model runs -----------------------------------------------
tabPanel(
  "Graphing model runs",
  titlePanel(
    h1("Graphing model runs", align = "center")
  ),
  sidebarPanel(fluidRow(h3("Figure legend:"),"Evolutionary dynamics of mean sperm trait, cryptic female choice trait, and sperm number and within-population coefficient of variation (CV) in trait values at a certain parameter combination chosen by using the buttons below. You can either display all populations at once (default), or select certain populations to visualize. When all populations are being displayed each line represents an individual run of the model (population), with blue being sperm trait and sperm number, red being cyptic female choice trait. When you select certain populations there will be a different color for a given population and whether the line is the male or female trait. Values are only displayed every 400 generations. Grey dashed line for sperm number is the sperm number predicted by the analytical model. Stabilizing selection runs are only available for 20 loci; moderate starting variation, and population sizes of 10,000; the grey dahsed line for average trait value is 50 representing the stabilizing selection optimum.",
    radioButtons(
      "means_r",
      h3("Show mean or CV?"),
      choices=list(
        "Mean",
        "CV"
      ),inline=T,
      selected="Mean"
    ),                    
    radioButtons(
      "alpha_r",
      h3("Strength of selection:"),
      choices = list(
        "50.0 (weak)" = "50.0",
        "12.5 (moderate)" = "12.5",
        "1.0 (strong)" = "1.0"
      ),inline=T,
      selected = "12.5"
    ),
    radioButtons(
      "rsc_r",
      h3("Risk of sperm competition:"),
      choices = list(
        "0.25" = "0.25",
        "0.5" = "0.5",
        "0.75" = "0.75",
        "1" = "1.0"
      ),inline=T,selected = "0.25"),
    radioButtons(
      "Tr_r",
      h3("Tradeoff:"),
      choices = list(
        "Tradeoff" = "true",
        "No Tradeoff" = "false"
      ),inline=T,selected = "true"),
    radioButtons(
      "type_r",
      h3("Type of selection:"),
      choices = list(
        "Cryptic female choice" = "Cryptic female choice",
        "Stabilizing selection" = "Stabilizing"
      ),inline=T,selected = "Cryptic female choice"),
    conditionalPanel(
      condition="input.type_r == 'Cryptic female choice'",
      radioButtons(
        "Loc_r",
        h3("Number of loci:"),
        choices = list(
          "2" = "2",
          "20" = "20"
        ),inline=T,selected = "20"),
      radioButtons(
        "Popsize_r",
        h3("Population size:"),
        choices = list(
          "N = 1000" = "1000",
          "N = 10,000" = "10000"
        ),inline=T,selected = "10000"),
      radioButtons(
        "V_r",
        h3("Starting variation:"),
        choices = list(
          "Small" = "LV",
          "Medium" = "MV",
          "Large"="HV"
        ),inline=T,
        selected = "MV")),
    radioButtons(
      "all_r",
      h3("Show all population runs?"),
      choices=list(
        "All Populations"="yes",
        "Only Selected Population(s)"="no"),inline=T,
      selected="yes"),
    conditionalPanel(
      condition="input.all_r == 'no'",
      checkboxGroupInput(
        "populations",
        h3("Populations to Show"),
        choices=list("1"=1,"2"=2,"3"=3,"4" =4,"5"=5,"6" = 6,"7" = 7, "8"=8,"9"= 9, "10"=10,"11"=11,"12"=12,"13"=13,"14" =14,"15"=15,"16" = 16,"17" = 17, "18"=18,"19"= 19, "20"=20,"21"=21,"22"=22,"23"=23,"24" =24,"25"=25,"26" = 26,"27" = 27, "28"=28,"29"= 29, "30"=30,"31"=31,"32"=32,"33"=33,"34" =34,"35"=35,"36" = 36,"37" = 37, "38"=38,"39"= 39, "40"=40,"41"=41,"42"=42,"43"=43,"44" =44,"45"=45,"46" = 46,"47" = 47, "48"=48,"49"= 49, "50"=50),
        selected=1,inline=TRUE     
      )
    )
  )
  ),
  mainPanel(plotOutput("runs",height=800))
),
# * 2.e Trait sensitivity ----------------------------------
    tabPanel(
        "Average trait changes",
  titlePanel(
    h1(
      "Average trait changes",
      align = "center"
   )),sidebarPanel(
     h3("Figure legend:"),"Evolutionary dynamics of the mean mean (above) and within-population coefficient of variation (CV; below) trait values across all 50 populations of a parameter combination. You can either display how the mean across 50 populations changes over time (default) or display boxplots of all 50 populations at specific generations. In the default ('Change over time') The solid line represents the mean and the shade represents the standard error in that mean. You can select which trait to show (cryptic female choice trait, sperm trait, or sperm number). The defualt setting shows the base model assumptions, but you can change the parameters (number of loci, population size, and starting variation) to show the results of the sensitivity analyses. Values are only displayed every 400 generations. Dashed lines for sperm number is the sperm number predicted by the analytical model at different risks of sperm competition. Stabilizing selection runs are only available for 20 loci; moderate starting variation, and population sizes of 10,000; the grey dahsed line for average trait value is 50 representing the stabilizing selection optimum.",
      radioButtons(
      "gensT",
      h3("What to show:"),
      choices = list(
        "Change over time" = "yes",
        "Specific generation" = "no"
      ),selected = "yes"),
     conditionalPanel(
       condition="input.gensT == 'no'",
       sliderInput(
         "gen_T",
         h3("Generation to show "),
         min = 0,
         max = 30000,
         value = 30000,
         step = 400
       )),
      radioButtons(
        "trait",
        h3("What trait to show:"),
        choices = list(
          "Cryptic female choice trait" = "female",
          "Sperm trait" = "male",
          "Sperm number"= "count"
        ),selected = "female"),
     radioButtons(
       "type_T",
       h3("Type of selection:"),
       choices = list(
         "Cryptic female choice" = "Cryptic female choice",
         "Stabilizing selection" = "Stabilizing"
       ),inline=T,selected = "Cryptic female choice"),
     conditionalPanel(
       condition="input.type_T == 'Cryptic female choice'",
     radioButtons(
       "Loc_T",
       h3("Number of Loci:"),
       choices = list(
         "2" = "2",
         "20" = "20"
       ),inline=T,selected = c("20")),
     radioButtons(
       "Popsize_T",
       h3("Population size:"),
       choices = list(
         "N = 1000" = "1000",
         "N = 10,000" = "10000"
       ),inline=T,selected = c("10000")),
     radioButtons(
       "V_T",
       h3("Starting Variation:"),
       choices = list(
         "Small" = "LV",
         "Medium" = "MV",
         "Large"="HV"
       ),inline=T,
       selected = c("MV")))
    ),
  mainPanel(plotOutput("runsTS",height=900))
),

# * 2.f Difference sensitivity -------------------------------------------
tabPanel(
  "Average trait deviations",
  titlePanel(
    h1(
      "Average trait deviations",
      align = "center"
    )),sidebarPanel(
      h3("Figure legend:"),"Evolutionary dynamics of deviations from predicted ejaculate investment (compared to analytical model) or sperm trait value (Compare to Figure 1B and Figure 3). You can either display how the mean across 50 populations changes over time (default) or display boxplots of all 50 populations at specific generations. In the default ('Change over time') The solid line represents the mean and the shade represents 95% confidence intervals in that mean. You can select which trait to show (ejaculate investment, or sperm trait). The defualt setting shows the base model assumptions, but you can change the parameters (number of loci, population size, and starting variation) to show the results of the sensitivity analyses. Values are only displayed every 400 generations. Dashed line represents zero or no deviation. Stabilizing selection runs are only available for 20 loci; moderate starting variation, and population sizes of 10,000.For cryptic female choice, deviation of sperm trait were calculated by subtracting mean sperm trait from mean cryptic female choice trait of a given population. For stabilizing selection, deviation of sperm trait was calculated by subtracting mean sperm trait from 50, the optimum set during those runs.",
      radioButtons(
        "gensD",
        h3("What to show:"),
        choices = list(
          "Change over time" = "yes",
          "Specific generation" = "no"
        ),selected = "yes"),
      conditionalPanel(
        condition="input.gensD == 'no'",
        sliderInput(
          "gen_D",
          h3("Generation to show "),
          min = 0,
          max = 30000,
          value = 30000,
          step = 400
        )),
      radioButtons(
        "trait_D",
        h3("What trait to show:"),
        choices = list(
          "Ejaculate investment"= "invest",
          "Sperm trait" = "male"
        ),selected = "invest"),
      radioButtons(
        "type_D",
        h3("Type of selection:"),
        choices = list(
          "Cryptic female choice" = "Cryptic female choice",
          "Stabilizing selection" = "Stabilizing"
        ),inline=T,selected = "Cryptic female choice"),
      conditionalPanel(
        condition="input.type_D == 'Cryptic female choice'",
        radioButtons(
          "Loc_D",
          h3("Number of Loci:"),
          choices = list(
            "2" = "2",
            "20" = "20"
          ),inline=T,selected = c("20")),
        radioButtons(
          "Popsize_D",
          h3("Population size:"),
          choices = list(
            "N = 1000" = "1000",
            "N = 10,000" = "10000"
          ),inline=T,selected = c("10000")),
        radioButtons(
          "V_D",
          h3("Starting Variation:"),
          choices = list(
            "Small" = "LV",
            "Medium" = "MV",
            "Large"="HV"
          ),inline=T,
          selected = c("MV")))
    ),
  mainPanel(plotOutput("runsDS",height=900))
),


# * 2.g Cors sensitivity ------------------------------------
tabPanel(
  "Correlations",
  titlePanel(
    h1(
      "Correlations",
      align = "center"
    )),sidebarPanel(
      h3("Figure legend:"),"Evolutionary dynamics of the mean genetic correlations of cryptic female choice trait and sperm trait. You can either display how these correlations change over time (default) or at specific generations display boxplots of all 50 populations (within-population correlations) and scatter plots of average male and female traits with line of best fit. In the default ('Change over time') genetic correlation graph the solid line represents the mean and the shade represents the standard error in that mean. The defualt setting shows the base model assumptions (i.e, Figure 2B), but you can change the parameters (number of loci, population size, and starting variation) to show the results of the sensitivity analyses. Stabilizing selection runs are only available for 20 loci; moderate starting variation, and population sizes of 10,000.",
      radioButtons(
        "gensC",
        h3("What to show:"),
        choices = list(
          "Change over time" = "yes",
          "Specific generation" = "no"
        ),selected = "yes"),
      conditionalPanel(
        condition="input.gensC == 'no'",
        sliderInput(
          "gen_C",
          h3("Generation to show "),
          min = 0,
          max = 30000,
          value = 30000,
          step = 400
        )),
      radioButtons(
        "type_C",
        h3("Type of selection:"),
        choices = list(
          "Cryptic female choice" = "Cryptic female choice",
          "Stabilizing selection" = "Stabilizing"
        ),inline=T,selected = "Cryptic female choice"),
      conditionalPanel(
        condition="input.type_C == 'Cryptic female choice'",
      radioButtons(
        "Loc_C",
        h3("Number of Loci:"),
        choices = list(
          "2" = "2",
          "20" = "20"
        ),inline=T,selected = c("20")),
      radioButtons(
        "Popsize_C",
        h3("Population size:"),
        choices = list(
          "N = 1000" = "1000",
          "N = 10,000" = "10000"
        ),inline=T,selected = c("10000")),
      radioButtons(
        "V_C",
        h3("Starting Variation:"),
        choices = list(
          "Small" = "LV",
          "Medium" = "MV",
          "Large"="HV"
        ),inline=T,
        selected = c("MV")))
    ),
  mainPanel(plotOutput("runsCS",height=900))
)
  )
)

# 3. Server side ----------------------------------------------------------
server <- function(input, output,session) {
  # Observes ----------------------------------------------------------------
  #Plotting runs
  observe({
    if(input$type_r=="Stabilizing"){
      updateSelectInput(session, "Loc_r",
                        selected = "20")
      updateSelectInput(session, "Popsize_r",
                        selected = "10000")
      updateSelectInput(session, "V_r",
                        selected = "MV")
    }})
  #average traits
  observe({
    if(input$type_T=="Stabilizing"){
      updateSelectInput(session, "Loc_T",
                        selected = "20")
      updateSelectInput(session, "Popsize_T",
                        selected = "10000")
      updateSelectInput(session, "V_T",
                        selected = "MV")
    }})
  
  #Devitonas
  observe({
    if(input$type_D=="Stabilizing"){
      updateSelectInput(session, "Loc_D",
                        selected = "20")
      updateSelectInput(session, "Popsize_D",
                        selected = "10000")
      updateSelectInput(session, "V_D",
                        selected = "MV")
    }})
  #Cors
  observe({
    if(input$type_C=="Stabilizing"){
      updateSelectInput(session, "Loc_C",
                        selected = "20")
      updateSelectInput(session, "Popsize_C",
                        selected = "10000")
      updateSelectInput(session, "V_C",
                        selected = "MV")
    }})  

# * 3.a Analtyical model -------------------------------------------------
#fun_plot_an
  output$fun_plot_an<-renderPlot({
    plot_function_an(1/input$alpha,input$beta)
  })
  
# * 3.b Selection function ----------------------------------------------
  output$fun_plot<-renderPlot({
    plot_function(input$fm,input$fs,input$cs,input$om)
  })
  

# * 3.c Plotting runs -----------------------------------------------------
  output$runs<-bindCache(renderPlot({
    runs_plot(input$means_r,input$alpha_r,input$all_r,input$populations,input$rsc_r,input$Loc_r,input$Popsize_r,input$Tr_r,input$V_r,input$type_r,run_data)
  }),input$means_r,input$alpha_r,input$all_r,input$populations,input$rsc_r,input$Loc_r,input$Popsize_r,input$Tr_r,input$V_r,input$type_r) 

# * 3.d Traits---------------------------------------------
  output$runsTS<-bindCache(renderPlot({
    if (input$gensT=="no"){#see if we are selecting specific generation
      gen<-input$gen_T
      if(input$gen_T==0){#generation zero in data is 1 so need to convert that
        gen<-1
      }

      if(input$trait=="male"){
        #mean plot        
        MT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
        run_data$Loci %in% input$Loc_T  &
          run_data$Var %in% input$V_T &
          run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
            x = A,
            y = MeanMale,
            color = Level,
            fill = Level
          )) +
      geom_boxplot(alpha = 0.5) +
      facet_grid(Tradeoff ~ .,
                 labeller = labeller(
                   Tradeoff = T.labs,
                   Var = Var.labs,
                   Loci = L.labs,
                   PopSize = P.labs
                 )) +
      scale_color_manual(name = "Risk of sperm competition", values =pal) +
      scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
      scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Strength of selection", y = "Mean male trait")
        if(input$type_T=="Stabilizing"){
          MT <-MT+geom_hline(aes(yintercept=50),linetype="dashed",alpha=0.5,color="black",size=2)
        }
    #standard deviation plot
    CVT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
                            run_data$Loci %in% input$Loc_T  &
                            run_data$Var %in% input$V_T &
                            run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
                              x = A,
                              y = CVMale,
                              color = Level,
                              fill = Level
                            )) +
      geom_boxplot(alpha = 0.5) +
      facet_grid(Tradeoff ~ .,
                 labeller = labeller(
                   Tradeoff = T.labs,
                   Var = Var.labs,
                   Loci = L.labs,
                   PopSize = P.labs
                 )) +
      scale_color_manual(name = "Risk of sperm competition", values =pal) +
      scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
      scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Strength of selection", y = "CV male trait") + theme(legend.position = "none")
    #put it all together
    return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }else if(input$trait=="female"){
        #mean plot
        MT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
                            run_data$Loci %in% input$Loc_T  &
                            run_data$Var %in% input$V_T &
                            run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
                              x = A,
                              y = MeanFemale,
                              color = Level,
                              fill = Level
                            )) +
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "Mean female trait")
        #standard deviation plot
        CVT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
                            run_data$Loci %in% input$Loc_T  &
                            run_data$Var %in% input$V_T &
                            run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
                              x = A,
                              y = CVFemale,
                              color = Level,
                              fill = Level
                            )) +
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "CV female trait") + theme(legend.position = "none")
        #put it all together
        return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }else {
        #mean plot
        MT<-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
                              run_data$Loci %in% input$Loc_T  &
                              run_data$Var %in% input$V_T &
                              run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
                                x = A,
                                y = MeanCount,
                                color = Level,
                                fill = Level
                              )) +
          geom_hline(data=PredAn,aes(yintercept =PredN,color=Level), linetype = "dashed",show.legend = T,alpha=0.5)+
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "Mean sperm number")
        #CV plot
        CVT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_T &
                            run_data$Loci %in% input$Loc_T &
                            run_data$Var %in% input$V_T &
                            run_data$Generation == gen& run_data$type %in% input$type_T, ], aes(
                              x = A,
                              y = CVCount,
                              color = Level,
                              fill = Level
                            )) +
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "CV sperm number") + theme(legend.position = "none")
        #put it all together
        return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }
      
    }else{#not showing specific generations
      if(input$trait=="male"){
        #mean plot
        MT<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                             SumRuns$Loci %in% input$Loc_T &
                             SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                   aes(
                     x = Generation,
                     y = MeanMale_mean,
                     color = Level,
                     fill = Level
                   )) +
          geom_line() +
          geom_ribbon(
            aes(
              ymin = MeanMale_mean - MeanMale_sd / sqrt(50),
              ymax = MeanMale_mean + MeanMale_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~ A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values = pal) +
          ylab("Mean male trait") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        if(input$type_T=="Stabilizing"){
          MT <-MT+geom_hline(aes(yintercept=50),linetype="dashed",alpha=0.5,color="black",size=2)
        }
        #standard deviation plot
        CVT <- ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                           SumRuns$Loci %in% input$Loc_T &
                           SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                 aes(
                   x = Generation,
                   y = CVMale_mean,
                   color = Level,
                   fill = Level
                 )) +
          geom_line() +
          geom_ribbon(
            aes(
              ymin = CVMale_mean - CVMale_sd / sqrt(50),
              ymax = CVMale_mean + CVMale_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~ A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values = pal) +
          scale_fill_manual(name = "Risk of sperm competition", values = pal) +
          ylab("CV male trait") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position = "none")
        #put it all together
      return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }else if(input$trait=="female"){
        #mean plot
        MT<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                             SumRuns$Loci %in% input$Loc_T &
                             SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                   aes(
                     x = Generation,
                     y = MeanFemale_mean,
                     color = Level,
                     fill = Level
                   )) +
          geom_line() +
          geom_ribbon(
            aes(
              ymin = MeanFemale_mean - MeanFemale_sd / sqrt(50),
              ymax = MeanFemale_mean + MeanFemale_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
           Tradeoff ~  A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values = pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean Female trait") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        #Standard deviation plot
        CVT <-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                           SumRuns$Loci %in% input$Loc_T &
                           SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                 aes(
                   x = Generation,
                   y = CVFemale_mean,
                   color = Level,
                   fill = Level
                 )) +
          geom_line() +
          geom_ribbon(
            aes(
              ymin = CVFemale_mean - CVFemale_sd / sqrt(50),
              ymax = CVFemale_mean + CVFemale_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~  A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("CV female trait") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position = "none")
        #put it all together
        return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }else {#sperm number
        #mean plot
        MT<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                             SumRuns$Loci %in% input$Loc_T &
                             SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                   aes(
                     x = Generation,
                     y = MeanCount_mean,
                     color = Level,
                     fill = Level
                   )) +
          geom_hline(data=PredAn,aes(yintercept =PredN,color=Level), linetype = "dashed",show.legend = T,alpha=0.5)+
          geom_line() +
          geom_ribbon(
            aes(
              ymin = MeanCount_mean - MeanCount_sd / sqrt(50),
              ymax = MeanCount_mean + MeanCount_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
           Tradeoff ~  A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean sperm number") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        #Standard deviation plot
        CVT <-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_T &
                           SumRuns$Loci %in% input$Loc_T &
                           SumRuns$Var %in% input$V_T& SumRuns$type %in% input$type_T, ],
                 aes(
                   x = Generation,
                   y = CVCount_mean,
                   color = Level,
                   fill = Level
                 )) +
          geom_line() +
          geom_ribbon(
            aes(
              ymin = CVCount_mean - CVCount_sd / sqrt(50),
              ymax = CVCount_mean + CVCount_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~ A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("CV sperm number") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position = "none")
        #put it all together
        return((MT/CVT&theme(legend.position ="bottom"))+plot_layout(guides="collect"))
      }
    }
  }),input$gensT,input$trait,input$Popsize_T ,input$Loc_T,input$V_T,input$gen_T,input$type_T) 
  
# * 3.e Devs---------------------------------------------
  output$runsDS<-bindCache(renderPlot({
    if (input$gensD=="no"){#see if we are selecting specific generation
      gen<-input$gen_D
      if(input$gen_D==0){#generation zero in data is 1 so need to convert that
        gen<-1
      }
      
      if(input$trait_D=="male"){
        #mean plot        
        MT <-ggplot(run_data[run_data$PopSize %in% input$Popsize_D &
                               run_data$Loci %in% input$Loc_D  &
                               run_data$Var %in% input$V_D &
                               run_data$Generation == gen& run_data$type %in% input$type_D, ], aes(
                                 x = A,
                                 y = DiffOpt,
                                 color = Level,
                                 fill = Level
                               ))  +
          geom_hline(aes(yintercept =0),color="black", linetype = "dashed",show.legend = T,alpha=0.5,size=2)+
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Difference from female trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "Mean male trait")
        if(input$type_T=="Stabilizing"){
          MT <-MT+ylab("Difference from 50")
        }
        #standard deviation plot
        #put it all together
        return((MT&theme(legend.position ="bottom")))
      }else {
        #mean plot
        MT<-ggplot(run_data[run_data$PopSize %in% input$Popsize_D &
                              run_data$Loci %in% input$Loc_D  &
                              run_data$Var %in% input$V_D &
                              run_data$Generation == gen& run_data$type %in% input$type_D, ], aes(
                                x = A,
                                y = relDif,
                                color = Level,
                                fill = Level
                              )) +
          geom_hline(aes(yintercept =0),color="black", linetype = "dashed",show.legend = T,alpha=0.5,size=2)+
          geom_boxplot(alpha = 0.5) +
          facet_grid(Tradeoff ~ .,
                     labeller = labeller(
                       Tradeoff = T.labs,
                       Var = Var.labs,
                       Loci = L.labs,
                       PopSize = P.labs
                     )) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
          scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Strength of selection", y = "Relative deviation \n from predicted investment")
        #put it all together
        return((MT&theme(legend.position ="bottom")))
      }
      
    }else{#not showing specific generations
      if(input$trait_D=="male"){
        #mean plot
        MT<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_D &
                             SumRuns$Loci %in% input$Loc_D&
                             SumRuns$Var %in% input$V_D& SumRuns$type %in% input$type_D, ],
                   aes(
                     x = Generation,
                     y = DiffOpt_mean,
                     color = Level,
                     fill = Level
                   )) +
          geom_hline(aes(yintercept =0),color="black", linetype = "dashed",show.legend = T,alpha=0.5,size=2)+
          geom_line() +
          geom_ribbon(
            aes(
              ymin = DiffOpt_mean - DiffOpt_sd / sqrt(50),
              ymax = DiffOpt_mean + DiffOpt_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~ A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values = pal) +
          ylab("Difference from female trait") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        if(input$type_T=="Stabilizing"){
          MT <-MT+ylab("Difference from 50")
        }
        #standard deviation plot
        #put it all together
        return((MT&theme(legend.position ="bottom")))
      }else {#sperm number
        #mean plot
        MT<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_D &
                             SumRuns$Loci %in% input$Loc_D &
                             SumRuns$Var %in% input$V_D& SumRuns$type %in% input$type_D, ],
                   aes(
                     x = Generation,
                     y = relDif_mean,
                     color = Level,
                     fill = Level
                   ))  +
          geom_hline(aes(yintercept =0),color="black", linetype = "dashed",show.legend = T,alpha=0.5,size=2)+
          geom_line() +
          geom_ribbon(
            aes(
              ymin = relDif_mean - relDif_sd / sqrt(50),
              ymax = relDif_mean + relDif_sd / sqrt(50)
            ),
            alpha = 0.3,
            color = NA
          ) +
          facet_nested(
            Tradeoff ~  A,
            labeller = labeller(
              Tradeoff = T.labs,
              A = A.labs
            )
          ) +
          scale_color_manual(name = "Risk of sperm competition", values =pal) +
          scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Relative deviation \n from predicted investment") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

        #put it all together
        return((MT&theme(legend.position ="bottom")))
      }
    }
  }),input$gensD,input$trait_D,input$Popsize_D ,input$Loc_D,input$V_D,input$gen_D,input$type_D) 
  
# * 3.f Cors ---------------------------------------------
  output$runsCS<-bindCache(renderPlot({
    if (input$gensC=="no"){#specific generation
      gen<-input$gen_C
      if(input$gen_C==0){#if zero generation (aka start is selected convert to 1)
        gen<-1
      }
      #within-population correlation
      wp<-ggplot(run_data[run_data$PopSize %in% input$Popsize_C &
                            run_data$Loci %in% input$Loc_C &
                            run_data$Var %in% input$V_C &
                            run_data$Generation == gen&
                            run_data$type %in% input$type_C, ], aes(
                              x = A,
                              y = cor,
                              color = Level,
                              fill = Level
                            )) +
        geom_boxplot(alpha = 0.5) +
        facet_grid(Tradeoff ~ .,
                   labeller = labeller(
                     Tradeoff = T.labs,
                     Var = Var.labs,
                     Loci = L.labs,
                     PopSize = P.labs
                   )) +
        scale_color_manual(name = "Risk of sperm competition", values =pal) +
        scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Mean male trait") +
        scale_x_discrete(labels = c("Weak", "Moderate", "Strong")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Strength of selection", y = "Within population correlation \n between male and female traits") +
        geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed")
      #Across-population correlations
      ap<-ggplot(run_data[run_data$PopSize %in% input$Popsize_C &
                            run_data$Loci %in% input$Loc_C &
                            run_data$Var %in% input$V_C &
                            run_data$Generation == gen&
                            run_data$type %in% input$type_C, ],
                 aes(
                   x = MeanFemale,
                   y = MeanMale,
                   color = Level,
                   fill = Level
                 )) +
        geom_point(show.legend = F) +
        geom_smooth(method = "lm", show.legend = F) +
        facet_grid(
          Tradeoff ~ A,
          labeller = labeller(
            Tradeoff = T.labs,
            Var = Var.labs,
            Loci = L.labs,
            PopSize = P.labs,
            A = A.labs
          )
        ) +
        scale_color_manual(name = "Risk of Sperm Competition", values =pal) +
        scale_fill_manual(name = "Risk of Sperm Competition", values =pal) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x ="Mean female", y = "Mean male")
      #put the plots together
      return((wp / ap &theme(legend.position = "bottom")) + plot_layout(guides = "collect"))
      
    }else{#show change over time
      #within-population correlation
      wp<-ggplot(SumRuns[SumRuns$PopSize %in% input$Popsize_C  &
                           SumRuns$Loci %in% input$Loc_C &
                           SumRuns$Var %in% input$V_C&
                           SumRuns$type %in% input$type_C, ],
                 aes(
                   x = Generation,
                   y = cor_mean,
                   color = Level,
                   fill = Level
                 )) +
        geom_line() +
        geom_ribbon(
          aes(
            ymin = cor_mean - cor_sd / sqrt(50),
            ymax = cor_mean + cor_sd / sqrt(50)
          ),
          alpha = 0.3,
          color = NA
        ) +
        facet_nested(
          Tradeoff ~  A,
          labeller = labeller(
            Tradeoff = T.labs,
            Var = Var.labs,
            Loci = L.labs,
            PopSize = P.labs,
            A = A.labs
          )
        ) +
        scale_color_manual(name = "Risk of sperm competition", values =pal) +
        scale_fill_manual(name = "Risk of sperm competition", values =pal) + ylab("Within-population correlation \n between male and female traits") +
        geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
        scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = guide_legend(show = FALSE))

      #put the plots together
        return((wp)+plot_layout(guides="collect"))
      
    }
  }),input$gensC,input$Popsize_C ,input$Loc_C,input$V_C,input$gen_C,input$type_C) 
  
  }

# Run app ----
shinyApp(ui, server)
