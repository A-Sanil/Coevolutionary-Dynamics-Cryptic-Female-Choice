#R script to make all the main figures and SI figures 4-11
#for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#please send any questions to mkustra@ucsc.edu
# 1. Loading up libraries -------------------------------------------------
library(tidyverse)
library(patchwork)
library(furrr)
library(parallel)
library(data.table)
library(ggh4x)
library(viridis)
library(forcats)
library(RColorBrewer)
library(wesanderson)
# 2. Default plotting themes/labels ---------------------------------------
mytheme <-
  theme_bw() + theme(
    legend.position = "bottom",
    #this puts legend on the bottom
    axis.title = (element_text(color = "black")),
    #Makes the axis line black and  thicker
    text = element_text(size = 16),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", size = 1)
  )
theme_set(mytheme)
##Labels
#vector of labels for strength of selection
A.labs<-c("Weak selection","Moderate selection","Strong selection")
#name the vector to match actual values in data frame
names(A.labs)<-c("50.0","12.5","1.0")
#vector of labels for sex traits
Sex.labs<-c("Cryptic choice trait","Sperm trait")
#name the vector to match actual values in data frame
names(Sex.labs)<-c("Female","Male")
#vector of labels for variation
Var.labs <- c("Small Var", "Medium Var", "Large Var")
#name the vector to match actual values in data frame
names(Var.labs) <- c("LV", "MV", "HV")
#vector of labels for tradeoff
T.labs <- c("NA", "No tradeoff", "Tradeoff")
#name the vector to match actual values in data frame
names(T.labs) <- c("-1", "false", "true")
#vector of labels for number of loci
L.labs <- c("2 Loci", "20 Loci")
#name the vector to match actual values in data frame
names(L.labs) <- c("2", "20")
#vector of labels for population size
P.labs <- c("N = 1000", "N = 10000")
#name the vector to match actual values in data frame
names(P.labs) <- c("1000", "10000")

#Make color palette
pal <- rev(plasma(8)[c(1, 3, 5, 7)])

# 3. Loading up/processing data -------------------------------------------

runs_rawS<-readRDS("Data/All_runs_shortened.rds")
#Mutate runs_raw S to make sure all factors in the proper order
runs_rawS <- runs_rawS %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0"), exclude = NULL),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000")),
    CVCount=SDCount/MeanCount,
    CVMale=SDMale/MeanMale,
    CVFemale=SDFemale/MeanFemale)
#Summarize Runs for phenotypic plots
#Takes a long time to run so optionally load up already processed file
#SumRuns <- runs_rawS %>%
#  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize) %>%
#  summarise_all(funs(mean = mean,sd = sd), na.rm = T) %>%
#  mutate(
#    A = factor(A, levels = c("50.0","12.5","1.0")),
#    Var = factor(Var, levels = c("LV", "MV", "HV")),
#    PopSize = factor(PopSize, levels = c("1000", "10000"))
#  )
#saveRDS(SumRuns,"SumRunsAll.rds")
SumRuns<-readRDS("Data/SumRunsAll.rds")%>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )
#Stabilizing selection
runs_rawStab<-readRDS("Data/Compiled_Shortened_DataS_MiddleVar_L_20.rds")
#Mutate runs_raw S to make sure all factors in the proper order
runs_rawStab <- runs_rawStab %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0"), exclude = NULL),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000")),
    CVCount=SDCount/MeanCount,
    CVMale=SDMale/MeanMale,
    CVFemale=SDFemale/MeanFemale)

#Summarize Runs for phenotypic plots
SumRunsStab <- runs_rawStab %>%
  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize) %>%
  summarise_all(funs(mean = mean,sd = sd), na.rm = T) %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )
#filter for the default scenario (20 number of loci, large population size, medium variation, and not fair raffle, )
SumRunsF20 <- SumRuns %>%
  filter(Loci == "20", PopSize == "10000", Var == "MV")
#filter last generation for inset plots
runs_rawS_LGF20 <- runs_rawS %>%
  filter(Loci == "20",
         PopSize == "10000",
         Var == "MV",
         Generation == 30000)

#Stabilizing
#filter for the default scenario (20 number of loci, large population size, medium variation, and not fair raffle, )
SumRunsF20Stab <- SumRunsStab %>%
  filter(Loci == "20", PopSize == "10000", Var == "MV")
#filter last generation for inset plots
runs_rawS_LGF20Stab <- runs_rawStab %>%
  filter(Loci == "20",
         PopSize == "10000",
         Var == "MV",
         Generation == 30000)


#make data frame of predicted investment and sperm number
#numbers directly taking from "AdaptiveDynamics.jl"
PredAn<-data.frame(Level=c("0.25","0.5","0.75","1.0"),Tradeoff=c("false","false","false","false","true","true","true","true"),PredN=c(16.10988593067346,23.66053786200191,28.9673742726416,33.18317009243252),PredT=c(16.10988593067346,23.66053786200191,28.9673742726416,33.18317009243252,805.4942965336729,1183.0268931000955,1448.36871363208,1659.1585046216262),PredS=c(16.10988593067346,23.66053786200191,28.9673742726416,33.18317009243252))

#Reading in last_2000 generation summary for all runs 
runs_rawSumLast<-readRDS("Data/All_sum_Last2000.rds")

#Reading in data from raw last 20000 data for main figures
LData<-readRDS("Data/Last2000_MiddleVar_L_20.rds")

#Reading in data for evolutionary lag across all sensitivity analyses
LagAll<-read_rds("Data/Lag_2000_all.rds")%>%
  mutate(Var = factor(Var, levels = c("LV", "MV", "HV")),
         Diff= (Lag0-MinDiff)/Lag0)

# 4. Figure 1. Traits mean and cv ------------------------------------------
#please note that facets were styalized manually in illustrator. 
###Male (sperm trait and number) and Female (cryptic choice trait) traits
runs_rawS_LGF20$CVCount<-runs_rawS_LGF20$SDCount/runs_rawS_LGF20$MeanCount
runs_rawS_LGF20L<-runs_rawS_LGF20%>%
  ungroup()%>%
  select(MeanMale,MeanFemale,MeanCount,CVMale,CVFemale,CVCount,Level,A,Tradeoff)%>%
  rename(
    Mean_Sperm = MeanMale,
    `Mean_Cryptic Choice` = MeanFemale,
    CV_Sperm = CVMale,
    `CV_Cryptic Choice` = CVFemale,
    Mean_Count=MeanCount,
    CV_Count=CVCount
  ) %>%
  pivot_longer(
    !c(Level,A,Tradeoff),
    names_to = c(".value", "Sex"),
    names_sep = "_"
  )%>%
  mutate(Sex=factor(Sex,levels=c("Cryptic Choice","Sperm","Count")))

#vector of labels for strength of selection
A.labs2<-c("Weak preference","Moderate preference","Strong preference")
#name the vector to match actual values in data frame
names(A.labs2)<-c("50.0","12.5","1.0")
(TM<-ggplot(runs_rawS_LGF20L, aes(
  x = Sex,
  y = Mean,
  fill = Sex,
  color = Sex,
  shape=Sex
)) +
    facet_nested(Tradeoff ~ A+Level, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs2))+
    ylab("Mean trait")+
    xlab("Risk of sperm competition (q)")+
    geom_point(alpha=0.3,show.legend = T,position=position_jitterdodge(dodge.width = 1),size=0.8)+
    geom_boxplot(alpha = 0.6,show.legend = F,outlier.shape = NA,width=0.5,position = position_dodge(1)) +
    scale_color_manual(name = "Trait", values = c("#F8766D","#5F9EA080","black")) +
    scale_fill_manual(name = "Trait", values = c("#F8766D","#5F9EA080","black"))+
    labs(shape="Trait")+guides(fill = guide_legend(override.aes = list(
      size = 5, alpha = 1
    )))+theme(legend.position = "top"))


(TCV<-ggplot(runs_rawS_LGF20L, aes(
  x = Sex,
  y = CV,
  fill = Sex,
  color = Sex,
  shape = Sex
)) +
    facet_nested(Tradeoff ~ A + Level, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs2))+
    ylab("Trait coefficient of variation")+
    xlab("Risk of sperm competition (q)")+
    geom_point(alpha=0.3,show.legend = F,position=position_jitterdodge(dodge.width = 1),size=0.8)+
    geom_boxplot(alpha = 0.6,show.legend = F,outlier.shape = NA,width=0.5,position = position_dodge(1)) +
    scale_color_manual(name = "Trait", values = c("#F8766D","#5F9EA080","black")) +
    scale_fill_manual(name = "Trait", values = c("#F8766D","#5F9EA080","black"))+
    labs(shape="Trait"))



((TM/TCV) & theme(text=element_text(size=18),axis.ticks.x = element_blank(),axis.text.x = element_blank())) + plot_annotation(tag_levels = "A", tag_suffix = ")")

#ggsave(
#  "Figure1.pdf",
#  height = 183 * 1.5,
#  width = 183 * 2,
#  units = "mm"
#)

# 5. Figure 2. Investment ----------------------------------------------------
runs_rawS_LGF20_pred<-merge(runs_rawS_LGF20,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=MeanMale-MeanFemale,investScale=ifelse(Tradeoff=="true",investment/50,investment))

#stabilizing selection
runs_rawS_LGF20_predStab<-merge(runs_rawS_LGF20Stab,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=MeanMale-50,investScale=ifelse(Tradeoff=="true",investment/50,investment))

#Merge everything
runs_rawS_LGF20_predStab$Select<-"Sperm competition only"
runs_rawS_LGF20_pred$Select<-"Cryptic female choice"
runs_rawS_LGF20_pred<-rbind(runs_rawS_LGF20_pred,runs_rawS_LGF20_predStab)


#difference in investment
(BP_Diff<- ggplot(runs_rawS_LGF20_pred, aes(
  x = Level,
  y = relDif,
  fill = Level,
  color = Level
)) +
    facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs))+
    ylab("Relative deviation \n from predicted investment")+
    xlab("Risk of sperm competition (q)")+
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
    geom_jitter(alpha=0.5,height=0,show.legend = F)+
    geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_continuous(n.breaks = 4,limits=c(-0.5,.25)))

#correlation
(CorCM<-ggplot(runs_rawS_LGF20_pred, aes(
  x = MeanMale,
  y = MeanCount,
  fill = Level,
  color = Level
)) +
    geom_hline(data=PredAn,aes(yintercept =PredN,color=Level), linetype = "dashed",show.legend = T,alpha=0.5)+
    geom_point(show.legend=F)+
    geom_smooth(show.legend = F)+
    facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs),scales="free_y")+
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal)+
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+
    labs(y="Sperm number",x="Sperm trait")+
    scale_x_continuous(n.breaks = 4,limits=c(20,160)))


#putting it all together
(((CorCM/BP_Diff)) &
    theme(panel.spacing.y = unit(0.5, "lines")) &
    guides(color = guide_legend(override.aes = list(
      size = 5, alpha = 1
    )), color = "none")
) + plot_annotation(tag_levels = "A", tag_suffix = ")") + plot_layout(guides = "collect")

#ggsave(
#  "Figure2.pdf",
#  height = 183 * 2,
#  width = 183 * 1.8,
#  units = "mm"
#)

# 6. Figure 3. Genetic Correlations----------------------------
#Within population genetic correlation
(wp20 <-
    ggplot(SumRunsF20,
           aes(
             x = Generation,
             y = cor_mean,
             color = Level,
             fill = Level
           )) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = cor_mean - cor_sd ,
        ymax = cor_mean + cor_sd
      ),
      alpha = 0.3,
      color = NA
    ) +
    facet_nested(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A = A.labs)) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal) +
    ylab("Genetic correlation between f and m") + 
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed") + scale_x_continuous(n.breaks = 5)+ scale_y_continuous(n.breaks = 4,limits=c(-0.01,0.3)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

#Example coevolution of phenotypic traits

#first filtering out the two extremes of stron, weak selection and low + high risk of sperm competition and main focus for main figures (Medium variation, 20 loci, pop size of 10,000)
runs_raw_p<-runs_rawS %>% filter(A%in%c("50.0","1.0"),
                                 Level %in%c("0.25","1.0"),
                                 Var == "MV",
                                 Loci == "20",
                                 PopSize == "10000")

#Select example populations that are minimum and maximum values
Selections<-runs_raw_p%>%group_by(A,Level,Tradeoff)%>%
  mutate(Max=max(MeanFemale),Min=min(MeanFemale),RepMin=Rep[which.min(MeanFemale)][1],RepMax=Rep[which.max(MeanFemale)][1])%>%
  mutate(Rep2=ifelse(Rep==RepMin,"Min",ifelse(Rep==RepMax,"Max",Rep)))

#want minimum and maximum population, plus two other random populations for the figure
vals<-c("Max","Min",10,20)

#make the evolutionary path plot

(zW <-Selections%>%filter(Rep2%in%vals,Tradeoff=="true")%>%
    arrange(Generation) %>%
    ggplot(aes(x = MeanFemale, y = MeanMale, color = factor(Rep2)))+
    geom_abline(slope=1,intercept = 0,color="black",alpha=0.3,linetype="dashed") +
    geom_path(aes(group = Rep2), alpha = 0.4) +
    geom_point(
      data = Selections %>% filter(Rep2%in%vals,
                                   Generation %in% c(30000),Tradeoff=="true"
      ),size=2
    ) +
    geom_point(
      data = Selections %>% filter(Rep2%in%vals,
                                   Generation %in% c(1),Tradeoff=="true",
      ),
      color = "black",
      shape = 15,size=2
    ) +
    scale_color_manual(values = wes_palette(n=4, name="GrandBudapest2"),breaks = c("Min","Max","10","20")) + theme(legend.position = "none") +
    labs(y = "Sperm trait (m)", x = "Cryptic choice trait (f)")+
    facet_nested(Level ~ A, labeller = labeller(A = A.labs),scales="free",independent = "all")
)


((wp20+theme(legend.position = c(0.25,0.75))|zW) & theme(text=element_text(size=18))&
    guides(fill = guide_legend(override.aes = list(
      size = 5, alpha = 1
    )), color = "none")
) + plot_annotation(tag_levels = "A", tag_suffix = ")")

#ggsave(
#  "Figure3.pdf",
#  height = 183,
#  width = 183 * 2,
#  units = "mm"
#)

# 7. Figure 4. Trait deviation + gamma----------------------------------------
(BP_OptDiff<- ggplot(runs_rawS_LGF20_pred, aes(
  x = Level,
  y = DiffOpt,
  fill = Level,
  color = Level
)) +
  facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs))+
  ylab("Difference from optimum sperm trait (m)")+
  xlab("Risk of sperm competition (q)")+
  geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
  geom_jitter(alpha=0.5,height=0,show.legend = F)+
  geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
  scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
  scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_continuous(n.breaks = 4,limits=c(-2,1)))

# * 7.b. Figure 4b. Gamma estimations ------------------------------------------------------------------
(BP_Gamma<- ggplot(runs_rawS_LGF20_pred, aes(
  x = Level,
  y = GMale,
  fill = Level,
  color = Level
)) +
  facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs))+
  ylab("Gamma sperm trait (m)")+
  xlab("Risk of sperm competition (q)")+
  geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
  geom_jitter(alpha=0.5,height=0,show.legend = F)+
  geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
  scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
  scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))

(((BP_OptDiff/BP_Gamma)) &
     theme(panel.spacing.y = unit(0.5, "lines")) &
     guides(color = guide_legend(override.aes = list(
       size = 5, alpha = 1
     )), color = "none")
  ) + plot_annotation(tag_levels = "A", tag_suffix = ")") + plot_layout(guides = "collect")
#ggsave(
#  "Figure4.pdf",
#  height = 183 * 2,
#    width = 183 * 1.8,
#    units = "mm"
#  )
# 8. Figure 5. Evolutionary Lags ---------------------------------------------------------
# * 8.a. Lag calculation function -----------------------------------------
#Difference function that was used to find maximum lags
#Not actually used in this code for computational time.
#To just replicate figures processed data is already loaded in the MaxAll data
Diff<-function(x,y,lag){
  Data<-rep(-1,lag)
  for(i in 0:lag){
    y2<-(lead(y,i)-mean(lead(y,i),na.rm=T))/sd(lead(y,i),na.rm=T)
    x<-(x-mean(x))/sd(x,na.rm=T)
    Data[i+1]<-mean(abs(x-y2),na.rm=T)
  }
  #return(Data)
  return(tibble(MinLag=which.min(Data)-1,MinDiff=min(Data),Lag0=Data[1]))
}

# * 8.b. Lag plots --------------------------------------------------------
#For the main figure of interest need to isolate loci = 20; medium variation, and popsize =10000
Max<-LagAll%>%filter(Loci=="20",Var=="MV",PopSize=="10000")

#vector of labels for strength of selection
A.labs2<-c("Weak","Moderate","Strong")
#name the vector to match actual values in data frame
names(A.labs2)<-c("50.0","12.5","1.0")

#lag that minimized absolute difference.
(BP_MinL<- ggplot(Max, aes(
  x = Level,
  y = MinLag,
  fill = Level,
  color = Level
)) +
    facet_grid(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs2))+
    ylab("Lag with minium difference (generation)")+
    xlab("Risk of sperm competition (q)")+
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
    geom_jitter(alpha=0.5,height=0,show.legend = F)+
    geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
    scale_color_manual(name = "Risk of sperm competition", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))

#Difference in best fit lag time from lag
(BP_Diff_opt<- ggplot(Max, aes(
  x = Level,
  y = Diff,
  fill = Level,
  color = Level
)) +
    facet_grid(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs2))+
    ylab("Relative improvement in MAE with lag")+
    xlab("Risk of sperm competition (q)")+
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
    geom_jitter(alpha=0.5,height=0,show.legend = F)+
    geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))

#Filtering data for example time series and plotting it
(LagP<-LData%>%
  filter(Rep==1,Tradeoff=="false",A=="1.0",Level=="0.25",Generation>29000,Generation <29800)%>%
  ggplot( aes(
    x = Generation,
    y = MeanFemale,
  )) +
  geom_line(color="#F8766D",size=1)+
  ylab("Mean trait")+
  xlab("Generation")+
  geom_line(aes(y=MeanMale), color="#5F9EA080",size=1,linetype="dashed"))

((BP_MinL)|BP_Diff_opt|(LagP))+ plot_annotation(tag_levels = "A", tag_suffix = ")")

#please note that arrow labels and lag annotations for 5C were added manually
#ggsave(
#  "Figure5.pdf",
#  height = 130,
#  width = 180*2,
#  units = "mm"
#)

# 9. SI Figure 5 Fair raffle comparison -----------------------------------
#Load up fair raffle results
runs_rawSfr<-read_rds("Data/Compiled_Data_Shortened_FairRaffle.rds")
runs_rawSfr_LGF20<-runs_rawSfr %>%
  filter(Generation == 30000)
#merge predictions to FR
runs_rawS_LGF20_predFR<-merge(runs_rawSfr_LGF20,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=MeanMale-MeanFemale,investScale=ifelse(Tradeoff=="true",investment/50,investment))
#Change selection naming for this graph
runs_rawS_LGF20_predStab$Select<-"Sperm competition only (stabilizing selection)"
runs_rawS_LGF20_predFR$Select<-"Fair raffle"

#Merge them together
runs_rawS_LGF20_predFR<-rbind(runs_rawS_LGF20_predFR%>%select(intersect(colnames(runs_rawS_LGF20_predFR),colnames(runs_rawS_LGF20_predStab))),runs_rawS_LGF20_predStab%>%select(intersect(colnames(runs_rawS_LGF20_predFR),colnames(runs_rawS_LGF20_predStab))))

#New labels for no selection
A.labs3<-c("No selection","Weak selection","Moderate selection","Strong selection")
#name the vector to match actual values in data frame
names(A.labs3)<-c("None","50.0","12.5","1.0")

(BP_DiffFR<- ggplot(runs_rawS_LGF20_predFR%>%mutate(A=factor(A,levels=c("None","50.0","12.5","1.0")))%>%filter(Tradeoff=="false"), aes(
  x = Level,
  y = relDif,
  fill = Level,
  color = Level
))+
    ylab("Relative deviation from predicted investment")+
    xlab("Risk of sperm competition (q)")+
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
    facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs3))+
    geom_jitter(alpha=0.5,height=0,show.legend = F)+
    geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
    scale_color_manual(name = "Risk of sperm competition", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))

#ggsave(
#  "FigureS5.pdf",
#  height = 183 ,
#  width = 183 * 1.2,
#  units = "mm"
#)
# 10. SI Figure 6: Directional selection ----------------------------------
#Read in directional selection simulation data
runs_rawSD<-read_rds("Data/Compiled_Shortened_Data_D.rds")
runs_rawSD <- runs_rawSD %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0"), exclude = NULL),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000")),
    CVCount=SDCount/MeanCount,
    CVMale=SDMale/MeanMale,
    CVFemale=SDFemale/MeanFemale)

runs_rawSD <- runs_rawSD %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0"), exclude = NULL),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000")),
    CVCount=SDCount/MeanCount,
    CVMale=SDMale/MeanMale,
    CVFemale=SDFemale/MeanFemale)

#Put in stabilizing selection
runs_rawStab<-readRDS("Data/Compiled_Shortened_DataS_MiddleVar_L_20.rds")
runs_rawS_LGF20D <- runs_rawSD %>%
  filter(Loci == "20",
         PopSize == "10000",
         Var == "MV",
         Generation == 30000)
runs_rawS_LGF20Stab <- runs_rawStab %>%
  filter(Loci == "20",
         PopSize == "10000",
         Var == "MV",
         Generation == 30000)

runs_rawS_LGF20_predD<-merge(runs_rawS_LGF20D,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=MeanMale-MeanFemale,investScale=ifelse(Tradeoff=="true",investment/50,investment))

#stabilizing selection
runs_rawS_LGF20_predStab<-merge(runs_rawS_LGF20Stab,PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT,DiffOpt=MeanMale-50,investScale=ifelse(Tradeoff=="true",investment/50,investment))

#Merge everything
runs_rawS_LGF20_predStab$Select<-"Sperm competition only"
runs_rawS_LGF20_predD$Select<-"Cryptic female choice"
runs_rawS_LGF20_predD<-rbind(runs_rawS_LGF20_predD%>%select(intersect(colnames(runs_rawS_LGF20_predD),colnames(runs_rawS_LGF20_predStab))),runs_rawS_LGF20_predStab%>%select(intersect(colnames(runs_rawS_LGF20_predD),colnames(runs_rawS_LGF20_predStab))))


#difference in investment
(BP_DiffD<- ggplot(runs_rawS_LGF20_predD, aes(
  x = Level,
  y = relDif,
  fill = Level,
  color = Level
)) +
    facet_nested(Tradeoff ~ Select+A, labeller = labeller(Tradeoff = T.labs, A =                                                 A.labs))+
    ylab("Relative deviation from predicted investment")+
    xlab("Risk of sperm competition (q)")+
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed")+
    geom_jitter(alpha=0.5,height=0,show.legend = F)+
    geom_boxplot(alpha = 0.3,show.legend = F,outlier.shape = NA,width=0.5) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_continuous(n.breaks = 4,limits=c(-0.5,.25)))

#Genetic correlations
SumRunsF20D <- runs_rawSD %>%
  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize) %>%
  summarise_all(funs(mean = mean,sd = sd), na.rm = T) %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )

(wp20D <-
    ggplot(SumRunsF20D,
           aes(
             x = Generation,
             y = cor_mean,
             color = Level,
             fill = Level
           )) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = cor_mean - cor_sd ,
        ymax = cor_mean + cor_sd
      ),
      alpha = 0.3,
      color = NA
    ) +
    facet_nested(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A = A.labs)) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal) +
    ylab("Genetic correlation between f and m") + 
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed") + scale_x_continuous(n.breaks = 5)+ scale_y_continuous(n.breaks = 4,limits=c(-0.01,0.3)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

(((BP_DiffD/wp20D)) &
     theme(panel.spacing.y = unit(0.5, "lines")) &
     guides(color = guide_legend(override.aes = list(
       size = 5, alpha = 1
     )), color = "none")
  ) + plot_annotation(tag_levels = "A", tag_suffix = ")") + plot_layout(guides = "collect")

#ggsave(
#  "FigureS6.png",
#  height = 183 * 2,
#  width = 183 * 1.8,
#  units = "mm"
#)

# * 10.b.  During directional selection -----------------------------------


runs_rawSD_F1000<-read_rds("Data/Compiled_First1000_D.rds")

SumRunsF1000 <- runs_rawSD_F1000 %>%
  group_by(Generation, Level, A, Tradeoff, Loci, Var, PopSize) %>%
  summarise_all(funs(mean = mean,sd = sd), na.rm = T) %>%
  mutate(
    A = factor(A, levels = c("50.0","12.5","1.0")),
    Var = factor(Var, levels = c("LV", "MV", "HV")),
    PopSize = factor(PopSize, levels = c("1000", "10000"))
  )

(corMD <-
    ggplot(SumRunsF1000%>%filter(Generation<50),
           aes(
             x = Generation,
             y = cor_mean,
             color = Level,
             fill = Level
           )) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = cor_mean - cor_sd,
        ymax = cor_mean + cor_sd
      ),
      alpha = 0.3,
      color = NA
    ) +
    geom_hline(data=SumRunsF20%>%filter(Generation==50),aes(yintercept=cor_mean,color=Level),linetype="dotted")+
    facet_nested(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A = A.labs)) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal) +
    ylab("Genetic correlation between f and m") + 
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))


(Bm <-
    ggplot(SumRunsF1000%>%filter(Generation<50),
           aes(
             x = Generation,
             y = BMale_mean,
             color = Level,
             fill = Level
           )) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = BMale_mean - BMale_sd,
        ymax = BMale_mean + BMale_sd
      ),
      alpha = 0.3,
      color = NA
    ) +
    facet_nested(Tradeoff ~ A, labeller = labeller(Tradeoff = T.labs, A = A.labs)) +
    scale_color_manual(name = "Risk of sperm competition (q)", values = pal) +
    scale_fill_manual(name = "Risk of sperm competition (q)", values = pal) +
    ylab("Beta selection estimate of m") + 
    geom_hline(aes(yintercept = 0), color ="black", linetype = "dashed") + scale_x_continuous(n.breaks = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))



(((Bm/corMD)) &
    theme(panel.spacing.y = unit(0.5, "lines")) &
    guides(color = guide_legend(override.aes = list(
      size = 5, alpha = 1
    )), color = "none")
) + plot_annotation(tag_levels = "A", tag_suffix = ")") + plot_layout(guides = "collect")

#ggsave(
#  "FigureS7.pdf",
#  height = 183 * 2,
#  width = 183 * 1.8,
#  units = "mm"
#)

# 11. SI Figure 8: Sensitivity analyse of deviations -------------------------------------------------------------
Runs_final<-runs_rawS %>%
  filter(Generation == 30000)%>%merge(PredAn,by=c("Level","Tradeoff"))%>%
  mutate(investment=ifelse(Tradeoff=="true",MeanCount*MeanMale,MeanCount),diff=MeanCount-PredN,diff2=investment-PredT,relDif=diff2/PredT)%>%
  group_by(Level, A, Tradeoff, PopSize, Var, Loci) %>%
  summarise(
    relDif_median = median(relDif))%>%
  mutate(Var=factor(Var,levels=c("LV","MV","HV")))

#Make plot
ggplot(Runs_final, aes(x = Level, y = A, fill = relDif_median)) +
  geom_tile(color = "black", size = 0.5) +
  facet_nested(
    Tradeoff + PopSize ~ Var + Loci,
    labeller = labeller(
      Tradeoff = T.labs,
      Var = Var.labs,
      Loci = L.labs,
      PopSize = P.labs
    )
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "right", aspect.ratio = 1) +
  scale_y_discrete(labels = c("Weak", "Moderate", "Strong"),
                   expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#E66100",
    mid="white",
    high="#5D3A9B",
    #limits = ,c(-0.01, .2),
    guide = guide_colorbar(
      frame.colour = "black",
      frame.linewidth = 2,
      ticks.linewidth = 2,
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 18
    )
  ) +
  labs(y = "Strength of cryptic preference", x = "Risk of sperm competition (q)", fill =
         "Relative difference \n from predicted investment ") +
  theme(legend.position = "top")


#Save figure
#ggsave(
#  "SIFigure8.pdf",
#  height = 200,
#  width = 180 * 2,
#  units = "mm"
#)

# 12. SI Figure 9: Sensitivity of within-population cors (SI Figure 5) ---------------------------------------------------------
#summarize data to get medians across 50 popultions of median within-population correlations over 2,000 generations for each parameter combination

runs_rawSumW <- runs_rawSumLast %>%
  group_by(Level, A, Tradeoff, PopSize, Var, Loci) %>%
  summarise(
    cor_median = median(cor_median))
#plot the figure
ggplot(runs_rawSumW, aes(x = Level, y = A, fill = cor_median)) +
  geom_tile(color = "black", size = 0.5) +
  facet_nested(
    Tradeoff + PopSize ~ Var + Loci,
    labeller = labeller(
      Tradeoff = T.labs,
      Var = Var.labs,
      Loci = L.labs,
      PopSize = P.labs
    )
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "right", aspect.ratio = 1) +
  scale_y_discrete(labels = c("Weak", "Moderate", "Strong"),
                   expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_distiller(
    palette = "Purples",
    na.value = "grey",
    direction = 1,
    limits = c(-0.01, .2),
    guide = guide_colorbar(
      frame.colour = "black",
      frame.linewidth = 2,
      ticks.linewidth = 2,
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 18
    )
  ) +
  labs(y = "Strength of cryptic preference", x = "Risk of sperm competition (q)", fill =
         "Within-population correlation between f and m") +
  theme(legend.position = "top")


#Save figure
#ggsave(
#  "SIFigure9.pdf",
#  height = 200,
#  width = 180 * 2,
#  units = "mm"
#)


# 13. SI Figure 10. Sensitivity of gamma estimates -------------------------------------------------------------
runs_rawSumLast %>%
  group_by(Level, A, Tradeoff, PopSize, Var, Loci) %>%
  summarise(
    GMale_medians = median(GMale_median))%>%
  mutate(Var=factor(Var,levels=c("LV","MV","HV")))%>%
  ggplot(aes(x = Level, y = A, fill = GMale_medians)) +
  geom_tile(color = "black", size = 0.5) +
  facet_nested(
    Tradeoff + PopSize ~ Var + Loci,
    labeller = labeller(
      Tradeoff = T.labs,
      Var = Var.labs,
      Loci = L.labs,
      PopSize = P.labs
    )
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "right", aspect.ratio = 1) +
  scale_y_discrete(labels = c("Weak", "Moderate", "Strong"),
                   expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#E66100",
    mid="white",
    high="#5D3A9B",
    #limits = ,c(-0.01, .2),
    guide = guide_colorbar(
      frame.colour = "black",
      frame.linewidth = 2,
      ticks.linewidth = 2,
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 18
    )
  ) +
  labs(y = "Strength of cryptic preference", x = "Risk of sperm competition (q)", fill =
         "Gamma sperm trait (m)") +
  theme(legend.position = "top")


#Save figure
#ggsave(
#  "SIFigure10.png",
#  height = 200,
#  width = 180 * 2,
#  units = "mm"
#)


# 14. SI Figure 11. Sensitivity of lags --------------------------------------------------
LagAll %>%
  group_by(Level, A, Tradeoff, PopSize, Var, Loci) %>%
  summarise(
    MinLag_median = median(MinLag))%>%
  mutate(Var=factor(Var,levels=c("LV","MV","HV")))%>%
  ggplot(aes(x = Level, y = A, fill = MinLag_median)) +
  geom_tile(color = "black", size = 0.5) +
  facet_nested(
    Tradeoff + PopSize ~ Var + Loci,
    labeller = labeller(
      Tradeoff = T.labs,
      Var = Var.labs,
      Loci = L.labs,
      PopSize = P.labs
    )
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "right", aspect.ratio = 1) +
  scale_y_discrete(labels = c("Weak", "Moderate", "Strong"),
                   expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#E66100",
    mid="white",
    high="#5D3A9B",
    #limits = ,c(-0.01, .2),
    guide = guide_colorbar(
      frame.colour = "black",
      frame.linewidth = 2,
      ticks.linewidth = 2,
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 18
    )
  ) +
  labs(y = "Strength of cryptic preference", x = "Risk of sperm competition (q)", fill =
         "Generational lag with minimum MAE") +
  theme(legend.position = "top")

#ggsave(
#  "SIFigure11.png",
#  height = 200,
#  width = 180 * 2,
#  units = "mm"
#)

