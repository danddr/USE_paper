#Violin plot figures
library(tidyverse)
library(ggplot2)

#---- 1. load 50 virtual species simulations results----
myVirtualSP_list=readRDS("50vs_prev1_occOnly_radius50km_2022-09-14.RDS")
myProc_out=do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp=myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set, exp_prevalence, 
                AUC, COR, Deviance,  
                TSS, nPredictorsModel, ModelType,RMSE,R2, Kappa, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         exp_prevalence=as.factor(exp_prevalence), 
         nPredictorsModel=as.factor(nPredictorsModel)) %>% 
  mutate(BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3) ,
         AUC=round(as.numeric(gsub(",", "", AUC)),3) , 
         COR=round(as.numeric(gsub(",", "", COR)),3) , 
         Deviance=round(as.numeric(gsub(",", "", Deviance)),3) , 
         TSS=round(as.numeric(gsub(",", "", TSS)),3) ,
         RMSE=round(as.numeric(gsub(",", "", RMSE)),3) ,
         R2=round(as.numeric(gsub(",", "", R2)),3) , 
         Kappa=round(as.numeric(gsub(",", "", Kappa)),3) , 
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  mutate(R2=ifelse(R2<0, 0, R2)) %>% 
  group_by(training_set, exp_prevalence,ModelType, nPredictorsModel) %>%
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, Kappa, TSS, Deviance, COR, R2,RMSE )) %>% 
  filter(name !="Deviance") %>% 
  mutate(training_set=factor(training_set, levels= c("myTrue_Abs", "mytrain_grid", "mytrain_rand", "mytrain_buf_in", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "Kappa", "TSS", "R2", "COR", "RMSE")))%>% 
  mutate(new_name_tr_set=dplyr::recode(training_set, myTrue_Abs = "True Absences", mytrain_grid= "Uniform", mytrain_rand= "Random", mytrain_buf_in= "BufferIn", mytrain_buf_out="BufferOut")) %>% 
  mutate(name=dplyr::recode(name, BoyceI = "CBI"))

#---- 2. Influence of the pseudo-absences sampling approach on HSMs predictive accuracy statistics ----
Npred=5
mySpecies=50
myPrev=1
mytitle=paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
unique(tmp$new_name_tr_set)[1]

p=tmp %>% 
  filter(training_set!="myTrue_Abs") %>% 
  filter(nPredictorsModel == Npred) %>% 
  filter(name!="R2", name!="RMSE",  name!="COR", name!="Kappa") %>% 
  ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+
  geom_violin()+
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray")+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "BufferOut" ),
                     values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
  facet_grid(ModelType~name, scales = "free_y")+
  ylim(-0.5,1)+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  # ggtitle(mytitle)+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(2, 'cm'))

p
outname=paste0(mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)

#---- 3. Sensitivity analysis on the sample prevalence ----
#prevalence 0.5
myVirtualSP_list=readRDS("10vs_prev05_occOnly_radius50km_2022-09-14.RDS")
myProc_out=do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp=myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set, exp_prevalence, 
                AUC, COR, Deviance,  
                TSS, nPredictorsModel, ModelType,RMSE,R2, Kappa, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         exp_prevalence=as.factor(exp_prevalence), 
         nPredictorsModel=as.factor(nPredictorsModel)) %>% 
  mutate(BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3) ,
         AUC=round(as.numeric(gsub(",", "", AUC)),3) , 
         COR=round(as.numeric(gsub(",", "", COR)),3) , 
         Deviance=round(as.numeric(gsub(",", "", Deviance)),3) , 
         TSS=round(as.numeric(gsub(",", "", TSS)),3) ,
         RMSE=round(as.numeric(gsub(",", "", RMSE)),3) ,
         R2=round(as.numeric(gsub(",", "", R2)),3) , 
         Kappa=round(as.numeric(gsub(",", "", Kappa)),3) , 
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  mutate(R2=ifelse(R2<0, 0, R2)) %>% 
  group_by(training_set, exp_prevalence,ModelType, nPredictorsModel) %>%
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, Kappa, TSS, Deviance, COR, R2,RMSE )) %>% 
  filter(name !="Deviance") %>% 
  mutate(training_set=factor(training_set, levels= c("myTrue_Abs", "mytrain_grid", "mytrain_rand", "mytrain_buf_in", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "Kappa", "TSS", "R2", "COR", "RMSE")))%>% 
  mutate(new_name_tr_set=dplyr::recode(training_set, myTrue_Abs = "True Absences", mytrain_grid= "Uniform", mytrain_rand= "Random", mytrain_buf_in= "BufferIn", mytrain_buf_out="BufferOut")) %>% 
  mutate(name=dplyr::recode(name, BoyceI = "CBI"))

Npred=5
mySpecies=10
myPrev=0.5
mytitle=paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p=tmp %>% 
  filter(training_set!="myTrue_Abs") %>% 
  filter(nPredictorsModel == Npred) %>% 
  filter(name!="R2", name!="RMSE",  name!="COR", name!="Kappa") %>% 
  ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+
  geom_violin()+
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray")+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "BufferOut" ),
                     values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
  facet_grid(ModelType~name, scales = "free_y")+
  ylim(-0.5,1)+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  # ggtitle(mytitle)+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(2, 'cm'))

p

outname=paste0(mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)

# prevalence 0.1 
myVirtualSP_list=readRDS("10vs_prev01_occOnly_radius50km_2022-09-14.RDS")
myProc_out=do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp=myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set, exp_prevalence, 
                AUC, COR, Deviance,  
                TSS, nPredictorsModel, ModelType,RMSE,R2, Kappa, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         exp_prevalence=as.factor(exp_prevalence), 
         nPredictorsModel=as.factor(nPredictorsModel)) %>% 
  mutate(BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3) ,
         AUC=round(as.numeric(gsub(",", "", AUC)),3) , 
         COR=round(as.numeric(gsub(",", "", COR)),3) , 
         Deviance=round(as.numeric(gsub(",", "", Deviance)),3) , 
         TSS=round(as.numeric(gsub(",", "", TSS)),3) ,
         RMSE=round(as.numeric(gsub(",", "", RMSE)),3) ,
         R2=round(as.numeric(gsub(",", "", R2)),3) , 
         Kappa=round(as.numeric(gsub(",", "", Kappa)),3) , 
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  mutate(R2=ifelse(R2<0, 0, R2)) %>% 
  group_by(training_set, exp_prevalence,ModelType, nPredictorsModel) %>%
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, Kappa, TSS, Deviance, COR, R2,RMSE )) %>% 
  filter(name !="Deviance") %>% 
  mutate(training_set=factor(training_set, levels= c("myTrue_Abs", "mytrain_grid", "mytrain_rand", "mytrain_buf_in", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "Kappa", "TSS", "R2", "COR", "RMSE")))%>% 
  mutate(new_name_tr_set=dplyr::recode(training_set, myTrue_Abs = "True Absences", mytrain_grid= "Uniform", mytrain_rand= "Random", mytrain_buf_in= "BufferIn", mytrain_buf_out="BufferOut")) %>% 
  mutate(name=dplyr::recode(name, BoyceI = "CBI"))

Npred=5
mySpecies=10
myPrev=0.1
mytitle=paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p=tmp %>% 
  filter(training_set!="myTrue_Abs") %>% 
  filter(nPredictorsModel == Npred) %>% 
  filter(name!="R2", name!="RMSE",  name!="COR", name!="Kappa") %>% 
  ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+
  geom_violin()+
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray")+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "BufferOut" ),
                     values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
  facet_grid(ModelType~name, scales = "free_y")+
  ylim(-0.5,1)+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  # ggtitle(mytitle)+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(2, 'cm'))

p

outname=paste0(mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)


