#Violin plot figures
library(tidyverse)
library(ggplot2)

#---- 1. load 50 virtual species simulations results----
myVirtualSP_list <- readRDS("50vs_prev1_occOnly_radius50km_2023-06-15.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp <- myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set,  
                AUC, TSS, nPredictorsModel, ModelType, rmse.distribution, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         nPredictorsModel=as.factor(nPredictorsModel),
         BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3),
         AUC=round(as.numeric(gsub(",", "", AUC)),3),
         TSS=round(as.numeric(gsub(",", "", TSS)),3),
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  group_by(training_set, ModelType, nPredictorsModel) %>% 
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, rmse.distribution, TSS )) %>% 
  mutate(training_set=factor(training_set, levels= c( "mytrain_grid", "mytrain_rand", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "TSS","rmse.distribution")),
         new_name_tr_set=dplyr::recode(training_set, mytrain_grid= "Uniform", mytrain_rand= "Random",  mytrain_buf_out="Buffer-out")) %>%
  mutate(name=dplyr::recode(name, BoyceI = "CBI"), 
         name=dplyr::recode(name,  rmse.distribution = "RMSE"))

#---- 2. Influence of the pseudo-absences sampling approach on HSMs predictive accuracy statistics ----
Npred <- 5
mySpecies <- 50
myPrev <- 1
mytitle <- paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
unique(tmp$new_name_tr_set)[1]

p <- tmp %>% 
      filter(nPredictorsModel == Npred) %>% 
      ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+ 
      geom_violin()+
      stat_summary(fun = median, geom = "point", size = 2) +
      scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "Buffer-out" ),
                         values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
      labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
      facet_grid(ModelType~name, scales = "free_y")+
      ylim(0,1)+
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
            legend.key.size = unit(1.5, 'cm'))

p
outname <- paste0(outdir, mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)

#---- 3. Sensitivity analysis on the sample prevalence ----
#prevalence 0.5
myVirtualSP_list <- readRDS("10vs_prev05_occOnly_radius50km_2023-06-14.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp <- myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set,  
                AUC, TSS, nPredictorsModel, ModelType, rmse.distribution, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         nPredictorsModel=as.factor(nPredictorsModel),
         BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3),
         AUC=round(as.numeric(gsub(",", "", AUC)),3),
         TSS=round(as.numeric(gsub(",", "", TSS)),3),
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  group_by(training_set, ModelType, nPredictorsModel) %>% 
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, rmse.distribution, TSS )) %>% 
  mutate(training_set=factor(training_set, levels= c( "mytrain_grid", "mytrain_rand", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "TSS","rmse.distribution")),
         new_name_tr_set=dplyr::recode(training_set, mytrain_grid= "Uniform", mytrain_rand= "Random",  mytrain_buf_out="Buffer-out")) %>%
  mutate(name=dplyr::recode(name, BoyceI = "CBI"), 
         name=dplyr::recode(name,  rmse.distribution = "RMSE"))

Npred <- 5
mySpecies <- 10
myPrev <- 0.5
mytitle <- paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- tmp %>% 
  filter(nPredictorsModel == Npred) %>% 
  ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+ 
  geom_violin()+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c( "Uniform", "Random",   "Buffer-out" ),
                     values=c( "#0072B2", "#E69F00",  "#CC79A7" ))+
  labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
  facet_grid(ModelType~name, scales = "free_y")+
  ylim(0,1)+
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
        legend.key.size = unit(1.5, 'cm'))

p

outname <- paste0(outdir, mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)

# prevalence 0.1 
myVirtualSP_list <- readRDS("10vs_prev01_occOnly_radius50km_2023-06-14.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp <- myProc_out %>% 
  as_tibble() %>%
  dplyr::select(training_set,  
                AUC, TSS, nPredictorsModel, ModelType, rmse.distribution, Sensitivity, Specificity, BoyceI) %>% 
  mutate(training_set=as.factor(training_set),
         nPredictorsModel=as.factor(nPredictorsModel),
         BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3),
         AUC=round(as.numeric(gsub(",", "", AUC)),3),
         TSS=round(as.numeric(gsub(",", "", TSS)),3),
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  group_by(training_set, ModelType, nPredictorsModel) %>% 
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, rmse.distribution, TSS )) %>% 
  mutate(training_set=factor(training_set, levels= c( "mytrain_grid", "mytrain_rand", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "TSS","rmse.distribution")),
         new_name_tr_set=dplyr::recode(training_set, mytrain_grid= "Uniform", mytrain_rand= "Random",  mytrain_buf_out="Buffer-out")) %>%
  mutate(name=dplyr::recode(name, BoyceI = "CBI"), 
         name=dplyr::recode(name,  rmse.distribution = "RMSE"))

Npred <- 5
mySpecies <- 10
myPrev <- 0.1
mytitle <- paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- tmp %>% 
  filter(nPredictorsModel == Npred) %>% 
  ggplot(aes(new_name_tr_set , value, color=new_name_tr_set))+ 
  geom_violin()+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "Buffer-out" ),
                     values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="", y="Predictive accuracy metrics value", color="Sampling method")+
  facet_grid(ModelType~name, scales = "free_y")+
  ylim(0,1)+
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
        legend.key.size = unit(1.5, 'cm'))
p

outname <- paste0(outdir, mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)


