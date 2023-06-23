# Dunn test for multiple comparison
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(multcomp)
library(data.table)
library(dunn.test)
#---- pairwise comparisons for model performance -----
myVirtualSP_list <- readRDS("50vs_prev1_occOnly_radius50km_2023-06-15.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp5 <- myProc_out %>% 
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
         new_name_tr_set=dplyr::recode(training_set, mytrain_grid= "Uniform", mytrain_rand= "Random",  mytrain_buf_out="BufferOut")) %>%
  mutate(name=dplyr::recode(name, BoyceI = "CBI"), 
         name=dplyr::recode(name,  rmse.distribution = "RMSE")) %>% 
  filter(nPredictorsModel == 5) %>%  #filter 5 predictors
  mutate(name=factor(name))

tmp5

OutList<-tmp5 %>% 
  as_tibble() %>% 
  mutate(Group=paste0(ModelType,'_',name)) %>% 
  dplyr::select(ModelType, name,new_name_tr_set,value,Group) %>% 
  split(.$Group) %>% 
  map(~dunn.test(x = .$value, g =.$new_name_tr_set , method = "holm",list=TRUE,wrap=TRUE)) 

pwc <- do.call(rbind,lapply(OutList,data.frame))
pwc$Model <- sapply(rownames(pwc),function(x) strsplit(x,'_')[[1]] [1])
pwc$Metric <- sapply(rownames(pwc),function(x) strsplit(x,'\\.|_')[[1]] [2])
head(pwc)

pwc.out<-pwc %>% 
  filter(str_detect(comparisons, 'Uniform')) %>% 
  dplyr::select(Model, Metric, comparisons, chi2, P.adjusted ) %>% 
  mutate(Sign=ifelse(P.adjusted< 0.05, "p<0.05", "ns"))

pwc.out %>% 
  group_by(comparisons,Sign) %>% 
  tally()

p1<-pwc.out %>% 
  group_by(comparisons, Sign) %>% 
  tally() %>% 
  mutate(newname=dplyr::recode(comparisons, #`Buffer IN - Uniform`  = "BufferIn",
                        `BufferOut - Uniform`  = "Buffer-out",
                        `Random - Uniform`  = "Random"), 
         newname=factor(newname, levels = c("Random", "Buffer-out"))) %>%
  ggplot(aes(x=newname, y=n, fill=Sign))+
  geom_bar(stat = "identity", position="fill", width=0.45)+
  scale_fill_manual(values = c("#98c1d9", "#ffce72"))+
  labs(x="", y="Relative proportion", fill="Significance")+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(1.5, 'cm'))

p1

p2<-pwc.out %>% 
  group_by(comparisons,Sign, Metric) %>% #Model 
  tally() %>% 
  mutate(newname=dplyr::recode(comparisons, 
                        `BufferOut - Uniform`  = "Buffer-out",
                        `Random - Uniform`  = "Random"),
         newname=factor(newname, levels = c("Random", "Buffer-out")),
         Metric=factor(Metric, levels=c("AUC", "CBI", "Sensitivity", "Specificity", "TSS", "RMSE" ))) %>% 
  ggplot(aes(x=newname, y=n, fill=Sign))+
  geom_bar(stat = "identity", position="fill", width=0.5)+
  scale_fill_manual(values = c("#98c1d9", "#ffce72"))+
  labs(x="", y="Relative proportion", fill="Significance")+
  facet_grid(~Metric)+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) ; p2


library(patchwork)
myLayout<-c("1##
            222")
pp<-(p1 / p2)+
  plot_layout(design = myLayout) +
  plot_annotation(tag_levels = 'A')

pp

outname <- paste(outdir, "dunnTest_pwc_", Sys.Date(),".png", sep="")
ggsave(pp, filename = outname, width = 14, height = 6, device='png', dpi=320)

outname <- paste(outdir, "dunnTest_pwc_", Sys.Date(),".csv", sep="")
pwc.out %>% 
  mutate(chi2=round(chi2,2), 
         P.adjusted=round(P.adjusted, 4), 
         P.adjusted=ifelse(P.adjusted< 0.001, "p<0.001",P.adjusted )) %>% 
  write_csv(outname)
