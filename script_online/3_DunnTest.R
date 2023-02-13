# Dunn test for multiple comparison
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(multcomp)
library(data.table)
library(dunn.test)

#---- pairwise comparisons for model performance -----
myVirtualSP_list <-readRDS('50vs_prev1_occOnly_radius50km_2022-09-14.RDS')
myProc_out=do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp5 =  myProc_out %>% 
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
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, Kappa, TSS, Deviance, COR, R2,RMSE)) %>% 
  filter(name !="Deviance") %>% 
  mutate(trainig_set=factor(training_set, levels= c("myTrue_Abs", "mytrain_grid", "mytrain_rand", "mytrain_buf_in", "mytrain_buf_out")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "Kappa", "TSS", "R2", "COR", "RMSE"))) %>% 
  mutate(new_name_tr_set=recode(training_set, myTrue_Abs = "True Absences",
                                mytrain_grid= "Uniform", mytrain_rand= "Random", 
                                mytrain_buf_in= "Buffer IN", mytrain_buf_out="Buffer OUT")) %>% 
  mutate(new_name_tr_set=factor(new_name_tr_set,levels = c("Uniform", 'Random',"Buffer IN","Buffer OUT"))) %>% 
  filter(nPredictorsModel == 5) %>%  #filter 5 predictors
  filter(name!= "Kappa",  name!= "COR" , name!= "R2" , name!= "RMSE") %>% 
  filter(new_name_tr_set !="True Absences") %>% 
  mutate(name=factor(name))

tmp5

tmp5 %>% 
  as_tibble() %>% 
  mutate(Group=paste0(ModelType,'_',name)) %>% 
  dplyr::select(ModelType, name,new_name_tr_set,value,Group) %>% 
  split(.$Group) %>%
  map(~dunn.test(x = .$value, g =.$new_name_tr_set , method = "holm",list=TRUE,wrap=TRUE))->OutList 

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
  group_by(comparisons,Sign) %>% 
  tally() %>% 
  mutate(newname=recode(comparisons, `Buffer IN - Uniform`  = "BufferIn",
                        `Buffer OUT - Uniform`  = "BufferOut",
                        `Random - Uniform`  = "Random")) %>% 
  ggplot(aes(x=newname, y=n, fill=Sign))+
  geom_bar(stat = "identity", position="fill", width=0.5)+
  scale_fill_manual(values = c("#98c1d9", "#ffce72"))+
  labs(x="", y="Relative proportion", fill="Significance")+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(2, 'cm'))

p1

p2<-pwc.out %>% 
  group_by(comparisons,Sign, Metric, Model) %>% 
  tally() %>% 
  mutate(newname=recode(comparisons, `Buffer IN - Uniform`  = "BufferIn",
                        `Buffer OUT - Uniform`  = "BufferOut",
                        `Random - Uniform`  = "Random")) %>% 
  ggplot(aes(x=newname, y=n, fill=Sign))+
  geom_bar(stat = "identity", position="fill", width=0.5)+
  scale_fill_manual(values = c("#98c1d9", "#ffce72"))+
  labs(x="", y="Relative proportion", fill="Significance")+
  facet_grid(Model~Metric)+
  theme_light()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(2, 'cm')) ; p2


library(patchwork)
myLayout<-c("1##
            222")
pp<-(p1 / p2)+
  plot_layout(design = myLayout) +
  plot_annotation(tag_levels = 'A')

pp

outname=paste("dunnTest_pwc_", Sys.Date(),".png", sep="")
ggsave(pp, filename = outname, width = 16, height = 8, device='png', dpi=320)

outname=paste("dunnTest_pwc_", Sys.Date(),".csv", sep="")
pwc.out %>% 
  mutate(chi2=round(chi2,2), 
         P.adjusted=round(P.adjusted, 4), 
         P.adjusted=ifelse(P.adjusted< 0.001, "p<0.001",P.adjusted )) %>% 
  write_csv(outname)