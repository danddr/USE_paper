# Daniele Da Re, Manuele Bazzichetto, Enrico Tordoni
# 20230616
#---- 1. Load libraries and external functions -----
library(raster)
library(sf)
library(USE)
library(virtualspecies)
library(dismo)
library(sdm)
library(ranger)
library(tidyverse)
library(foreach)
library(parallel)
library(doParallel)
library(ecospat)
library(Metrics)

source("script/d_squared.R")
#---- 2. Download bioclimatic rasters -----
#download bioclimatic raster
# Worldclim<-raster::getData('worldclim', var='bio', res=10)
Worldclim <- geodata::worldclim_global(var='bio', res=10, path=getwd()) 
envData <- terra::crop(Worldclim, terra::ext(-12, 25, 36, 60))
rm(Worldclim)
rpc <- USE::rastPCA(envData, nPC = 2, naMask = TRUE, stand = TRUE)
dt <- na.omit(terra::as.data.frame(rpc$PCs, xy=TRUE))
dt <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))
myRes<-USE::optimRes(sdf=dt, grid.res=c(1:12), perc.thr = 10, showOpt = TRUE, cr=4)

#---- 3. Generate virtual species, create pseudo-absences using different sampling approaches, SDMs exercise and models statistics -----
nVirtspecies <- 10
percTesting <- 30
Nreps <- 5
myFold <- 5
myCores <-5
BuffSize <- c(50*10^3, 100*10^3, 200*10^3) #in meters
myCRS <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
Vsprev <- 0.25 # species geographic prevalence 
bgk_prev <- 1  # sample prevalence 

envData <- raster::stack(envData) # convert to raster object, otherwise virtualspecies R pkg doesn't work

#run the whole framework
myVirtualSP_list <- list()

for(myVs in 1:nVirtspecies){
  # myVs=1
  
  # create widespread virtual species
  myRandNum <- sample(1:19,size=5, replace = FALSE)
  random.sp <- virtualspecies::generateRandomSp(envData[[myRandNum]], 
                                                convert.to.PA = FALSE, 
                                                species.type = "additive",
                                                realistic.sp = TRUE, 
                                                plot = FALSE)
  
  #reclassify suitability raster using a probability conversion rule
  new.pres <- convertToPA(x=random.sp,
                          beta="random",
                          alpha = -0.05, plot = FALSE , 
                          species.prevalence = Vsprev)
  
  # Sample true occurrences
  presence.points <- sampleOccurrences(new.pres,
                                       n = 300, # The number of points to sample
                                       type = "presence only",
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  
  plot(new.pres$pa.raster, main=paste("species #",   myVs))
  
  myPres <- presence.points$sample.points[, c( "x", "y",  "Observed")]
  coordinates(myPres)<- ~x+y
  raster::crs(myPres)<- myCRS
  
  # loop on two sets of environmental dataset: one with 5 bioclimatic variables and one with 19 bioclimatic variables
  env_datasets <- list(envData, envData[[random.sp$details$variables]])
  myOut_list <- list()
  
  for(k in 1:length(env_datasets)){
    # k=2
    # Create pseudo absences using kernel and environmental distances
    set.seed(666)
    
    #grid.psAbs_sampling in grid
    myGrid.psAbs <- USE::paSampling(env_datasets[[k]], 
                                    pres=myPres, 
                                    thres=0.75,  
                                    H=NULL,
                                    grid.res=myRes$Opt_res, 
                                    prev=bgk_prev, sub.ts=FALSE, 
                                    plot_proc=FALSE)
    
    #only keep coords
    myGrid.psAbs <- data.frame(x = myGrid.psAbs$x, y = myGrid.psAbs$y)
    coordinates(myGrid.psAbs) <- ~x+y
    raster::crs(myGrid.psAbs) <- myCRS
    # plot(myGrid.psAbs, add=T)
    obs_prev <- round(nrow(myPres)/length(myGrid.psAbs),2)
    
    # create pseudoabs from buffer
    buf<-lapply(BuffSize, function(x){circles(myPres, d = x, lonlat = TRUE)})

    # buffer OUT, constrained by convex hull
    chull <- myPres %>% st_as_sf() %>% 
      summarise( geometry = st_combine( geometry ) ) %>%
      st_convex_hull()
    chull <- as(st_geometry(chull), "Spatial")
    
    buf.mask <-lapply(buf, function(x){mask(envData[[1]], rgeos::gUnaryUnion(x@polygons),inverse=TRUE)})
    buf.mask <- lapply(buf.mask, function(x){mask(x, chull)})
    myPseudo_buf <- lapply(buf.mask, function(x){as.data.frame(sampleRandom(x, 
                                                                            size = round(length(myPres)/bgk_prev) , 
                                                                            xy = TRUE))})
    #create the training datasets
    b <- as(extent(envData), "SpatialPolygons")
    myPseudo_rand <- spsample(b, n =  round(length(myPres)/bgk_prev), type = "random")
    myPseudo_rand <- SpatialPointsDataFrame(myPseudo_rand,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
      
    myPseudo_buf.50km <- myPseudo_buf[[1]]
    coordinates(myPseudo_buf.50km)<-~x+y
    myPseudo_buf.50km <- SpatialPointsDataFrame(myPseudo_buf.50km,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    
    myPseudo_buf.100km <- myPseudo_buf[[2]]
    coordinates(myPseudo_buf.100km)<-~x+y
    myPseudo_buf.100km <- SpatialPointsDataFrame(myPseudo_buf.100km,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    
    myPseudo_buf.200km <- myPseudo_buf[[3]]
    coordinates(myPseudo_buf.200km)<-~x+y
    myPseudo_buf.200km <- SpatialPointsDataFrame(myPseudo_buf.200km,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    
    myGrid.psAbs$Observed <- rep(0, nrow(myGrid.psAbs@coords))
    
    proj4string(myPseudo_rand) <- raster::crs(myPres)
    proj4string(myPseudo_buf.50km) <- raster::crs(myPres)
    proj4string(myPseudo_buf.100km) <- raster::crs(myPres)
    proj4string(myPseudo_buf.200km) <- raster::crs(myPres)
    proj4string(myGrid.psAbs) <- raster::crs(myPres)
    
    myPAlist <- list(
      "mytrain_rand" = myPres+myPseudo_rand, 
      "mytrain_buf_out50" = myPres+myPseudo_buf.50km,
      "mytrain_buf_out100" = myPres+myPseudo_buf.100km,
      "mytrain_buf_out200" = myPres+myPseudo_buf.200km,
      "mytrain_grid" = myPres+myGrid.psAbs
      )
    
    #dataset for predictions
    mypredstack <- env_datasets[[k]]
    
    #count proportion of correct absences
    prop_corr_abs <- lapply(myPAlist, function(x){y=subset(x, Observed==0); z=raster::extract(new.pres$pa.raster, y); 
    return(length(which(z==0))/length(z))})
    
    #modelling
    myOut_tmp <- data.frame() 
    for(i in 1:length(myPAlist)){
      # i=1
      cat(paste("\nModel ",i,"/",length(myPAlist),"\n"))
      
      #extract enviromental data
      d <- sdmData(formula=Observed~., train=myPAlist[[i]], predictors=mypredstack)
      
      #model          
      m1 <- sdm(Observed~.,data=d,
                methods=c('glm', 'rf' ,'gam' , 'brt', 'Maxent'),
                replication='sub',test.percent=percTesting,n=Nreps, parallelSettings=list("parallel", ncore=myCores))
      
      #store info of the models statistics (internal validation on the subsets)
      p <- length(which(myPAlist[[i]]$Observed==1))
      pseabs <- length(which(myPAlist[[i]]$Observed==0))
      tmp_out <- setNames(data.frame(matrix(ncol = 21, nrow = 0)), c("Species", "trainig_set", "exp_prevalence",  "obs_prevalence", "pres",'abs', "prop_corr_abs", "nDrivingFactors", 
                                                                     "ModelType", "nPredictorsPCA", "nPredictorsModel", 
                                                                     "AUC", "COR", "Deviance", "TSS", "Kappa", 
                                                                     "Sensitivity", "Specificity", "BoyceI", 'RMSE',"R2"   ))
      mynames <- names(tmp_out)
      
      if(is.null(unlist(lapply(m1@models$Observed$glm, function(x){x@errorLog})))){
        glm_rows=rbind(c(myVs, names(myPAlist[i]),bgk_prev,obs_prev,  p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "GLM", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), 
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@statistics$COR[1]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@statistics$Deviance) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$Kappa[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){ mean(ecospat.boyce(fit =  x@evaluation$test.dep@predicted, 
                                                                                                    obs=   x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)], 
                                                                                                    nclass=0,  
                                                                                                    window.w="default", res=100, PEplot = FALSE,  method = 'pearson')$cor) 
                         })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){
                           Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$glm, function(x){Dsquared(x@object,adjust = TRUE)})), na.rm=T)))
      } else {glm_rows=c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "GLM", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,9))}
      
      
      if(is.null(unlist(lapply(m1@models$Observed$gam, function(x){x@errorLog})))){
        gam_rows=rbind(c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "GAM", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@statistics$COR[1]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@statistics$Deviance) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@threshold_based$Kappa[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){ mean(ecospat.boyce(fit =  x@evaluation$test.dep@predicted, 
                                                                                                    obs=   x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)], 
                                                                                                    nclass=0,  
                                                                                                    window.w="default", res=100, PEplot = FALSE,  method = 'pearson')$cor) 
                         })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){
                           Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$gam, function(x){summary(x@object)$r.sq})), na.rm=T)))
      } else {gam_rows=c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "GAM", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,10))}
      
      if(is.null(unlist(lapply(m1@models$Observed$rf, function(x){x@errorLog})))){ 
        rf_rows=rbind(c(myVs, names(myPAlist[i]),bgk_prev,obs_prev, p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "RF", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), 
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@statistics$COR[1]) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@statistics$Deviance) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@threshold_based$Kappa[2]) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(ecospat.boyce(fit =  x@evaluation$test.dep@predicted, 
                                                                                                  obs=   x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)], 
                                                                                                  nclass=0,  
                                                                                                  window.w="default", res=100, PEplot = FALSE, method = 'pearson')$cor) 
                        })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ sqrt(mean(x@object$mse)) })), na.rm=T),
                        mean(unlist(lapply(m1@models$Observed$rf, function(x){ mean(x@object$rsq) })), na.rm=T)))
      } else {rf_rows=c(myVs, names(myPAlist[i]),bgk_prev,obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "RF", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,10))}
      
      if(is.null(unlist(lapply(m1@models$Observed$brt, function(x){x@errorLog})))){ 
        brt_rows=rbind(c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "BRT", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@statistics$COR[1]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@statistics$Deviance) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$Kappa[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){ mean(ecospat.boyce(fit =  x@evaluation$test.dep@predicted, 
                                                                                                    obs=   x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)], 
                                                                                                    nclass=0,  
                                                                                                    window.w="default", res=100, PEplot = FALSE,  method = 'pearson')$cor) 
                         })), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){
                           Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                         mean(unlist(lapply(m1@models$Observed$brt, function(x){
                           round(cor(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)^2,3)})), na.rm=T)))
      } else {brt_rows=c(myVs, names(myPAlist[i]),bgk_prev,obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "BRT", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,10))}
      
      if(is.null(unlist(lapply(m1@models$Observed$maxent, function(x){x@errorLog})))){ 
        maxent_rows=rbind(c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "Maxent", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@statistics$COR[1]) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@statistics$Deviance) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@threshold_based$Kappa[2]) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){ mean(ecospat.boyce(fit =  x@evaluation$test.dep@predicted, 
                                                                                                          obs=   x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)], 
                                                                                                          nclass=0,  
                                                                                                          window.w="default", res=100, PEplot = FALSE, method = 'pearson')$cor) 
                            })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){
                              Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){
                              round(cor(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)^2,3)})), na.rm=T)))
      } else {maxent_rows=c(myVs, names(myPAlist[i]),bgk_prev,obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "Maxent", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,10))}
      
      
      
      tmp_out <- rbind(tmp_out, glm_rows, gam_rows, rf_rows, brt_rows, maxent_rows)
      names(tmp_out) <- mynames
      
      #compute spatial predictions and rmse
      message("\n Compute spatial predictions and accuracy metrics on the observed VS distribution \n")
      m1.preds <- predict(m1, newdata=mypredstack, overwrite=TRUE)
      m1.preds <- na.omit(raster::as.data.frame(raster::stack(new.pres$suitab.raster, m1.preds)))
      rmse.distribution <- apply(m1.preds[,2:ncol(m1.preds)], 2, function(x){round(Metrics::rmse(m1.preds$layer,x),3)})
      cor.distribution <- apply(m1.preds[,2:ncol(m1.preds)], 2, function(x){round(cor(m1.preds$layer,x),3)})
      rmse.distribution <- data.frame(ModelType=attributes(rmse.distribution)$names, 
                                      rmse.distribution=rmse.distribution, 
                                      cor.distribution=cor.distribution)
      rmse.distribution$ModelType <- gsub(sapply(rmse.distribution$ModelType, function(x){strsplit(x, "_")[[1]][[4]]}), 
                                          pattern=".re", replacement = "")
      
      rmse.distribution <- rmse.distribution %>% 
        group_by(ModelType) %>% 
        summarise(rmse.distribution=mean(rmse.distribution, na.rm=TRUE), 
                  cor.distribution=mean(cor.distribution, na.rm=TRUE)) %>% 
        mutate(ModelType=toupper(ModelType), 
               ModelType=ifelse(ModelType=="MAXENT", "Maxent", ModelType))
      
      tmp_out<-tmp_out %>% 
        left_join(rmse.distribution, by="ModelType")
      
      myOut_tmp <- rbind(myOut_tmp, tmp_out)
      
    }
    #store myOut_tmp for different environmental datasets used for modelling
    myOut_list[[k]] <- myOut_tmp
  }
  #combine list and store it by species
  myOut_list <- do.call(rbind, myOut_list)
  myVirtualSP_list[[myVs]] <- myOut_list
}

length(myVirtualSP_list)
myVirtualSP_list[[1]]

outname <- paste0("simulations/paper_sim/June2023/",nVirtspecies, "vs_prev", sub("\\.", "", as.character(bgk_prev)), "_occOnly_", "radiusSensitivity_", Sys.Date(), ".RDS")
saveRDS(myVirtualSP_list, outname)

#---- 1. load 50 virtual species simulations results----
myVirtualSP_list <- readRDS("10vs_prev1_occOnly_radiusSensitivity_2023-06-15.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

#pivot longer and boxplot
tmp <-   myProc_out %>%
  as_tibble() %>%
  dplyr::select(trainig_set,  
                AUC, TSS, nPredictorsModel, ModelType, rmse.distribution, Sensitivity, Specificity, BoyceI) %>% 
  mutate(trainig_set=as.factor(trainig_set),
         nPredictorsModel=as.factor(nPredictorsModel),
         BoyceI=round(as.numeric(gsub(",", "", BoyceI)),3),
         AUC=round(as.numeric(gsub(",", "", AUC)),3),
         TSS=round(as.numeric(gsub(",", "", TSS)),3),
         Sensitivity=round(as.numeric(gsub(",", "", Sensitivity)),3) , 
         Specificity=round(as.numeric(gsub(",", "", Specificity)),3)) %>% 
  group_by(trainig_set, ModelType, nPredictorsModel) %>% 
  pivot_longer(cols = c(AUC, BoyceI, Sensitivity, Specificity, rmse.distribution, TSS )) %>% 
  mutate(trainig_set=factor(trainig_set, levels= c( "mytrain_grid", "mytrain_rand", "mytrain_buf_out50", "mytrain_buf_out100", "mytrain_buf_out200")), 
         ModelType=factor(ModelType, levels= c("GLM", "GAM", "RF", "BRT", "Maxent")), 
         name=factor(name, levels= c("AUC", "BoyceI", "Sensitivity", "Specificity", "TSS","rmse.distribution")),
         new_name_tr_set=dplyr::recode(trainig_set, mytrain_grid= "Uniform", mytrain_rand= "Random",  mytrain_buf_out50="Buffer-out50km", mytrain_buf_out100="Buffer-out100km", mytrain_buf_out200="Buffer-out200km")) %>%
  mutate(name=dplyr::recode(name, BoyceI = "CBI"), 
         name=dplyr::recode(name,  rmse.distribution = "RMSE"))

#---- 2. Influence of the pseudo-absences sampling approach on HSMs predictive accuracy statistics ----
Npred <- 5
mySpecies <- 10
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
  scale_color_manual(breaks = c("Uniform", "Random","Buffer-out50km", "Buffer-out100km", "Buffer-out200km" ),
                     values=c( "#0072B2", "#E69F00",  "#CC79A7", "#824E6B", "#422837"))+
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
outname <- paste0(outdir, mySpecies, "sp_prev", sub("\\.", "", as.character(myPrev)), "_", Npred,"pred_RadiusSens_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 20, height = 15, device='png', dpi=320)
