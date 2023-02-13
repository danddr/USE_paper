# Daniele Da Re, Manuele Bazzichetto, Enrico Tordoni
# 20230210

#---- 1. Load libraries and external functions -----
library(USE)
library(iSDM)
library(raster)
library(virtualspecies)
library(sdm)
library(ranger)
library(dismo)
library(rgeos)
library(sf)
library(Metrics)
library(ggplot2)
library(tidyverse)
library("foreach")
library(parallel)
library(doParallel)
library(RStoolbox)
library(data.table)
library(ecospat)
source("d_squared.R")

#---- 2. Download bioclimatic rasters -----
#download bioclimatic raster
Worldclim<-raster::getData('worldclim', var='bio', res=10) 
envData<-crop(Worldclim, extent(-12, 25, 36, 60))
envData
rm(Worldclim)
rpc<-rasterPCA(envData,spca = TRUE)
dt <- na.omit(data.table(as.data.frame(rpc$map[[c("PC1", "PC2")]], xy = TRUE)))
dt=st_as_sf(dt, coords = c("PC1", "PC2"))
myRes=USE::optimRes(sdf=dt, grid.res=c(1:12), perc.thr = 10, showOpt = TRUE, cr=4)

#---- 3. Generate virtual species, create pseudo-absences using different sampling approaches, SDMs exercise and models statistics -----
nVirtspecies=50
percTesting=30
Nreps=5
myFold=5
myCores=5
BuffSize<-50000 #in meters
myCRS<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
Vsprev=0.25 #species geographic prevalence 
bgk_prev= 1 #0.5 #sample prevalence 

#run the whole framework
myVirtualSP_list=list()

for(myVs in 1:nVirtspecies){
  # myVs=1
  
  # create widespread virtual species
  myRandNum=sample(1:19,size=5, replace = FALSE)
  random.sp <- virtualspecies::generateRandomSp(envData[[myRandNum]], 
                                                convert.to.PA = FALSE, 
                                                species.type = "additive",
                                                realistic.sp = TRUE, 
                                                plot = FALSE)
  
  #reclassify suitability raster using a probability conversion rule
  new.pres<-convertToPA(x=random.sp,
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
  
  myPres=presence.points$sample.points[, c( "x", "y",  "Observed")]
  coordinates(myPres)<-~x+y
  raster::crs(myPres)<-myCRS
  
  # loop on two sets of environmental dataset: one with 5 bioclimatic variables and one with 19 bioclimatic variables
  env_datasets=list(envData, envData[[random.sp$details$variables]])
  myOut_list=list()
  
  for(k in 1:length(env_datasets)){
    # k=1
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
    raster::crs(myGrid.psAbs)<-myCRS
    plot(myGrid.psAbs, add=T)
    obs_prev=round(nrow(myPres)/length(myGrid.psAbs),2)
    
    # create pseudoabs from buffer
    buf <- circles(myPres, d = BuffSize, lonlat = TRUE)
    
    # buffer OUT, constrained by convex hull
    chull <-myPres %>% st_as_sf() %>% 
      summarise( geometry = st_combine( geometry ) ) %>%
      st_convex_hull()
    chull <- as(st_geometry(chull), "Spatial")
    buf.mask <- mask(envData[[1]], rgeos::gUnaryUnion(buf@polygons),inverse=TRUE)
    buf.mask<- mask(buf.mask, chull)
    myPseudo_buf<- as.data.frame(sampleRandom(buf.mask, 
                                              size = round(length(myPres)/bgk_prev) , 
                                              xy = TRUE))
    coordinates(myPseudo_buf)<-~x+y
    
    # buffer IN
    buf.mask <- mask(envData[[1]], rgeos::gUnaryUnion(buf@polygons),inverse=FALSE)
    innerBuf=as.data.frame(sampleRandom(buf.mask, 
                                        size =  round(length(myPres)/bgk_prev), 
                                        xy = TRUE))
    coordinates(innerBuf)<-~x+y
    
    #create the training datasets
    b <- as(extent(envData), "SpatialPolygons")
    myPseudo_rand<- spsample(b, n =  round(length(myPres)/bgk_prev), type = "random")
    myPseudo_rand=SpatialPointsDataFrame(myPseudo_rand,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    
    myPseudo_buf=SpatialPointsDataFrame(myPseudo_buf,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    innerBuf=SpatialPointsDataFrame(innerBuf,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    myGrid.psAbs$Observed <- rep(0, nrow(myGrid.psAbs@coords))
    
    proj4string(myPseudo_rand)=raster::crs(myPres)
    proj4string(myPseudo_buf)=raster::crs(myPres)
    proj4string(innerBuf)=raster::crs(myPres)
    proj4string(myGrid.psAbs)=raster::crs(myPres)
    
    myPAlist=list(
      "mytrain_rand"=myPres+myPseudo_rand, 
      "mytrain_buf_out"=myPres+myPseudo_buf,
      "mytrain_buf_in"=myPres+innerBuf,
      "mytrain_grid" = myPres+myGrid.psAbs)
    
    #dataset for predictions
    mypredstack=env_datasets[[k]]
    
    #count proportion of correct absences
    prop_corr_abs=lapply(myPAlist, function(x){y=subset(x, Observed==0); z=raster::extract(new.pres$pa.raster, y); 
    return(length(which(z==0))/length(z))})
    
    #modelling
    myOut_tmp=data.frame() 
    for(i in 1:length(myPAlist)){
      # i=1
      cat(paste("\nModel ",i,"/",length(myPAlist),"\n"))
      
      #extract enviromental data
      d <- sdmData(formula=Observed~., train=myPAlist[[i]], predictors=mypredstack)
      
      #model          
      m1 <- sdm(Observed~.,data=d,
                methods=c('glm', 'rf' ,'gam' , 'brt', 'Maxent'),
                replication='sub',test.percent=percTesting,n=Nreps, parallelSettings=list("parallel", ncore=myCores))
     
      #store info of the models statistics 
      p<-length(which(myPAlist[[i]]$Observed==1))
      pseabs<-length(which(myPAlist[[i]]$Observed==0))
      tmp_out=setNames(data.frame(matrix(ncol = 21, nrow = 0)), c("Species", "trainig_set", "exp_prevalence",  "obs_prevalence", "pres",'abs', "prop_corr_abs", "nDrivingFactors", 
                                                                  "ModelType", "nPredictorsPCA", "nPredictorsModel", 
                                                                  "AUC", "COR", "Deviance", "TSS", "Kappa", 
                                                                  "Sensitivity", "Specificity", "BoyceI", 'RMSE',"R2"   ))
      mynames=names(tmp_out)
      
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
                                                                                                    window.w="default", res=100, PEplot = FALSE)$Spearman.cor) 
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
                                                                                                    window.w="default", res=100, PEplot = FALSE)$Spearman.cor) 
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
                                                                                                  window.w="default", res=100, PEplot = FALSE)$Spearman.cor) 
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
                                                                                                    window.w="default", res=100, PEplot = FALSE)$Spearman.cor) 
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
                                                                                                          window.w="default", res=100, PEplot = FALSE)$Spearman.cor) 
                            })), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){
                              Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                            mean(unlist(lapply(m1@models$Observed$maxent, function(x){
                              round(cor(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)^2,3)})), na.rm=T)))
      } else {maxent_rows=c(myVs, names(myPAlist[i]),bgk_prev,obs_prev,p,pseabs, prop_corr_abs[[i]], length(random.sp$details$variables), "Maxent", nlayers(env_datasets[[k]]), nlayers(env_datasets[[k]]), rep_len(NA,10))}
      
      
      
      tmp_out=rbind(tmp_out, glm_rows, gam_rows, rf_rows, brt_rows, maxent_rows)
      names(tmp_out)=mynames
      myOut_tmp=rbind(myOut_tmp, tmp_out)
    }
    #store myOut_tmp for different environmental datasets used for modelling
    myOut_list[[k]]=myOut_tmp
  }
  #combine list and store it by species
  myOut_list=do.call(rbind, myOut_list)
  myVirtualSP_list[[myVs]]=myOut_list
}

length(myVirtualSP_list)

myVirtualSP_list[[1]]

outname=paste0("simulations/",nVirtspecies, "vs_prev", sub("\\.", "", as.character(bgk_prev)), "_occOnly_radius", BuffSize/1000, "km" , Sys.Date(), ".RDS")
saveRDS(myVirtualSP_list, outname)
