
#---- 1. Load libraries and external functions -----
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
require(hypervolume)
library(tidyverse)
library(USE)

#---- 2. Download bioclimatic rasters -----
#download bioclimatic raster
Worldclim <- raster::getData('worldclim', var='bio', res=10) 
envData <-crop(Worldclim, extent(-12, 25, 36, 60))
envData
rm(Worldclim)
rpc <-rasterPCA(envData,spca = TRUE)
dt <- na.omit(data.table(as.data.frame(rpc$map[[c("PC1", "PC2")]], xy = TRUE)))
dt <- st_as_sf(dt, coords = c("PC1", "PC2"))
myRes <- USE::optimRes(sdf=dt, grid.res=c(1:12), perc.thr = 10, showOpt = TRUE, cr = 4)


#---- 3. Generate virtual species, create background points using different methodologies, SDMs exercise and models statistics -----
nVirtspecies <- 10
percTesting <- 30
Nreps <- 5
myFold <- 5
myCores <- 5
BuffSize <- 50000 #in meters
myCRS <-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
Vsprev <- 0.25 #species geographic prevalence 
bgk_prev <- 1 #sample prevalence 

#run the whole framework
myVirtualSP_list <- list()

for(myVs in 1:nVirtspecies){
  # myVs=1
  message(paste0('\nSpecies ',myVs,'\n'))
  #create virtual species
  myRandNum=sample(1:19,size=5, replace = FALSE)
  #widespread species, at the moment 
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
  
  #Sample true occurrences
  presence.points <- sampleOccurrences(new.pres,
                                       n = 300, 
                                       type = "presence only",
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  
  plot(new.pres$pa.raster, main=paste("species #",   myVs))
  
  myPres=presence.points$sample.points[, c( "x", "y",  "Observed")]
  coordinates(myPres)<-~x+y
  raster::crs(myPres)<-myCRS

  #loop on two sets of environmental dataset
  env_datasets=list(envData, envData[[random.sp$details$variables]])
  myOut_list=list()
  
  for(k in 1:length(env_datasets)){
    set.seed(666)
    
    #grid.psAbs_ sampling in grid
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
    
    #create pseudoabs from buffer
    buf <- circles(myPres, d = BuffSize, lonlat = TRUE)
    
    #buffer OUT, constrained by convex hull
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
    
    #buffer IN
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
      "mytrain_grid" = myPres+myGrid.psAbs 
    )
    
    #dataset for predictions
    mypredstack=env_datasets[[k]]
    
    # count proportion of correct absences
    prop_corr_abs=lapply(myPAlist, function(x){y=subset(x, Observed==0); z=raster::extract(new.pres$pa.raster, y); 
    return(length(which(z==0))/length(z))})
    
    #modelling
    myOut_tmp=data.frame() 
    
    for(i in 1:length(myPAlist)){
      #  i=1
      cat(paste("\nSampling strategy ",i,"/",length(myPAlist),"\n"))
      df=na.omit(data.frame(pa=myPAlist[[i]]$Observed,terra::extract(y=data.frame(myPAlist[[i]])[,2:3],x=rast(mypredstack))[,-1]))
      PCAenv <- data.frame(pa=df$pa,princomp(scale(df[,-1]))$scores[,1:2])
      hypAbs <- hypervolume(PCAenv[PCAenv$pa==0,-1],verbose=FALSE)
      hypPres <- hypervolume(PCAenv[PCAenv$pa==1,-1],verbose=FALSE)
      ovrlp = get_volume(hypervolume_set(hypPres,hypAbs, check.memory=FALSE,verbose=FALSE))[[3]] #3 is intersection
      
      #store info of the models statistics 
      p<-length(which(myPAlist[[i]]$Observed==1))
      pseabs<-length(which(myPAlist[[i]]$Observed==0))
      tmp_out=data.frame(matrix(ncol = 10, nrow = 1,
                                dimnames = list(NULL,c("Species", "training_set", "exp_prevalence",  "obs_prevalence", "pres",'abs', "prop_corr_abs", 
                                                                  "nDrivingFactors","nPredictorsPCA", "Overlap"))))
      mynames=names(tmp_out)
      tmp_out[1,]=c(myVs, names(myPAlist[i]),bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], 
              length(random.sp$details$variables),nlayers(env_datasets[[k]]),ovrlp)
              
     
      myOut_tmp=rbind(myOut_tmp, tmp_out)
    }
    #store myOut_tmp for different environmental datasets used for modelling
    myOut_list[[k]]=myOut_tmp
  }
  #combine list and store it by species
  myOut_list=do.call(rbind, myOut_list)
  myVirtualSP_list[[myVs]]=myOut_list
}

outname=paste0(nVirtspecies, "vs_prev", sub("\\.", "", as.character(bgk_prev)), "_occOnly_radius_NicheOverlap_", BuffSize/1000, "km_" , Sys.Date(), ".RDS")
saveRDS(myVirtualSP_list, outname)

#---- 4. Plotting the results----
require(hypervolume)
library(tidyverse)

myVirtualSP_list <- readRDS("10vs_prev1_occOnly_radius_NicheOverlap_newKernel_50km_2023-06-14.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

tmp <-  myProc_out %>%
  as_tibble() %>%
  dplyr::select(training_set, exp_prevalence,
                nPredictorsPCA, Overlap) %>% 
    filter(training_set!= "mytrain_buf_in") %>% 
    mutate(training_set=as.factor(training_set),
         exp_prevalence=as.factor(exp_prevalence),
         nPredictorsModel=as.factor(nPredictorsPCA),
         Overlap=round(as.numeric(gsub(",", "", Overlap)),3), 
         trainig_set=factor(training_set, levels= c("mytrain_grid", "mytrain_rand",  "mytrain_buf_out")),
        new_name_tr_set=recode(training_set,  mytrain_grid= "Uniform", mytrain_rand= "Random", mytrain_buf_out="BufferOut"),
        new_name_tr_set=factor(new_name_tr_set,levels = c("Uniform", 'Random',"BufferIn", "BufferOut")))

tmp
summary(tmp)

Npred <- 5
mySpecies <- 10
myPrev <- 1
mytitle <- paste0("N. species = ",mySpecies, "; Prevalence = ", myPrev, "; N. predictors= ", Npred)

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
unique(tmp$new_name_tr_set)[1]

p <- tmp %>%
  filter(nPredictorsModel == Npred) %>%
  ggplot(aes(y=Overlap, x= new_name_tr_set, color=new_name_tr_set))+
  geom_boxplot()+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "BufferOut" ),
                     values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="",  y="Overlap", color="Training datasets")+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme_light()+
  theme(aspect.ratio = 1.5, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), 
        legend.title = element_text(size=14),
        legend.key.size = unit(1.5, 'cm'))
p
outname<- paste(outdir, "classOverlap_boxplot", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 16, height = 8, device='png', dpi=320)

#---- 5. ANOVA ----
library(car)
Npred<-5
aov.df <- tmp %>%
  filter(nPredictorsModel == Npred)

overlap.aov <- aov(data = aov.df, Overlap~new_name_tr_set)
car::qqPlot(resid(overlap.aov))
summary(overlap.aov)
TukeyHSD(overlap.aov, conf.level=.95)


