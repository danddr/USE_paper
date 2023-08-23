setwd("/home/ddare/working_files/SDM_pseudoAbsence/")
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
library(Hmisc)

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

envData <- raster::stack(envData) # convert to raster object, otherwise virtualspecies R pkg doesn't work

#run the whole framework
myVirtualSP_list <- list()

for(myVs in 1:nVirtspecies){
  # myVs=1
  message(paste0('\nSpecies ',myVs,'\n'))
  #create virtual species
  myRandNum <- sample(1:19,size=5, replace = FALSE)
  #widespread species, at the moment 
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
  
  #Sample true occurrences
  presence.points <- sampleOccurrences(new.pres,
                                       n = 300, 
                                       type = "presence only",
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  
  plot(new.pres$pa.raster, main=paste("species #",   myVs))
  
  myPres <- presence.points$sample.points[, c( "x", "y",  "Observed")]
  coordinates(myPres) <- ~x+y
  raster::crs(myPres) <- myCRS

  #loop on two sets of environmental dataset
  env_datasets <- list(envData, envData[[random.sp$details$variables]])
  myOut_list <- list()
  
  for(k in 1:length(env_datasets)){
    #k=2
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
    raster::crs(myGrid.psAbs) <- myCRS
    plot(myGrid.psAbs, add=T)
    obs_prev <- round(nrow(myPres)/length(myGrid.psAbs),2)
    
    #create pseudoabs from buffer
    buf <- circles(myPres, d = BuffSize, lonlat = TRUE)
    
    #buffer OUT, constrained by convex hull
    chull <- myPres %>% st_as_sf() %>% 
      summarise( geometry = st_combine( geometry ) ) %>%
      st_convex_hull()
    
    chull <- as(st_geometry(chull), "Spatial")
    buf.mask <- raster::mask(envData[[1]], rgeos::gUnaryUnion(buf@polygons), inverse=TRUE)
    buf.mask <- raster::mask(buf.mask, chull)
    myPseudo_buf <- as.data.frame(sampleRandom(buf.mask, 
                                              size = round(length(myPres)/bgk_prev) ,
                                              xy = TRUE))
    coordinates(myPseudo_buf)<-~x+y
   
    #create the training datasets
    b <- as(extent(envData), "SpatialPolygons")
    myPseudo_rand <- spsample(b, n =  round(length(myPres)/bgk_prev), type = "random")
    myPseudo_rand <- SpatialPointsDataFrame(myPseudo_rand,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    myPseudo_buf <- SpatialPointsDataFrame(myPseudo_buf,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
    myGrid.psAbs$Observed <- rep(0, nrow(myGrid.psAbs@coords))
    
    proj4string(myPseudo_rand) <- raster::crs(myPres)
    proj4string(myPseudo_buf) <- raster::crs(myPres)
    proj4string(myGrid.psAbs) <- raster::crs(myPres)
    
    myPAlist <- list(
      "mytrain_rand" = myPseudo_rand, 
      "mytrain_buf_out" =  myPseudo_buf,
      "mytrain_grid" =  myGrid.psAbs 
    )
    PA_name <- c( "Random", "BufferOut", "Uniform")
    
    #dataset for predictions
    mypredstack<-env_datasets[[k]]
    
    # count proportion of correct absences
    prop_corr_abs<-lapply(myPAlist, function(x){y=subset(x, Observed==0); z=raster::extract(new.pres$pa.raster, y); 
    return(length(which(z==0))/length(z))})
    
    #modelling
    rpc <- USE::rastPCA(env.rast = mypredstack, nPC = 2, stand = TRUE)
    pcstack <- c(rpc$PCs$PC1, rpc$PCs$PC2)
    
    myOut_tmp<-data.frame() 
    
    for(i in 1:length(myPAlist)){
       # i=1
      tmp<-terra::vect(myPAlist[[i]])
      tmp<-na.omit(terra::extract(pcstack, tmp, df=TRUE))
      tmp$group<-PA_name[[i]]
     
      rangePC1 <- round(range(tmp$PC1)[2] - range(tmp$PC1)[1],3)
      rangePC2 <- round(range(tmp$PC2)[2] - range(tmp$PC2)[1],3)
       
      #store info of the models statistics 
      p <- length(which(myPAlist[[i]]$Observed==1))
      pseabs <- length(which(myPAlist[[i]]$Observed==0))
     
      tmp_out <- data.frame(matrix(ncol = 11, nrow = 1,
                                dimnames = list(NULL,c("Species", "training_set", "exp_prevalence",  "obs_prevalence", "pres",'abs', "prop_corr_abs", 
                                                       "nDrivingFactors","nPredictorsPCA", "range", "PCs"))))
      mynames <- names(tmp_out)
     
      tmp_out[1:2,] <- rbind(
        c(myVs, PA_name[i],bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], 
          length(random.sp$details$variables),nlayers(env_datasets[[k]]),  rangePC1, 1) , 
        c(myVs, PA_name[i],bgk_prev, obs_prev,p,pseabs, prop_corr_abs[[i]], 
          length(random.sp$details$variables),nlayers(env_datasets[[k]]),  rangePC2, 2)
      )
   
      myOut_tmp<-rbind(myOut_tmp, tmp_out)
    }
    
    #store myOut_tmp for different environmental datasets used for modelling
    myOut_list[[k]]<-myOut_tmp
  }
  #combine list and store it by species
  myOut_list=do.call(rbind, myOut_list)
  myVirtualSP_list[[myVs]]=myOut_list
}

outdir<-"YOUR_OUTDIR"
outname<-paste0(outdir, nVirtspecies, "vs_prev", sub("\\.", "", as.character(bgk_prev)), "_PCs_ranges_", BuffSize/1000, "km_" , Sys.Date(), ".RDS")
saveRDS(myVirtualSP_list, outname)

#---- 4. Plotting the results----
library(tidyverse)
myVirtualSP_list <- readRDS( "10vs_prev1_PCs_ranges_50km_2023-08-09.RDS")
myProc_out <- do.call(rbind, myVirtualSP_list)

myProc_out <- myProc_out %>%
  as_tibble() %>%
  dplyr::select(training_set, exp_prevalence,
                nPredictorsPCA, range,  PCs) %>%   
  mutate(training_set = as.factor(training_set),
         exp_prevalence = as.factor(exp_prevalence),
         nPredictorsModel = as.factor(nPredictorsPCA),
         range = as.numeric(range), 
         # PCs = as.numeric(PCs), 
         PCs = ifelse(PCs=="1", "PC1", "PC2"),
         training_set = dplyr::recode(training_set,
                                 BufferOut = "Buffer-out"),
         training_set = factor(training_set,levels = c("Uniform", 'Random', "Buffer-out")))

unique_combo <- unique(myProc_out[, c("training_set", "PCs")])

rangeOut <- as.data.frame(do.call(rbind, Map(function(x,y){
  df <- myProc_out[myProc_out$training_set == x & myProc_out$PCs == y, ]
  df$range <- as.numeric(df$range)
  ourRange <- Hmisc::smean.cl.boot(df$range, B=2000) 
}, x=unique_combo$training_set, y=unique_combo$PCs)
))

rangeOut$PCs <- unique_combo$PCs
rangeOut$training_set <- unique_combo$training_set

# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- rangeOut %>%
  ggplot(aes(y=Mean, x= training_set, color=training_set))+
  geom_point()+
  geom_errorbar(aes(ymin=Lower, ymax=Upper))+
  scale_color_manual(breaks = c("True Absences", "Uniform", "Random",  "BufferIn", "Buffer-out" ),
  values=c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7" ))+
  labs(x="",  y="Total range (95% CI)", color="Training datasets")+
  # scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  facet_wrap(~PCs)+
  theme_light()+
  theme(aspect.ratio = 1.5, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), 
        legend.title = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))
p

outname<- paste("figures_paper/June2023/", "locatioBiasRange_plot", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 16, height = 8, device='png', dpi=320)

#---- 5. ANOVA ----
library(car)
# Npred<-5

with(subset(myProc_out, PCs == "PC1" ), leveneTest(range, training_set))
#cannot use ANOVA

library(dunn.test)
myProc_out %>% 
  dplyr::select(training_set, PCs, range) %>% 
  split(.$PCs) %>%   
  map(~dunn.test(x = .$range, g =.$training_set, altp=TRUE, method = "holm",list=TRUE,wrap=TRUE)) 


