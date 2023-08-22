library(USE)
library(raster)
library(virtualspecies)
library(tidyverse)
library(ggplot2)
library(dismo)
library(cowplot)
library(remotes)
library(viridis)
library(BuenColors)
library(RColorBrewer)
source("color_density.R")

cl <- colorRampPalette(c("#3E49BB", "#3498DB", "yellow", "orange", "red", "darkred"))
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#---- 1. download bioclimatic raster ----
Worldclim <- raster::getData('worldclim', var='bio', res=10) #Valid resolutions are 0.5, 2.5, 5, and 10 (minutes of a degree). In the case of res=0.5, you must also provide a lon and lat argument for a tile
my.stack <- crop(Worldclim, extent(-12, 25, 36, 60))
b <- as(extent(my.stack), "SpatialPolygons")
myCRS <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
Vsprev <- 0.25
bgk_prev <- 1
myRes <- list()
myRes$Opt_res <- 6
BuffSize <- 50000

#create virtual species
set.seed(666)
myRandNum <- sample(1:19,size=5, replace = FALSE)
envData <- my.stack[[myRandNum]]

#add environmental density background
rpc <- USE::rastPCA(env.rast = envData, nPC = 2, stand = TRUE)
pcstack <- c(rpc$PC1, rpc$PC2)

#plot pca in the environmental space
env_pca <- na.omit(as.data.frame(pcstack))
env_pca$density <- get_density(env_pca$PC1, env_pca$PC2, n = 100)
env_pca %>% 
  ggplot( aes(PC1, PC2))+
  geom_point(data=env_pca, aes(PC1, PC2, color = density), alpha=0.7)+
  scale_color_viridis()+
  labs(color="Density of PC-scores")+
  # xlim(-16,10)+ylim(-10,10)+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14),  legend.text=element_text(size=10))

#---- 2. Create virtual species ----
random.sp <- virtualspecies::generateRandomSp(envData, 
                                              convert.to.PA = FALSE, 
                                              species.type = "additive",
                                              realistic.sp = TRUE, 
                                              plot = FALSE)

#reclassify suitability raster using a probability conversion rule
new.pres <- convertToPA(x=random.sp, 
                      # beta=0.55,
                      beta="random",
                      alpha = -0.05, plot = FALSE, 
                      species.prevalence = Vsprev) 

#Sample true occurrences
presence.points <- sampleOccurrences(new.pres,
                                     n = 300, # The number of points to sample
                                     type = "presence only",
                                     plot = TRUE)
myPres <- presence.points$sample.points[, c( "x", "y",  "Observed")]
coordinates(myPres) <- ~x+y
raster::crs(myPres) <- myCRS

#---- 3. Create pseudo-absences datasets ----
myGrid.psAbs <- USE::paSampling(envData,
                                pres=myPres,
                                thres=0.75,
                                H=NULL,
                                grid.res=myRes$Opt_res,
                                prev=bgk_prev,
                                sub.ts=FALSE, 
                                plot_proc=FALSE)

pcscores <- data.frame(st_coordinates(myGrid.psAbs))
class(pcscores)
names(pcscores) <- c("PC1", "PC2")

#only keep coords
myGrid.psAbs <- data.frame(x = myGrid.psAbs$x, y = myGrid.psAbs$y)
coordinates(myGrid.psAbs) <- ~x+y
raster::crs(myGrid.psAbs)<-myCRS
obs_prev=round(nrow(myPres)/length(myGrid.psAbs),2)

# create pseudoabs from buffer
buf <- circles(myPres, d = BuffSize, lonlat = TRUE)

# buffer OUT, constrained by convex hull
chull <-myPres %>% st_as_sf() %>% 
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull()
chull <- as(st_geometry(chull), "Spatial")
buf.mask <- mask(envData[[1]], rgeos::gUnaryUnion(buf@polygons),inverse=TRUE)
buf.mask <- mask(buf.mask, chull)
myPseudo_buf <- as.data.frame(sampleRandom(buf.mask, 
                                          size = round(length(myPres)/bgk_prev) , #length(myGrid.psAbs)
                                          xy = TRUE))
coordinates(myPseudo_buf)<-~x+y

# random
myPseudo_rand <- sampleRandom(envData$bio11, size =  round(length(myPres)/bgk_prev), sp=TRUE)
names(myPseudo_rand) <-"Observed"
myPseudo_rand$Observed <-0

#create the training datasets
myPseudo_buf <- SpatialPointsDataFrame(myPseudo_buf,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
innerBuf <- SpatialPointsDataFrame(innerBuf,data.frame("Observed"=rep(0,  round(length(myPres)/bgk_prev))))
myGrid.psAbs$Observed <- rep(0, nrow(myGrid.psAbs@coords))

proj4string(myPseudo_rand) <- raster::crs(myPres)
proj4string(myPseudo_buf) <- raster::crs(myPres)
proj4string(myGrid.psAbs) <- raster::crs(myPres)

my_meth <- list( myPseudo_rand, innerBuf, myPseudo_buf)
meth_name <- c( "Random", "BufferOut")

#---- 4. extract PC_scores from the geographic space ----
pcstack
myOut_list <- list()
for(i in 1:length(my_meth)){
  # i=1
  tmp=subset(my_meth[[i]], Observed==0)
  tmp=raster::extract(pcstack, tmp, df=TRUE)
  tmp$group=meth_name[[i]]
  myOut_list[[i]]=tmp
  rm(tmp)  
}

myOut_list <- do.call(rbind.data.frame, myOut_list)
pcscores$group <- "Uniform"

myOut_list <- myOut_list %>%
  dplyr::select(PC1, PC2, group) %>% 
  drop_na() %>%
  bind_rows(pcscores)

names(myOut_list) 
myOut_list$group <- factor(myOut_list$group, levels = c( "Uniform", "Random", "BufferOut"))

#---- 5. plot bivariate geom dens----
myOut_list %>% 
  ggplot( aes(PC1, PC2))+
  geom_density_2d(col= "darkgray")+
  xlim(-16,10)+ylim(-10,10)+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  facet_wrap(~group)+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14))

#add environmental density background
p <- myOut_list %>% 
  ggplot( aes(PC1, PC2))+
  geom_point(data=env_pca, aes(PC1, PC2, color = density), alpha=0.7)+
  scale_color_viridis()+
  geom_density_2d(col= "red1", bins=20, linewidth=0.5)+
  labs(color="Density of PC-scores")+
  xlim(-6,6)+ylim(-6,6)+
  facet_wrap(~group, nrow=4)+
  theme_light()+
  theme(plot.title = element_text(size=14,face = 'bold'),
        aspect.ratio = 1,
        legend.background=element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size=12,face = 'bold'),
        legend.position = 'none',
         legend.text = element_text(size=9,angle = 30))+
  guides(color = guide_colorbar(title.position="top", title.vjust =0.5))

p

p2 <- myOut_list %>% 
  pivot_longer(c("PC1", "PC2")) %>% 
  drop_na() %>% 
  ggplot(aes(value, fill=name))+
  geom_histogram(alpha=0.6)+
  ylim(0,60)+
  labs(x="PC-scores", y="Frequency", fill="")+
  facet_grid(group~name)+ #
  theme_light()+
  theme(plot.title = element_text(size=14,face = 'bold'),
        legend.background=element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1,
        strip.text.x = element_text(size=12,face = 'bold'),
        strip.text.y = element_text(size=12,face = 'bold'), 
        legend.position = 'right')+
        guides(color = guide_colorbar(title.position="top", title.vjust =0.5))

p2

#combine plots
pp <-plot_grid(p,p2, ncol=2,labels = "AUTO")
outname <- paste("bivariate_dens_plot_UESampling_", Sys.Date(),".pdf", sep="")
ggsave(pp, filename = outname, width = 16, height = 8, device='pdf', dpi=320)

#suppl mat
p<-myOut_list %>% 
  filter( group=="Uniform"  ) %>%
  ggplot( aes(PC1, PC2))+
  geom_point(data=env_pca, aes(PC1, PC2, color = density), alpha=0.7)+
  scale_color_viridis()+
  geom_density_2d(col= "red1", bins=10, size=0.5)+
  labs(color="Density of PC-scores")+
  xlim(-6,6)+ylim(-6,6)+
  facet_wrap(~group, nrow=4)+
  theme_light()+
  theme(plot.title = element_text(size=14,face = 'bold'),
        aspect.ratio = 1,
        legend.background=element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size=12,face = 'bold'),
        legend.position = 'none',
        legend.text = element_text(size=9,angle = 30))+
  guides(color = guide_colorbar(title.position="top", title.vjust =0.5))

p