# Code for the analyses of Fagus sylvatica: an example about how to use the USE package to uniformly sample
# pseudo-absences within the environmental space
# Author: Manuele Bazzichetto
# Date: June, 16th, 2023
# For inquiries about the analyses presented below contact manuele.bazzichetto@gmail.com

#R libraries
#handle spatial data
library(sf)
library(terra)
library(raster)
#library(geodata) #as an alternative to using raster::getData to download climatic data
#plot data
library(mapview)
library(ggplot2)
library(ggpubr)
library(tmap) #"World" shapefile
library(GGally)
#modelling tools
library(car)
library(caret) #findCorrelation
library(performance) #performance function
library(ranger)
#install USE
library(devtools)
library(USE)
#working in parallel
library(parallel)

#Some important info:

#1) to-install USE from GitHub repo - last instal.: June, 16th 2023 (version 0.1.1.)
#devtools::install_github("danddr/USE")

#2) data for replicating the analyses can be found at the following sources:
#WorldClim -> automatically doewnloaded using geodata::worldclim_global
#sPlotOpen -> https://idata.idiv.de/ddm/Data/ShowData/3474?version=54
#EU-Forest -> https://figshare.com/articles/dataset/Occurrences_location_shapefile/3497891?backTo=/collections/A_high-resolution_pan-European_tree_occurrence_dataset/3288407

##Get climatic data -------------------------------------------------------------
#2.5 (degrees at the Equator) resolution. Date of last download: June, 16th (2023)
#BioClimateData <- geodata::worldclim_global(var = 'bio', res = 2.5, path = "~/Documents/UEsampling/CleanCode/worldclim_cleancode/") #here, put a directory where to store climatic data
BioClimateData <- raster::getData(name = "worldclim", download = T, res = 2.5, var = "bio", path = "~/Documents/UEsampling/CleanCode/worldclim_cleancode/") #last download: June, 16th (2023)

class(BioClimateData)

BioClimateData <- as(BioClimateData, "SpatRaster")

crs(BioClimateData, proj = T)

#crop and mask BioClimateData to cover the area of interest (i.e. Italy, France, Spain)
#bbox used for cropping World allows first getting rid of France's out of seas territories (e.g., French Guiana)
data("World")
st_crs(World) #EPSG:4326

#subset for IFS -> check as there's still some line within the (dissolved) polygon
IFS_poly <- st_crop(st_union(World[World$name %in% c("Italy", "France", "Spain"), ]),
                    xmin = -9.264307, ymin = 35.507217, xmax = 18.439300, ymax = 50.995800)

st_crs(BioClimateData); st_crs(IFS_poly)

#crop and mask BioClimateData to get BioClimateData.Eu
BioClimateData.Eu <- terra::crop(BioClimateData, ext(st_bbox(IFS_poly)))

BioClimateData.Eu <- terra::mask(BioClimateData.Eu, vect(IFS_poly))

#Bioclimatic variables from Bio1 to Bio11 must be scaled by 10 or 100 (Bio3 & Bio4)

for(i in seq_len(11)) {
  if(i %in% c(3, 4)) {
    BioClimateData.Eu[[i]] <- BioClimateData.Eu[[i]]/100
  } else {
    BioClimateData.Eu[[i]] <- BioClimateData.Eu[[i]]/10
  }
}

#get rid of BioClimateData
rm(BioClimateData)

plot(BioClimateData.Eu$bio1)

## Find optimal resolution for sampling European environmental (here climatic) space ----
PCA_Eu <- USE::rastPCA(env.rast = BioClimateData.Eu, nPC = 2, stand = TRUE)

PCstack <- c(PCA_Eu$PCs$PC1, PCA_Eu$PCs$PC2)

#Get dataframe from PCstack
PCstack.df <- as.data.frame(PCstack, xy = T, na.rm = T)

#make the PCstack.df spatial
PCstack.sp <- st_as_sf(PCstack.df, coords = c("PC1", "PC2"))

#set the clusters for parallelised version of USE::optimRes
cr7 <- makeCluster(7)
#clusterExport is to export variables to the clusters
clusterExport(cr7, "PCstack.sp")
#clusterEvalQ is to export packages to the clusters
clusterEvalQ(cr7, library(sf))

#get optimal resolution of the grid for sampling European environmental space
#stopClusters is called internally
Optres_Eu <- optimRes(sdf = PCstack.sp, grid.res = seq(1, 15), cr = cr7, showOpt = T)

##Get Fagus sylvatica data from sPlot ------------------------------------------ 
#note that these data will be used as a supplementary testing dataset to assess the
#predictive performance of modelling algorithms (GLM and Random Forest) calibrated on Fagus sylvatica occurrence (presences/pseudo-absences) data
#sampled using the uniform environmental sampling approach

#import sPlotOpen matrices from sPlotOpen.RData
#all matrices + an R function for automatically generating a list of data references will be imported
load("~/Documents/UEsampling/sPlotOpen.RData") #use the directory where sPlotOpen is stored

#these are all "tbl_df", "tbl", "data.frame" -> coerce to data.frame
vapply(list(header.oa, DT2.oa, metadata.oa, CWM_CWV.oa, reference.oa), class, FUN.VALUE = character(3))

header.oa <- as.data.frame(header.oa)
DT2.oa <- as.data.frame(DT2.oa)
metadata.oa <- as.data.frame(metadata.oa)
CWM_CWV.oa <- as.data.frame(CWM_CWV.oa)
reference.oa <- as.data.frame(reference.oa)

#check
sapply(list(header.oa, DT2.oa, metadata.oa, CWM_CWV.oa, reference.oa), class)

#PlotObservationID is unique
sum(duplicated(header.oa$PlotObservationID)) #good

#keep records in Italy, Spain and France for Resample 1
c("Italy", "France", "Spain") %in% header.oa$Country #TRUE*3
class(header.oa$Resample_1)

#select data from Italy, France and Spain and belonging to the resample #1
header.IFS.R1 <- header.oa[(header.oa$Country %in% c("Italy", "France", "Spain")) & header.oa$Resample_1, ]

#make the selection spatial
header.IFS.R1.sp <- st_as_sf(header.IFS.R1, coords = c("Longitude", "Latitude"))

#set CRS
st_crs(header.IFS.R1.sp) <- 4326

#count number of records for Fagus sylvatica in DT2.oa
sum(DT2.oa$Species %in% "Fagus sylvatica") #3761

#get data on Fagus sylvatica
Fagus_DT2 <- DT2.oa[which(DT2.oa$Species %in% "Fagus sylvatica"), ]

#count number of observations for Fagus sylvatica in Italy, France, Spain
nrow(header.IFS.R1[header.IFS.R1$PlotObservationID %in% Fagus_DT2$PlotObservationID, ]) #392

#create a dataframe with sPlot data for Fagus sylvatica from Italy, France and Spain
Fagus.sPlot <- header.IFS.R1[header.IFS.R1$PlotObservationID %in% Fagus_DT2$PlotObservationID, ]

#add a PA (presence/absence) column and get rid of information (e.g., columns) that will not be used
Fagus.sPlot$PA <- 1
Fagus.sPlot <- Fagus.sPlot[c("PA", "Longitude", "Latitude")]

#retrieve absence data
Fagus.abs.sPlot <- header.IFS.R1[!header.IFS.R1$PlotObservationID %in% Fagus_DT2$PlotObservationID, ]
Fagus.abs.sPlot$PA <- 0
Fagus.abs.sPlot <- Fagus.abs.sPlot[c("PA", "Longitude", "Latitude")]

#put presence/absence data together
Fagus.sPlot <- rbind(Fagus.sPlot, Fagus.abs.sPlot)
colnames(Fagus.sPlot)[-1] <- c("long", "lat") 

#get bioclimatic values for presence/absence data
Fagus.sPlot <- data.frame(Fagus.sPlot, terra::extract(BioClimateData.Eu, Fagus.sPlot[c("long", "lat")]))

#get rid of NAs
Fagus.sPlot <- na.omit(Fagus.sPlot)

table(Fagus.sPlot$PA) #367 presences, 4162 absences

##Get Fagus sylvatica data from EU-forest --------------------------------------

#import EU-forest data
EU_forest <- read.csv(file = "~/Documents/UEsampling/EU_forest_occ/EUForestspecies.csv", sep = ",",
                      header = T, stringsAsFactors = F) #use the directory where EUForest.csv is stored

#check NAs, but note that for some variables NA is coded as -9.9999
anyNA(EU_forest) #F

unique(EU_forest$COUNTRY)

"Fagus sylvatica" %in% EU_forest$SPECIES.NAME

#make sure there's only Fagus sylvatica
unique(grep(pattern = "Fagus", x = EU_forest$SPECIES.NAME, value = T)) #there's only Fagus sylvatica

#geographic subset -> this will only include data for Fagus sylvatica from Italy, France and Spain
EU_forest.IFS.Fag <- EU_forest[EU_forest$COUNTRY %in% c("Italy", "France", "Spain") & (EU_forest$SPECIES.NAME %in% "Fagus sylvatica"), ]

#remove EU-forest to save memory
rm(EU_forest)

#make the subset spatial
EU_forest.IFS.Fag.sp <- st_as_sf(EU_forest.IFS.Fag, coords = c("X", "Y"))

#EPSG from metadata EU-Forest seems wrong -> using EPSG:3035 (ETRS89-extended LAEA) 
#set CRS
st_crs(EU_forest.IFS.Fag.sp) <- 3035

#have a look at the data
mapview(x = EU_forest.IFS.Fag.sp, zcol = "COUNTRY")

#transform spatial subset to EPSG:4326
EU_F_IFS.Fag.geo <- st_transform(EU_forest.IFS.Fag.sp, crs = 4326)

#have a look at the data
mapview(raster(BioClimateData.Eu[[1]])) + mapview(EU_F_IFS.Fag.geo, zcol = "COUNTRY")

#there are 12.444 obs for Fagus sylvatica (presences)
#uniformSampling (from USE) will be used to subset both presences and pseudo-absences within the env. space
#the same will be done for getting testing presences and pseudo-absences

##Subset presences within the environmental space with uniformSampling --------------

#extract PC1/2 values from presence points
#PCstack was created using BioClimateData.EU
EU_Fag.PC12 <- terra::extract(PCstack, st_coordinates(EU_F_IFS.Fag.geo))

#these columns will serve as coordinates of the presences within the environmental space
EU_F_IFS.Fag.geo$PC1 <- EU_Fag.PC12[, 1]
EU_F_IFS.Fag.geo$PC2 <- EU_Fag.PC12[, 2]

#get back to a non-spatial object
EU_Fag.pres <- EU_F_IFS.Fag.geo
st_geometry(EU_Fag.pres) <- NULL

#re-attach geographical coordinates to use them later
EU_Fag.pres <- data.frame(EU_Fag.pres, st_coordinates(EU_F_IFS.Fag.geo))

#count rows with NAs
length(unique(which(is.na(EU_Fag.pres), arr.ind = T)[, 1])) #171

#check which observation have NAs for PC1
mapview(raster(PCstack[[1]])) + mapview(EU_F_IFS.Fag.geo[which(is.na(EU_F_IFS.Fag.geo$PC1)), ]) #+ mapview(EU_Fag.pres)

#assign new geometry column
EU_Fag.pres <- st_as_sf(na.omit(EU_Fag.pres), coords = c("PC1", "PC2"))

#subset presences using uniformSampling
#100 presences will be (whenever possible) sampled from each cell of the sampling grid
#same thing will be done to sample testing presences
set.seed(17534)
Fag.pres.tr.ts <- USE::uniformSampling(sdf = EU_Fag.pres, grid.res = Optres_Eu$Opt_res,
                                         n.tr = 100, sub.ts = T, n.ts = 100, plot_proc = T)

#count number of training and testing presences for country
table(Fag.pres.tr.ts$obs.tr$COUNTRY)
table(Fag.pres.tr.ts$obs.ts$COUNTRY)

#make training and testing (presence) datasets spatial 
Fag.pres.tr.ts.sp <- lapply(Fag.pres.tr.ts, function(x) {
  spdt <- x
  spdt$geometry <- NULL
  spdt.sp <- st_as_sf(spdt, coords = c("X", "Y"))
  st_crs(spdt.sp) <- 4326
  return(spdt.sp)
})

mapview(raster(PCstack[[1]])) + mapview(Fag.pres.tr.ts.sp$obs.tr, color = "blue") + mapview(Fag.pres.tr.ts.sp$obs.ts, color = "red")

##Get pseudo-absences within the environmental space with paSampling --------------

#get pseudo-absences using paSampling
#here we use all Fagus s. presences (12.444 - 171 NAs for PC1 and PC2) to inform the kernel-based filter, so that we can use all information available about locations suitable for the species
#this will allow us safely excluding as many pseudo-absences located in areas of the environmental space where conditions may actually be suitable for Fagus sylvatica (presences) as possible
#however, as we want prevalence to be = 1 in the training set, we set the prev arg as if we were using Fag.pres.tr.ts$obs.tr
#(so if we had 1.827 presences)
#basically we 'inflate' prevalence so that we can finally get a number of pseudo-absences comparable with Fag.pres.tr.ts$obs.tr

(12444-171)/(nrow(Fag.pres.tr.ts$obs.tr)) #our prevalence is 6.71757

all.equal(which(is.na(EU_F_IFS.Fag.geo$PC1)), which(is.na(EU_F_IFS.Fag.geo$PC2)))

EU_F_IFS.Fag.geo <- EU_F_IFS.Fag.geo[-which(is.na(EU_F_IFS.Fag.geo$PC1)), ]

set.seed(17893)
Fag.abs.UE <- USE::paSampling(env.rast = BioClimateData.Eu,
                                      pres = EU_F_IFS.Fag.geo, n.tr = 20,
                                      grid.res = Optres_Eu$Opt_res, prev = 6.71757,
                                      sub.ts = T, n.ts = 150, plot_proc = T)

#make training and testing (absence) datasets spatial 
Fag.abs.tr.ts.sp <- lapply(Fag.abs.UE, function(x) {
  spdt <- x
  spdt$geometry <- NULL
  spdt.sp <- st_as_sf(spdt, coords = c("x", "y"))
  st_crs(spdt.sp) <- 4326
  return(spdt.sp)
})

mapview(raster(PCstack[[1]])) + mapview(Fag.abs.tr.ts.sp$obs.tr, color = "blue") +
  mapview(Fag.abs.tr.ts.sp$obs.ts, color = "red")

#extract climatic data and join training and testing datasets
FagusEU_data <- Map(function(x, y) {
  x <- data.frame(PA = 1, st_coordinates(x))
  y <- data.frame(PA = 0, st_coordinates(y))
  colnames(x)[-1] <- colnames(y)[-1] <- c("long", "lat")
  df <- rbind(x, y)
  df <- data.frame(df, terra::extract(BioClimateData.Eu, df[c(2, 3)]))
  return(df)
}, x = Fag.pres.tr.ts.sp, y = Fag.abs.tr.ts.sp)

#length(FagusEU_data)

#set names of the list including training and testing datasets
names(FagusEU_data) <- c("Tr_data", "Ts_data")

#n° of presences and absences
table(FagusEU_data$Tr_data$PA)
table(FagusEU_data$Ts_data$PA)

#make prevalence in the testing set = 1 (i.e. reduce number of testing presences)
Subs.ts.pres <- FagusEU_data$Ts_data
Subs.ts.pres <- Subs.ts.pres[Subs.ts.pres$PA == 1, ]
Subs.ts.pres <- Subs.ts.pres[sample(x = nrow(Subs.ts.pres), size = sum(FagusEU_data$Ts_data$PA == 0), replace = F), ]

#substitute testing dataset in FagusEU_data
FagusEU_data$Ts_data <- rbind(Subs.ts.pres, FagusEU_data$Ts_data[FagusEU_data$Ts_data$PA == 0, ])

#check
table(FagusEU_data$Ts_data$PA)

#check NAs
sapply(FagusEU_data, anyNA) #F F

##Modelling (GLM + Random forest) occurrence of Fagus sylvatica ----------------
GGally::ggcorr(cor(FagusEU_data$Tr_data[-c(1, 2, 3, 4)]), geom = "text")

#retrieve highly correlated variables in training and testing datasets
#these variables will be excluded to avoid multicollinearity issues
Vars_to_remove <- caret::findCorrelation(cor(FagusEU_data$Tr_data[-c(1, 2, 3, 4)]), cutoff = .6, names = T)

#variables to keep are: 
paste0("bio", seq_len(19))[!paste0("bio", seq_len(19)) %in% Vars_to_remove]
#bio 6: min. temperature of the coldest month;
#bio 7: temperature annual range;
#bio 8: mean temperature wettest quarter.

#get rid of correlated variables from training and testing datasets
FagusEU_data <- lapply(FagusEU_data, function(x) {
  x <- x[!colnames(x) %in% Vars_to_remove]
  return(x)
})

#check distribution of values of predictors for presence and absence data
ggarrange(plotlist = lapply(colnames(FagusEU_data$Tr_data[-c(1, 2, 3, 4)]), function(nm) {
  ggplot(FagusEU_data$Tr_data, aes_string(x = nm)) +
    geom_density(aes(fill = as.factor(PA)), alpha = .3) +
    scale_fill_viridis_d(name = paste("Pres (1)", "Abs (0)", sep = "\n")) +
    theme_classic()
  }),
  nrow = 2, ncol = 2)

#check distribution of predictors
hist(FagusEU_data$Tr_data$bio6, xlab = "bio6", main = NULL)
hist(FagusEU_data$Tr_data$bio7, xlab = "bio7", main = NULL)
hist(FagusEU_data$Tr_data$bio8, xlab = "bio8", main = NULL)

#binary Generalized Linear Model
Mod_FagusEU <- glm(PA ~ bio6 + bio7 + bio8 + lat, family = binomial, data = FagusEU_data$Tr_data)

car::S(Mod_FagusEU) #Wald's test
car::Anova(Mod_FagusEU) #Likelihood ratio test
car::vif(Mod_FagusEU)
car::marginalModelPlots(Mod_FagusEU)

#some diagnostics
car::residualPlots(Mod_FagusEU) #include all poly(x, 2) terms
car::influenceIndexPlot(Mod_FagusEU)
car::outlierTest(Mod_FagusEU)

#introduce 2nd order polynomials for variables
Mod_FagusEU.2 <- glm(PA ~ poly(bio6, 2) + poly(bio7, 2) + poly(bio8, 2) +
                       poly(lat, 2), 
                     family = binomial, data = FagusEU_data$Tr_data)


car::S(Mod_FagusEU.2)

#compare nested model without poly for bio8
anova(glm(PA ~ poly(bio6, 2) + poly(bio7, 2) + bio8 + poly(lat, 2), 
          family = binomial, data = FagusEU_data$Tr_data), Mod_FagusEU.2, test = "Chisq")

#get rid of second-order term for for bio8
Mod_FagusEU.3 <- glm(PA ~ poly(bio6, 2) + poly(bio7, 2) + bio8 + poly(lat, 2), 
                     family = binomial, data = FagusEU_data$Tr_data)

car::marginalModelPlots(Mod_FagusEU.3) #fitted model seems much more adequate

#GLM predictive performance

#1) within-sample cross-validation
#set seed for reproducibility (although results are consistent across CV replicates)
set.seed(23635)
CV.GLM <- colMeans(cross_val(df = FagusEU_data$Tr_data, pa_col = "PA", formula = formula(Mod_FagusEU.3),
                   mod_type = "GLM", folds = 5L))

#2) out-of-sample prediction (OOS)
OOS.fitted <- predict(Mod_FagusEU.3, newdata = FagusEU_data$Ts_data, type = "response")
OOS.tss <- ecospat.max.tss(Pred = OOS.fitted, Sp.occ = FagusEU_data$Ts_data$PA)$max.TSS
OOS.boyce_ind <- ecospat.boyce(fit = OOS.fitted, obs = OOS.fitted[which(FagusEU_data$Ts_data$PA == 1)], nclass = 0,
                           window.w = "default", res = 100, PEplot = F)$Spearman.cor #check later details args

#3) predictive performance on sPlot data - note that Fagus.sPlot has 367 presences and 4162 absences (prevalence = 0.09)
sPlot.fitted <- predict(Mod_FagusEU.3, newdata = Fagus.sPlot, type = "response")
sPlot.tss <- ecospat.max.tss(Pred = sPlot.fitted, Sp.occ = Fagus.sPlot$PA)$max.TSS
sPlot.boyce_ind <- ecospat.boyce(fit = sPlot.fitted, obs = sPlot.fitted[which(Fagus.sPlot$PA == 1)], nclass = 0,
                                 window.w = "default", res = 100, PEplot = F)$Spearman.cor #check later details args

#goodness-of-fit
performance::model_performance(Mod_FagusEU.3) #Tjur's R2 = 0.355
(Mod_FagusEU.3$null.deviance - Mod_FagusEU.3$deviance)/Mod_FagusEU.3$null.deviance #Deviance-based  R2 0.3
#same as Deviance-based R2
modEvA::Dsquared(Mod_FagusEU.3) #0.3

# Random Forest
#for multicollinearity issues see: https://stats.stackexchange.com/questions/141619/wont-highly-correlated-variables-in-random-forest-distort-accuracy-and-feature

set.seed(45372)

RF_FagusEU <- ranger::ranger(formula = PA ~  bio6 + bio7 + bio8 + lat,
                               importance = 'permutation',
                               data = FagusEU_data$Tr_data)

ranger::importance(RF_FagusEU)

#check range of (within-sample) predictions is within the unit interval [0, 1]
range(predict(object = RF_FagusEU, data = FagusEU_data$Tr_data, type = "response")$predictions)
hist(predict(object = RF_FagusEU, data = FagusEU_data$Tr_data, type = "response")$predictions)

#RF predictive performance
#1) within-sample cross-validation
set.seed(89347)
CV.RF <- colMeans(cross_val(df = FagusEU_data$Tr_data, pa_col = "PA", formula = formula(RF_FagusEU),
                   mod_type = "RF", folds = 5L, importance = "permutation"))

#2) out-of-sample prediction (OOS)
OOS.fitted.RF <- predict(object = RF_FagusEU, data = FagusEU_data$Ts_data, type = "response")$predictions
OOS.tss.RF <- ecospat.max.tss(Pred = OOS.fitted.RF, Sp.occ = FagusEU_data$Ts_data$PA)$max.TSS
OOS.boyce_ind.RF <- ecospat.boyce(fit = OOS.fitted.RF, obs = OOS.fitted.RF[which(FagusEU_data$Ts_data$PA == 1)], nclass = 0,
                                  window.w = "default", res = 100, PEplot = F)$Spearman.cor #check later details args

#3) predictive performance on sPlot data - notice that Fagus.sPlot has 367 presences and 4162 absences (prevalence = 0.09)
sPlot.fitted.RF <- predict(object = RF_FagusEU, data = Fagus.sPlot, type = "response")$predictions
sPlot.tss.RF <- ecospat.max.tss(Pred = sPlot.fitted.RF, Sp.occ = Fagus.sPlot$PA)$max.TSS
sPlot.boyce_ind.RF <- ecospat.boyce(fit = sPlot.fitted.RF, obs = sPlot.fitted.RF[which(Fagus.sPlot$PA == 1)], nclass = 0,
                                  window.w = "default", res = 100, PEplot = F)$Spearman.cor #check later details args

#goodness-of-fit
RF_FagusEU$r.squared #0.66

##Maps of predicted prob. of occurrence ----------------------------------------

#create a raster layer of latitude values
Lat_grid <- terra::init(BioClimateData.Eu$bio1, "y")
#mask the latitude layer as done for the bioclimatic variables
Lat_grid <- terra::mask(Lat_grid, mask = vect(IFS_poly))
#name the layer using the same name of corresponding predictor in the GLM and RF
names(Lat_grid) <- "lat"
#add the layer to the stack
Fagus_vars.stack <- c(BioClimateData.Eu, Lat_grid)

##Compare prediction maps RF and GLM ------------------------------------------------------
RF_FagusEU.pred <- predict(stack(Fagus_vars.stack), RF_FagusEU, type = "response",
                           predict.all = FALSE, na.rm = T, progress = "text",
                           fun = function(model, ...) predict(model, ...)$predictions)

par(mfrow = c(1, 2))
plot(RF_FagusEU.pred)

#note that color scales differ between figures (i.e range of predictions is much shorter for GLM)
plot(raster::predict(stack(Fagus_vars.stack), model = Mod_FagusEU.3, type = "response"))

##Functions used ---------------------------------------------------------------------------
#function for cross validation - now implemented for both GLM & RF
cross_val <- function(df, pa_col, formula, mod_type = c("GLM", "RF"), folds = 5L, ...) {
  require(ecospat)
  require(ranger)
  stopifnot(exprs = {
    inherits(df, "data.frame")
    inherits(pa_col, "character") && (length(pa_col) == 1)
    inherits(formula, "formula")
    inherits(folds, "integer")
  })
  if(!pa_col %in% colnames(df)) stop("PA col not found in df")
  mod_type <- match.arg(mod_type)
  n_pr <- length(which(df[[pa_col]] == 1))
  n_abs <- nrow(df) - n_pr
  id_vec <- c(sample(folds, n_pr, replace = TRUE),
              sample(folds, n_abs, replace = TRUE))
  metric_df <- do.call(rbind, lapply(unique(id_vec), function(id) {
    train <- df[id_vec != id, ]
    test <- df[id_vec == id, ]
    message(paste("Model fitted on:", sum(test$PA == 1), "presences and", sum(test$PA == 0), "absences"))
    fit_v <- switch(mod_type, "GLM" = {
      mod.glm <- glm(formula = formula, data = train, family = binomial)
      predict(mod.glm, newdata = test, type = "response")
      }, 
      "RF" = {
        mod.rf <- ranger(formula = formula, data = train, ...) #check with Dani
        predict(object = mod.rf, data = test, type = "response")$predictions
        })
    tss_stat <- ecospat.max.tss(Pred = fit_v, Sp.occ = test$PA)$max.TSS
    boy_index <- ecospat.boyce(fit = fit_v, obs = fit_v[which(test$PA == 1)], nclass = 0,
                               window.w = "default", res = 100, PEplot = F)$Spearman.cor #check later details args
    return(data.frame(TSS = tss_stat, Boyce = boy_index))
  })
  )
  return(metric_df)
}
