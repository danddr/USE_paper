## Overall description
Script and examples used for the USE paper

## Example of USE application with *Fagus sylvatica*
This repo includes i) the R scripts (Analysis) for replicating the analyses described in the paper [Da Re et al., (2023)](https://doi.org/10.1111/2041-210X.14209), and ii) the R script (Example) for replicating the analyses of the distribution of Fagus sylvatica across west Europe (Italy, France and Spain). 

The R script of the analyses provides a step-by-step guide on how to use the [USE](https://github.com/danddr/USE) R package (Daniele Da Re, Enrico Tordoni and Manuele Bazzichetto) to uniformly sample the environmental space and collect pseudo-absences for species distribution modelling.

Data used: [WorldClim](https://www.worldclim.org/data/index.html) bioclimatic variables (through the function raster::getData), sPlotOpen ([Sabatini et al., 2021](https://doi.org/10.1111/geb.13346)), and EU-Forest ([Mauri et al., 2017](https://doi.org/10.1038/sdata.2016.123)). The sPlotOpen and EU-Forest datasets can be found [here](https://idata.idiv.de/ddm/Data/ShowData/3474?version=54) and [here](https://figshare.com/articles/dataset/Occurrences_location_shapefile/3497891?backTo=/collections/A_high-resolution_pan-European_tree_occurrence_dataset/3288407), respectively. Notice that these datasets may be updated in the future, so results presented in the R script may change accordingly. Please, mail manuele.bazzichetto@gmail.com to get the same data version used for the analysis of Fagus. Also, notice that these data are intended as a step-by-step description on how to use the UEsampling package. Finally, WorldClim, sPlotOpen and EU-Forest are open datasets, but please check their data usage license and follow the guidelines included in the related publications for citation.

The R script is roughly divided in three sections:
1) Data from the above-mentioned sources are downloaded (WorldClim) and cleaned for the analyses.
2) The uniform sampling of the environmental space is implemented (using the UEsampling package) to i) apply kind of a spatial thinning procedure for subsetting presence data of Fagus sylvatica, and ii) collecting background points within the environmental space.
3) The obtained training and testing datasets are used to calibrate and validate, respectively, species ditribution models for Fagus sylvatica fitted using binary generalised linear models and random forest.
