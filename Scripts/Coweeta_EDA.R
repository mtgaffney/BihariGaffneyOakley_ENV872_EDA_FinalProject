#Coweeta data Initial EDA
#Enikoe Bihari, Michael Gaffney, Cal Oakley

#check working directory; should be the project folder; used setwd if not.
getwd()

#load packages
library(tidyverse)
library(readxl)
library(sf)
library(mapview)

#load the species data for 108 coweeta plots (rare species already thinned)
coweeta.species <- read_excel("./Data/Raw/chl_sppdataBA.xlsx", sheet = "dataBA")
#load the environmental variables for the same 108 plots
coweeta.env <- read_excel("./Data/Raw/chl_sitedata.xlsx")
#load the geospatial data for the sites (sheet CHL has the sites)
coweeta.xy <- read_excel("./Data/Raw/chl_XYdata.xlsx", sheet = "CHL")

#sweet, let's get the XY data into a spatial format; UTM 17N is the projection
coweeta.sf <- coweeta.xy %>% 
  st_as_sf(coords = c('UTME','UTMN'),
           crs=32617)
#check out where this stuff shows up on the map
mapView(coweeta.sf, col.regions = "darkgreen", map.types = "CartoDB.Positron", legend = FALSE)

#epxlore the species data a little; try a random plot of some of the species (Acer Rubrum, Rhododendron maximum, Carya tomentosa) against elevation
#set ggplot data source to null because we're using two different datasets
ggplot(NULL, aes(x = coweeta.env$Elevation)) +
  geom_point(aes(y = coweeta.species$ACRU), color = "darkgreen") +
  geom_point(aes(y = coweeta.species$RHMA), color = "darkred") +
  geom_point(aes(y = coweeta.species$CATO), color = "darkblue") +
  labs(x = "Elevation", y = "Basal Area")
  
