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
coweeta.species <- as.data.frame(read_excel("./Data/Raw/chl_sppdataBA.xlsx", sheet = "dataBA"))
#load the environmental variables for the same 108 plots
coweeta.env <- as.data.frame(read_excel("./Data/Raw/chl_sitedata.xlsx"))
#load the geospatial data for the sites (sheet CHL has the sites)
coweeta.xy <- as.data.frame(read_excel("./Data/Raw/chl_XYdata.xlsx", sheet = "CHL"))

#sweet, let's get the XY data into a spatial format; UTM 17N is the projection
coweeta.sf <- coweeta.xy %>% 
  st_as_sf(coords = c('UTME','UTMN'),
           crs=32617)
#check out where this stuff shows up on the map
mapView(coweeta.sf, col.regions = "red", map.types = "Esri.WorldImagery", legend = FALSE)

#epxlore the species data a little; try a random plot of some of the species (Acer Rubrum, Rhododendron maximum, Carya tomentosa) against elevation
#set ggplot data source to null because we're using two different datasets
plot1 <- ggplot(NULL, aes(x = coweeta.env$Elevation)) +
  geom_point(aes(y = coweeta.species$ACRU), color = "darkgreen") +
  geom_point(aes(y = coweeta.species$RHMA), color = "darkred") +
  geom_point(aes(y = coweeta.species$CATO), color = "darkblue") +
  labs(x = "Elevation", y = "Basal Area")
plot1

#display the distribution of basal areas by species
#first pivot_longer & subset to canopy species
coweeta.species.canopy <- coweeta.species %>% 
  select(Plot:ACRU, BELE:CATO, LITU:NYSY, PIRI:QUVE, ROPS) %>% 
  pivot_longer(!Plot, names_to = 'speciesCode', values_to = 'BasalArea')

#If you want to do this for understory species as well........................
# coweeta.species.understory <- coweeta.species %>%
#   select(Plot,AMAR, COFL:KALA, OXAR, RHMA) %>%
#   pivot_longer(!Plot, names_to = 'speciesCode', values_to = 'BasalArea')

#plot the mean basal area for each species across all plots to get sense of 
#landscape level distributions
plot.BA.sum <- ggplot(coweeta.species.canopy,
                  aes(x =speciesCode, y=BasalArea, fill = speciesCode)) +
  stat_summary(fun = sum, geom = 'bar') +
  labs(y = expression(paste('Basal Area (', "ft"^2, ')')), 
       x = 'Species Code',
       title = 'Mean basal area for canopy species at Coweeta LTERS',
       fill = "Species Code") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot.BA.sum

plot.BA.mean <- ggplot(coweeta.species.canopy,
                  aes(x =speciesCode, y=BasalArea, fill = speciesCode)) +
  stat_summary(fun = mean, geom = 'bar') +
  labs(y = expression(paste('Basal Area (', "ft"^2, ')')), 
       x = 'Species Code',
       title = 'Mean basal area for canopy species at Coweeta LTERS',
       fill = "Species Code") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot.BA.mean

#plot one species' (ACRU) basal area for each plot 
#first subset just for desired species
coweeta.ACRU <- coweeta.species.canopy %>% 
  filter(speciesCode == "ACRU")

#then plot
plot.ACRU <-  ggplot(coweeta.ACRU, aes(x=BasalArea)) + 
  geom_area(aes(y=..count.., fill = speciesCode, group = speciesCode), 
            stat = 'bin') +
  geom_vline(xintercept = mean(coweeta.ACRU$BasalArea), color = "black") +
  theme(legend.position = 'none') +
  annotate(geom = 'text', x = 7, y = 10, label = paste('mean = ', 
    round(
      mean(coweeta.ACRU$BasalArea), 
      digits = 3), 
    expression(ft^2))) +
  # xlim(0,20) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Frequency distribution of basal areas for plots at Coweeta LTRS")
plot.ACRU

#QUPR
coweeta.QUPR <- coweeta.species.canopy %>% 
  filter(speciesCode == "QUPR")

plot.QUPR <-  ggplot(coweeta.QUPR, aes(x=BasalArea)) + 
  geom_area(aes(y=..count.., fill = speciesCode, group = speciesCode), 
            stat = 'bin') +
  geom_vline(xintercept = mean(coweeta.QUPR$BasalArea), color = "black") +
  theme(legend.position = 'none') +
  annotate(geom = 'text', x = 4.5, y = 15, label = paste('mean = ', 
                                                       round(
                                                         mean(coweeta.QUPR$BasalArea), 
                                                         digits = 3), 
                                                       expression(ft^2))) +
  # xlim(0,20) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Frequency distribution of basal areas for plots at Coweeta LTRS")
plot.QUPR

#QURU
coweeta.QURU <- coweeta.species.canopy %>% 
  filter(speciesCode == "QURU")

plot.QURU <-  ggplot(coweeta.QURU, aes(x=BasalArea)) + 
  geom_area(aes(y=..count.., fill = speciesCode, group = speciesCode), 
            stat = 'bin') +
  geom_vline(xintercept = mean(coweeta.QURU$BasalArea), color = "black") +
  theme(legend.position = 'none') +
  annotate(geom = 'text', x = 6, y = 30, label = paste('mean = ', 
                                                         round(
                                                           mean(coweeta.QURU$BasalArea), 
                                                           digits = 3), 
                                                         expression(ft^2))) +
  # xlim(0,20) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Frequency distribution of basal areas for plots at Coweeta LTRS")
plot.QURU

#check out correlations in the environmental variables; there will be a bunch of these
names(coweeta.env)
coweeta.env.cor <- data.frame(cor(coweeta.env[,-1]))
print(coweeta.env.cor,digits=3)
#save correlations for presentation
write.csv(coweeta.env.cor, './Data/Outputs/CoweetaVarCors.csv', row.names = TRUE)

#check a correaltion directly betwene two variables.
cor.test(coweeta.env$Acidity, coweeta.env$C)


