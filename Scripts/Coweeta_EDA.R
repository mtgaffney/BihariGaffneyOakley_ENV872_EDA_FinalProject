#Coweeta data Initial EDA
#Enikoe Bihari, Michael Gaffney, Cal Oakley

#check working directory; should be the project folder; used setwd if not.
getwd()

#load packages
library(tidyverse)
library(readxl)
library(sf)
library(mapview)
library(cowplot)

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

# #epxlore the species data a little; try a random plot of some of the species (Acer Rubrum, Rhododendron maximum, Carya tomentosa) against elevation
# #set ggplot data source to null because we're using two different datasets
# plot1 <- ggplot(NULL, aes(x = coweeta.env$Elevation)) +
#   geom_point(aes(y = coweeta.species$ACRU), color = "darkgreen") +
#   geom_point(aes(y = coweeta.species$RHMA), color = "darkred") +
#   geom_point(aes(y = coweeta.species$CATO), color = "darkblue") +
#   labs(x = "Elevation", y = "Basal Area")
# plot1

################################################################################
################################################################################
################################################################################

# ggplot stuff (from Eni + Cal)

######################## BA for all species ######################## 

#display the distribution of basal areas by species

#first pivot_longer & subset to canopy species
coweeta.species.canopy <- coweeta.species %>% 
  select(Plot:ACRU, BELE:CATO, LITU:NYSY, PIRI:QUVE, ROPS) %>% 
  pivot_longer(!Plot, names_to = 'speciesCode', values_to = 'BasalArea')

#If you want to do this for understory species as well
coweeta.species.understory <- coweeta.species %>%
  select(Plot,AMAR, COFL:KALA, OXAR, RHMA) %>%
  pivot_longer(!Plot, names_to = 'speciesCode', values_to = 'BasalArea')

#plot the mean basal area for each species across all plots to get sense of 
#landscape level distributions
# plot.BA.sum <- ggplot(coweeta.species.canopy,
#                   aes(x =speciesCode, y=BasalArea, fill = speciesCode)) +
#   stat_summary(fun = sum, geom = 'bar') +
#   labs(y = expression(paste('Basal Area (', "ft"^2, ')')), 
#        x = 'Species Code',
#        title = 'Mean basal area for canopy species at Coweeta LTERS',
#        fill = "Species Code") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# plot.BA.sum
# 
# coweeta.species.canopy.means = 
#   coweeta.species.canopy %>% 
#   group_by(speciesCode) %>% 
#   summarize(meanBA = mean(BasalArea))

plot.BA.mean.canopy <- 
  ggplot(coweeta.species.canopy,
         aes(x =speciesCode, 
             y=BasalArea)) +
  stat_summary(fun = mean, 
               geom = 'bar',
               alpha = 0.6,
               fill = 'seagreen1',
               color = 'seagreen',
               size = 1) +
  labs(y = expression(paste('Basal Area (', "ft"^2, ')')), 
       x = 'Species Code',
       title = 'Mean basal area per plot for canopy species at Coweeta LTERS',
       fill = "Species Code") +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
plot.BA.mean.canopy

plot.BA.mean.ustory <- 
  ggplot(coweeta.species.understory,
         aes(x =speciesCode,
             y=BasalArea)) +
  stat_summary(fun = mean, 
               geom = 'bar', 
               alpha = 0.6,
               fill = 'lightsalmon',
               color = 'lightsalmon3',
               size = 1) +
  labs(y = expression(paste('Basal Area (', "ft"^2, ')')), 
       x = 'Species Code',
       title = 'Mean basal area per plot for understory species at Coweeta LTERS',
       fill = "Species Code") +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
plot.BA.mean.ustory

######################## BA for most common species ######################## 

#plot one species' basal area for each plot 

# ACRU 

#first subset just for desired species
coweeta.ACRU <- coweeta.species.canopy %>% 
  filter(speciesCode == "ACRU")

mean.ACRU = as.character(round(mean(coweeta.ACRU$BasalArea), digits = 3))
mean.ACRU

#then plot
plot.ACRU <-  
  ggplot(coweeta.ACRU, aes(x=BasalArea)) + 
  geom_area(aes(y=..count..), 
            stat = 'bin',
            binwidth = 3,
            fill = '#C0FF3E',
            alpha = .6) +
  geom_vline(xintercept = mean(coweeta.ACRU$BasalArea), 
             color = '#6B8E23',
             size = 1) +
  theme(legend.position = 'none') +
  annotate(geom = 'text', 
           x = 12, 
           y = 30, 
           color = '#6B8E23',
           label = expression(paste("mean = ", 4.442, " ft"^2))) +
  xlim(0,45) +
  ylim(0,75) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Red maple (ACRU)")
plot.ACRU

#QUPR

coweeta.QUPR <- coweeta.species.canopy %>% 
  filter(speciesCode == "QUPR")

mean.QUPR = as.character(round(mean(coweeta.QUPR$BasalArea), digits = 3))
mean.QUPR

plot.QUPR <-  
  ggplot(coweeta.QUPR, aes(x=BasalArea)) + 
  geom_area(aes(y=..count..), 
            stat = 'bin',
            binwidth = 3,
            fill = 'lightsalmon',
            alpha = .6) +
  geom_vline(xintercept = mean(coweeta.QUPR$BasalArea), 
             color = "lightsalmon3",
             size = 1) +
  theme(legend.position = 'none') +
  annotate(geom = 'text', 
           x = 15, 
           y = 25, 
           color = "lightsalmon3",
           label = expression(paste("mean = ", 7.087, " ft"^2))) +
  xlim(0,45) +
  ylim(0,75) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Chestnut oak (QUPR)")
plot.QUPR

#QURU

coweeta.QURU <- coweeta.species.canopy %>% 
  filter(speciesCode == "QURU")

mean.QURU = as.character(round(mean(coweeta.QURU$BasalArea), digits = 3))
mean.QURU

plot.QURU <-  
  ggplot(coweeta.QURU, aes(x=BasalArea)) + 
  geom_area(aes(y=..count..), 
            stat = 'bin',
            binwidth = 3,
            fill = '#FBDB0C',
            alpha = .4) +
  geom_vline(xintercept = mean(coweeta.QURU$BasalArea), 
             color = "#CDAD00",
             size = 1) +
  theme(legend.position = 'none') +
  annotate(geom = 'text', 
           x = 11, 
           y = 50, 
           color = "#CDAD00",
           label = expression(paste("mean = ", 3.276, " ft"^2))) +
  xlim(0,45) +
  ylim(0,75) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Northern red oak (QURU)")
plot.QURU

#QUCO

coweeta.QUCO <- coweeta.species.canopy %>% 
  filter(speciesCode == "QUCO")

mean.QUCO = as.character(round(mean(coweeta.QUCO$BasalArea), digits = 3))
mean.QUCO

plot.QUCO <-  
  ggplot(coweeta.QUCO, aes(x=BasalArea)) + 
  geom_area(aes(y=..count..), 
            stat = 'bin',
            binwidth = 3,
            fill = 'seagreen1',
            alpha = .4) +
  geom_vline(xintercept = mean(coweeta.QUCO$BasalArea), 
             color = "seagreen",
             size = 1) +
  theme(legend.position = 'none') +
  annotate(geom = 'text', 
           x = 10, 
           y = 40, 
           color = "seagreen",
           label = expression(paste("mean = ", 2.564, " ft"^2))) +
  xlim(0,45) +
  ylim(0,75) +
  labs(y="Plot Count", x=expression(paste('Basal Area (', "ft"^2, ')')),
       title = "Scarlet Oak (QUCO)")
plot.QUCO

# cowplot of all

# make title
title <- ggdraw() + 
  draw_label(
    "Frequency distribution of areas for dominant canopy trees in plots at Coweeta LTRS",
    x = 0,
    hjust = 0,
    size = 14,) 

  # theme(
  #   # add margin on the left of the drawing canvas,
  #   # so title is aligned with left edge of first plot
  #   plot.margin = margin(0, 0, 0, 7)
  # )

plots_grid <-
  plot_grid(plot.QUPR,
            plot.QURU,
            plot.ACRU,
            plot.QUCO,
            nrow = 2, 
            align = 'hv')

plots_grid.title = 
  plot_grid(title, 
            plots_grid, 
            nrow = 2, 
            align = 'hv',
            rel_heights = c(0.1, 2))
print(plots_grid.title)

# all on one plot - DID NOT LOOK VERY GOOD

# coweeta.dom <- coweeta.species.canopy %>% 
#   filter(speciesCode %in% c("QUCO", "ACRU", "QURU", "QUPR"))
# 
# plot.BA.histograms.dom <-  
#   
#   ggplot(coweeta.dom, 
#          aes(x=BasalArea)) + 
#   
#   # # all at once
#   # geom_area(aes(y=..count.., fill = speciesCode, group = speciesCode),
#   #           stat = 'bin',
#   #           binwidth = 3,
#   #           alpha = 0.5) +
#   
#   # ACRU
#   geom_area(data = coweeta.ACRU,
#             aes(y=..count..),
#             stat = 'bin',
#             binwidth = 3,
#             alpha = 0.5,
#             color = "yellow",
#             fill = NA) +
#   geom_vline(xintercept = mean(coweeta.ACRU$BasalArea), color = "yellow") +
#   # annotate(geom = 'text',
#   #          x = 6,
#   #          y = 30,
#   #          color = "yellow",
#   #          label = paste('mean = ', round(mean(coweeta.QUCO$BasalArea),digits = 3),expression(ft^2))) +
#   
#   # QURU
#   geom_area(data = coweeta.QURU,
#             aes(y=..count..),
#             stat = 'bin',
#             binwidth = 3,
#             alpha = 0.5,
#             color = "blue",
#             fill = NA) +
#   geom_vline(xintercept = mean(coweeta.QURU$BasalArea), color = "blue") +
#   # annotate(geom = 'text',
#   #          x = 6,
#   #          y = 30,
#   #          color = "blue",
#   #          label = paste('mean = ', round(mean(coweeta.QUCO$BasalArea),digits = 3),expression(ft^2))) +
#   
#   # QUCO
#   geom_area(data = coweeta.QUCO,
#             aes(y=..count..),
#             stat = 'bin',
#             binwidth = 3,
#             alpha = 0.5,
#             color = "red",
#             fill = NA) +
#   geom_vline(xintercept = mean(coweeta.QUCO$BasalArea), color = "red") +
#   # annotate(geom = 'text',
#   #          x = 6,
#   #          y = 30,
#   #          color = "red",
#   #          label = paste('mean = ', round(mean(coweeta.QUCO$BasalArea),digits = 3),expression(ft^2))) +
# 
#   # QUPR
#   geom_area(data = coweeta.QUPR,
#             aes(y=..count..),
#             stat = 'bin',
#             binwidth = 3,
#             alpha = 0.5,
#             color = "green",
#             fill = NA) +
#   geom_vline(xintercept = mean(coweeta.QUPR$BasalArea), color = "green") +
#   # annotate(geom = 'text',
#   #          x = 6,
#   #          y = 30,
#   #          color = "green",
#   #          label = paste('mean = ', round(mean(coweeta.QUCO$BasalArea),digits = 3),expression(ft^2))) +
# 
#   xlim(0,45) +
#   theme(legend.position = 'left') +
#   labs(y="Plot Count", 
#        x=expression(paste('Basal Area (', "ft"^2, ')')),
#        title = "Frequency distribution of basal areas for most common canopy tree species \nin plots at Coweeta LTRS")
# 
# plot.BA.histograms.dom

######################## Correlation stuff ######################## 

#check out correlations in the environmental variables; there will be a bunch of these
names(coweeta.env)
coweeta.env.cor <- data.frame(cor(coweeta.env[,-1]))
print(coweeta.env.cor,digits=3)
#save correlations for presentation
write.csv(coweeta.env.cor, './Data/Outputs/CoweetaVarCors.csv', row.names = TRUE)

#check a correaltion directly betwene two variables.
cor.test(coweeta.env$Acidity, coweeta.env$C)
