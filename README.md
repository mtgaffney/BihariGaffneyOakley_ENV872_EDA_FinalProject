# BihariGaffneyOakley_ENV872_EDA_FinalProject README


## Summary
This project folder contains data and scripts for a final project for a Spring 2022 Nicholas School of the Environment class--Environmental Data Analytics--using Coweeta LTER data. The central question of this project was to determine what environmental factors were sorting tree communities in Coweeta. In order to answer this question, we conducted Nonmetric Multidimensional Scaling (NMS) on the tree species data, and then overlayed environmental variables on the ordination. Our findings suggest that tree species in Coweeta tend to sort on soil moisture and elevation, primarily. 

## Investigators

Eni Bihari (enikoe.bihari@duke.edu), Michael Gaffney (michael.gaffney@duke.edu), and Cal Oakley (daniel.oakley@duke.edu)

## Keywords

LTER; Southern Appalachia; Coweeta; Nonmetric Multidimensional Scaling

## Database Information

Three datasets in total were used for this analysis, all provided by Dean Urban and acquired between 1997 and 1999. Data collection was organized through 108 sample quadrats (20x20m) arranged in clusters of 3 to 4 along transects that covered an elevation gradient. Species data consists of the basal area of 18 species across 108 plots in Coweeta; relatively uncommon species were removed prior to analysis. Environmental data consists of 20 environmental variables across the same 108 sites. Finally, we were provided geospatial data for the sites.


## Folder structure, file formats, and naming conventions 

Data = contains all excel sheet data used in this project
Docs = contains output PDF from knitting
Scripts = contains all scripts and code written for this analysis

.R script files were used for drafting code and for helper function. .RMD was used for the final report; .xslx was the original format of the data.

Main script files were named with coweeta at the beginning and then the section of the project for which the script was written.

## Metadata

1. chl_sppdataBA contains basal area of these species per each 108 plot:

Species Code |Scientific Name           |Common Name  
-------------|--------------------------|--------------------------
ACRU	       |Acer rubrum	              |Red Maple
AMAR	       |Amelanchier arborea	      |Common Serviceberry
BELE	       |Betula lenta	            |Sweet Birch
CAGL	       |Carya glabra	            |Pignut Hickory
CATO	       |Carya tomentosa	          |Mockernut Hickory
COFL	       |Cornus florida	          |Flowering Dogwood
KALA	       |Kalmia latifolia	        |Mountain Laurel
LITU	       |Liriodendron tulipifera	  |Tulip Tree
NYSY	       |Nyssa sylvatica	          |Black Tupelo
OXAR	       |Oxydendron arboreum	      |Sourwood
PIRI	       |Pinus rigida	            |Pitch Pine
QUAL	       |Quercus alba	            |White Oak
QUCO	       |Quercus coccinea	        |Scarlet Oak
QUPR	       |Quercus prinus	          |Chestnut Oak
QURU	       |Quercus rubra	            |Red Oak
QUVE	       |Quercus velutina	        |Black Oak
RHMA	       |Rhododendron maximum	    |Great Laurel
ROPS	       |Robinia pseudoacacia	    |Black Locust


2. chl_sitedata contains environmental data per each 108 plot.

Variable       |Description
---------------|---------------------------------
Elevation	     |Elevation (M)
SlopeDEM	     |Slope
TAspectD       |Transformed Aspect
TCI	           |Topograhpic Convergence Index
xDepth	       |Mean soil depth (cm)
sDepth	       |Standard deviation of depth
pH	           |pH
Acidity	       |Acidity
Ca	           |Calcium (cmol(+)/kg)
K	             |Potassium (cmol(+)/kg)
Mg	           |Magnesium (cmol(+)/kg)
P	             |Phosphorous (g/g)
C	             |Soil carbon (%)
N	             |Nitrogen (%)
C:N	           |Carbon:Nitrogen ratio
ECEC	         |Effective Cation Exchange Capacity
BS	           |Base Saturation
Clay	         |Clay (%)
Silt	         |Silt (%)
Sand	         |Sand (%)

3. chl_XYdata contains XY coordinates for each 108 plot.

## Scripts and code

Coweeta_EDA.R contains the initial data exploration. Other Coweeta_EDA scripts contain draft code.

Coweeta_NMS.R contains the NMS analysis. Other Coweeta_NMS scripts contain draft code.

EDA_ProjectFinal_Michael_Cal_Eni.rmd contaisn the combined analysis and EDA and the is the final script.

screen_cor.R is a script provided by Dean Urban that correlates variables in a datatset and returns only those correlations above a certain threshold.

reldata.R is a script provided by Dean Urban that relativizes data, which is useful in cases where the basal area is substantially higher for some species.

## Quality assurance/quality control

QA procedures taken by Dean Urban before providing the dataset.