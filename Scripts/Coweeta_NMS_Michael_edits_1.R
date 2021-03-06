#Coweeta -- Nonmetric Multidimensional Scaling
#Enikoe Bihari, Michael Gaffney, Cal Oakley

#MAKE SURE TO RUN COWEETA_EDA.R BEFORE RUNNING THIS SCRIPT
library(vegan)

#check working directory; should be the project folder; used setwd if not.
getwd()

#Load in species data from the Coweeta_EDA script first; check to make sure these datasets are already here
View(coweeta.env)
View(coweeta.species)

#first, we need to cull the environmental variables that are correlated
#we will use a .7 cutoff for this; ie variables above a .7 correlation will be tossed (we are only doing this for visual purposes)
#load this convenient script from Urban to return a list of variables with high correlations
source("./scripts/screen_cor.R")
screen.cor(coweeta.env[,-1])

#remove some of the correlated variables
coweeta.env.11 <- coweeta.env %>%
  select(-c("N", "ECEC", "BS", "Silt", "Acidity", "Mg"))

#next, we need to relativize the species data to account for the fact that some species are far more common than others
#source the helper function that comes directly from an R script, not a package (thanks again to Urban)
source("./scripts/reldata.R")
#now we can relativize the data using this function from the script
#relativization method: column maxima and row sums, sometimes called "Wisconsin double"
coweeta.species.rel <- reldata(coweeta.species[,-1], byrow=T, bycol=T, rowfirst=F)
#new dataset as relative values of all species per plot, similar to the original set

#now convert the realtivized values to a distance matrix using the Bray-Curtis method
#load ecodist library
library(ecodist)
coweeta.species.bcd <- distance(as.matrix(coweeta.species.rel), method="bray")

#now that's we've built a matrix of ecological distances, we need to check whether the distance measure is saturated
#distance matrix size
length(coweeta.species.bcd)
#number of saturated plots?
length(coweeta.species.bcd[coweeta.species.bcd==1])
#saturation ratio
length(coweeta.species.bcd[coweeta.species.bcd==1])/length(coweeta.species.bcd)
#the distance matrix is not saturated, so we don't need to do anything else.

#now that we have a distance matrix that will actually work, let's run the NMS.
#NB: this is a first run, and it will create 60 ordinations (stepdown, 1:6 dims, 10 iterations each). 
#We will use these to figure out what the best number of axes is.

coweeta.species.nmsstep <- nmds(coweeta.species.bcd, nits=10, mindim=1, maxdim=6)
attributes(coweeta.species.nmsstep)

#Now we need to plot the stress values and R2 values for each dimension (# of axes)
plot.nmds(coweeta.species.nmsstep) 
#second plot shows the R^2 explained by # of dimensions; 2 axes gets us about 60% of the variance, so we can try that.
#run the NMS again for 2 dimensions, with 20 repetitions
coweeta.species.nmds <- nmds(coweeta.species.bcd, mindim=2, maxdim=2, nits=20)
#which of these 20 repetitions did the best? NB: NMS is stochastic, so every run is a little different
(s.min <- which.min(coweeta.species.nmds$stress))
#this returns the "stress" value for the best one--need to check what that value actually means
coweeta.species.nmds$stress[s.min]
#grab the best run; this explains about 60% of the variance, which is solid but not fantastic
coweeta.species.nms <- nmds.min(coweeta.species.nmds,dims=2)

#Now we need to rotate the NMS so that we can actually understand what the hell it did. 
#Use a principle component analysis to force most of the variance into axis 1, leaving the rest for axis 2
nms.pca <- princomp(coweeta.species.nms)
print(nms.pca)
summary(nms.pca)
#now replace the original NMS scores with the rotated PCA scores
coweeta.species.nms <- nms.pca$scores
#rename the columns for purposes of clarity
colnames(coweeta.species.nms) <- c("NMS1", "NMS2")

#the whole point of NMS is to get the sample points that are ecologically similiar to be close together on a graph
#A Shepard Diagram let's us know how well we did this.
#compute the distances between points in the NMS
nms2.xod <- dist(coweeta.species.nms)
#plot the NMS distances against the original Extended Bray Curtis Distances
plot(nms2.xod, coweeta.species.bcd, pch="*",xlab="Ordination Distance", ylab="Extended B-C Distance")
abline(0,1,col="red") # put in the 1:1 line (intercept=0, slope=1)
title("Shepard Diagram")
box(lwd=2)
#This came out alright, but not fantastic (we'd want the points to hug the line better). Check with DLU to see if work is required on this.
#We can also run a Mantel test to see what the correlations actually are; but this code isn't playing nice right now
#ecodist::mantel(nms2.xod ~ coweeta.species.bcd,nperm=10000,nboot=0)

#now we want to get the R^2 for the NMS
#This is done by differencing for NMS, since it doesn't work like a linear model at all
nms.od1 <- dist(coweeta.species.nms[,1])
nms.od2 <- dist(coweeta.species.nms[,1:2])
#check axis 1
r1 <- cor(coweeta.species.bcd,nms.od1)
r2.1 <- r1^2; r2.1
# axis 2 is 2-D minus 1-D solution:
r2 <- cor(coweeta.species.bcd,nms.od2)
r2.2 <- r2^2; r2.2-r2.1; r2.2
#first number is the amount captured by the second axis; second is total variance

#The final stage for NMS is to actually interpret what we have in a plot.
#first, corrleate the NMS axes (which come from species data) with the environmental variables
coweeta.species.nms.cor2m <- as.matrix(cor2m(coweeta.env.11[,-1], coweeta.species.nms))
#this is a table we'll probably want at some point, so make a note of it

#now we need to get the species into NMS space; we do this by creating weighted averages of them.
coweeta.species.nms.wa <- vegan::wascores(coweeta.species.nms, coweeta.species[,-1])

#now let's plot the NMS axes. I'm using base R here but we probably want to transfer this over to ggplot for aestheic purposes.
#First, this gives us the base plot; NMS axes for x and y, and the 108 plots in NMS space
plot(coweeta.species.nms[,2:1],pch=19, xlab="NMS 2",ylab="NMS 1")
#now, let's add the species names on top of this.
text(coweeta.species.nms.wa[,2:1], rownames(coweeta.species.nms.wa),cex=0.8, col = "blue")
#this shows us roughly where the species themselves, apart from the plots, show up in ordination space
#now, let's add in correlation vectors for the environmental variables; this is how we know what the axes mean
#use the vector fitting function from the ecodist package for this
coweeta.species.nms.vf <- vf(coweeta.species.nms[,2:1], coweeta.env.11[,-1],nperm=1000)
#now plot the vectors; note that we're only plotting the highly statistically significant vectors
plot.vf(coweeta.species.nms.vf, pval=0.01, col="red", lwd=2, length=0.067)
# (length specifies the size of the arrowheads, in inches, and it annoys R)
box(lwd=2)
#NMS 1 has something to do with soil content: Carbon Nitrogen ratio is positively correlated, potassium negatively
#NMS 2 seems to be a very weak elevation gradient; given that elevation isn't super important here, that's not a surprise

#We can plot this again to better iullstrate how the species are sorting in NMS space. This gets messsy visually, so two plots is better
#First, create a plot that illustrates the soil moisture gradient
#plot NMS space
plot(coweeta.species.nms[,2:1],pch=3, xlab="NMS 2",ylab="NMS 1")
#color species by relative abundance!
#red maple is everywhere; I use it as a baselayer
points(coweeta.species.nms[,2:1],pch=19,col="gray", cex=6.0*coweeta.species.rel$ACRU)
#DLU says that scarlet oak is following the sopil moisture gradient in this space, so let's plot that
points(coweeta.species.nms[,2:1],pch=19,col="purple", cex=6.0*coweeta.species.rel$QUCO)
#DLU also said mountain laurels are a pretty good guess here, so let's do that
points(coweeta.species.nms[,2:1],pch=19,col="red", cex=6.0*coweeta.species.rel$KALA)
#pitch pine = low soil moisture
points(coweeta.species.nms[,2:1],pch=19,col="blue", cex=6.0*coweeta.species.rel$PIRI)
#DLU also says that birches would sort on the low end of the soil moisutre gradient, so add those in as well
points(coweeta.species.nms[,2:1],pch=19,col="green", cex=6.0*coweeta.species.rel$BELE)
#tulip seems to sort on high moisture
points(coweeta.species.nms[,2:1],pch=19,col="orange", cex=6.0*coweeta.species.rel$LITU)
#add in the legend names and colors:
legend.txt <- c("ACRU","QUCO","KALA","PIRO","BELE", "LITU")
legend.col <- c("gray","purple", "red","blue","green","orange")
legend("topright",horiz=F,legend=legend.txt,text.col=legend.col,bty="n")
box(lwd=2)

#let's make a second plot to illustrate the elevation gradient on NMS 2, same method
plot(coweeta.species.nms[,2:1],pch=3, xlab="NMS 2",ylab="NMS 1")
#red maple baselayer
points(coweeta.species.nms[,2:1],pch=19,col="gray", cex=6.0*coweeta.species.rel$ACRU)
#white oak seems to do very in low elevation
points(coweeta.species.nms[,2:1],pch=19,col="lightgreen", cex=6.0*coweeta.species.rel$QUAL)
#dogwoods do very well in low elevation too
points(coweeta.species.nms[,2:1],pch=19,col="darkgreen", cex=6.0*coweeta.species.rel$COFL)
#great laurels clump in the heigher elevations.
points(coweeta.species.nms[,2:1],pch=19,col="darkred", cex=6.0*coweeta.species.rel$RHMA)
#Black locusts also work well in this case.
points(coweeta.species.nms[,2:1],pch=19,col="red", cex=6.0*coweeta.species.rel$ROPS)

#add in the legend names and colors:
legend.txt <- c("ACRU","QUAL","COFL","RHMA","ROPS")
legend.col <- c("gray","lightgreen", "darkgreen","darkred","red")
legend("topright",horiz=F,legend=legend.txt,text.col=legend.col,bty="n")
box(lwd=2)

#that's all folks