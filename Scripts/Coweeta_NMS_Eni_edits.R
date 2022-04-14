install.packages("installr")

library(installr)

updateR()#Coweeta -- Nonmetric Multidimensional Scaling
#Enikoe Bihari, Michael Gaffney, Cal Oakley

#MAKE SURE TO RUN COWEETA_EDA.R BEFORE RUNNING THIS SCRIPT
install.packages("ggvegan")
# install.packages("ggpubr")
library(vegan)
library(ggplot2)
library(ggvegan)
library(ggpubr)

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
#becuase the distance matrix is saturdated, we're going to use extended bray curtis here. Might not be necessary.
coweeta.species.xbcd <- stepacross(coweeta.species.bcd)

#now that we have a distance matrix that will actually work, let's run the NMS.
#NB: this is a first run, and it will create 60 ordinations (stepdown, 1:6 dims, 10 iterations each). 
#We will use these to figure out what the best number of axes is.

coweeta.species.nmsstep <- nmds(coweeta.species.xbcd, nits=10, mindim=1, maxdim=6)
attributes(coweeta.species.nmsstep)

#Now we need to plot the stress values and R2 values for each dimension (# of axes)
plot.nmds(coweeta.species.nmsstep) 
#second plot shows the R^2 explained by # of dimensions; 2 gets us about 60% of the variance, so we can try that.
#run the NMS again for 2 dimensions, with 20 repetitions
coweeta.species.nmds <- nmds(coweeta.species.xbcd, mindim=2, maxdim=2, nits=20)
#which of these 20 repetitions did the best? NB: NMS is stochastic, so every run is a little different
(s.min <- which.min(coweeta.species.nmds$stress))
#this returns the "stress" value for the best one--need to check what that value actually means
coweeta.species.nmds$stress[s.min]
#grab the best run; this explains about 60% of the variance, which is solid but not fantastic
coweeta.species.nms <- nmds.min(coweeta.species.nmds,dims=2)

fort = fortify(coweeta.species.nms)

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
plot(nms2.xod, coweeta.species.xbcd, pch="*",xlab="Ordination Distance", ylab="Extended B-C Distance")
abline(0,1,col="red") # put in the 1:1 line (intercept=0, slope=1)
title("Shepard Diagram")
box(lwd=2)
#This came out alright, but not fantastic (we'd want the points to hug the line better). Check with DLU to see if work is required on this.
#We can also run a Mantel test to see what the correlations actually are; but this code isn't playing nice right now
#ecodist::mantel(nms2.xod ~ coweeta.species.xbcd,nperm=10000,nboot=0)

#now we want to get the R^2 for the NMS
#This is done by differencing for NMS, since it doesn't work like a linear model at all
nms.od1 <- dist(coweeta.species.nms[,1])
nms.od2 <- dist(coweeta.species.nms[,1:2])
#check axis 1
r1 <- cor(coweeta.species.xbcd,nms.od1)
r2.1 <- r1^2; r2.1
# axis 2 is 2-D minus 1-D solution:
r2 <- cor(coweeta.species.xbcd,nms.od2)
r2.2 <- r2^2; r2.2-r2.1; r2.2
#first number is the amount captured by the second axis; second is total variance

#The final stage for NMS is to actually interpret what we have in a plot.
#first, corrleated the NMS axes (which come from species data) with the environmental variables
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
coweeta.species.nms.vf <- vf(coweeta.species.nms[,2:1], coweeta.env.11[,-1], nperm=1000)
print(coweeta.species.nms.vf)
#now plot the vectors; note that we're only plotting the highly statistically significant vectors
plot.vf(coweeta.species.nms.vf, pval=0.01, col="red", lwd=2, length=0.067)
# (length specifies the size of the arrowheads, in inches, and it annoys R)
box(lwd=2)
#NMS 1 has something to do with soil content: Carbon Nitrogen ratio is positively correlated, potassium negatively
#NMS 2 seems to be a very weak elevation gradient; given that elevation isn't super important here, that's not a surprise

#We can plot this again to better iullstrate how the species are sorting in NMS space
plot(coweeta.species.nms[,2:1],pch=3, xlab="NMS 2",ylab="NMS 1")
#color species by relative abundance!
#red maple is everywhere so it should be on here too
points(coweeta.species.nms[,2:1],pch=19,col="darkgreen", cex=6.0*coweeta.species.rel$ACRU)
#tulip seems to sort nice so we can add that in too
points(coweeta.species.nms[,2:1],pch=19,col="orange", cex=6.0*coweeta.species.rel$LITU)
#pitch pine sorts well with the soil gradient component so let's do that too.
points(coweeta.species.nms[,2:1],pch=19,col="blue", cex=6.0*coweeta.species.rel$PIRI)
#add red oak in to fill out the corner
points(coweeta.species.nms[,2:1],pch=19,col="lightgreen", cex=6.0*coweeta.species.rel$QURU)
#chestnut oak is very common, so it's worth putting that in as well
points(coweeta.species.nms[,2:1],pch=19,col="cyan", cex=6.0*coweeta.species.rel$QUPR)
#great laurels are everywhere
points(coweeta.species.nms[,2:1],pch=19,col="darkred", cex=6.0*coweeta.species.rel$RHMA)
#black locust illustrates the elevation gradient on NMS 2, so we want that
points(coweeta.species.nms[,2:1],pch=19,col="purple", cex=6.0*coweeta.species.rel$ROPS)
#add in the legend names and colors:
legend.txt <- c("ACRU","LITU","PIRI","QURU","QUPR", "RHMA","ROPS")
legend.col <- c("darkgreen","orange","blue","lightgreen","cyan","darkred","purple")
legend("topright", horiz=F,legend=legend.txt,text.col=legend.col,bty="n")
box(lwd=2)

#that's all folks

################################################################################

# Eni's stuff

# example code

# spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
# spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
# 
# p <- ggplot(scrs) +
#   geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Group)) +
#   coord_fixed() + ## need aspect ratio of 1!
#   geom_segment(data = spp.scrs,
#                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#                arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#   geom_text(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species),
#             size = 3)

# make the environmental vector into a data frame with only the statistically significant vectors
coweeta.species.nms.vf

NMS1 = coweeta.species.nms.vf[,2]
NMS2 = coweeta.species.nms.vf[,1]
pval = coweeta.species.nms.vf[,4]
r2 = coweeta.species.nms.vf[,3]

coweeta.species.nms.vf.df = data.frame(NMS1 = NMS1, NMS2 = NMS2, p.val = pval, r.2 = r2) %>% 
  mutate(env.var = rownames(.), NMS1.scaled = NMS1 * r2, NMS2.scaled = NMS2 * r2) %>% 
  filter(p.val <= 0.01) 
coweeta.species.nms.vf.df

########################### SOIL MOISTURE GRADIENT ########################### 

# get species labels, positioned according to their weighted average vectors

# get just the important species from the weighted average dataset (these will be labelled)
coweeta.species.nms.wa.df.sm = as.data.frame(coweeta.species.nms.wa) %>% 
  filter(row.names(.) %in% c('ACRU', 'QUCO', 'KALA', 'PIRI', 'BELE', 'LITU'))

# get dominant trees + their BA for each plot (for nicer points on the graph)

# get just the important species from the relativized species data set
coweeta.species.rel.subset.sm = coweeta.species.rel %>% 
  select(ACRU, QUCO, KALA, PIRI, BELE, LITU)

# make nms result into a dataframe and bind the relativized dataset onto it
coweeta.species.nms.df.sm = as.data.frame(coweeta.species.nms) %>%
  cbind(., coweeta.species.rel.subset.sm)

# calculate the species with the highest BA in each plot 
highest_vals = c()
highest_specs = c()

for (i in 1:nrow(coweeta.species.nms.df.sm)) { 
  values = coweeta.species.nms.df.sm[i, -c(1,2)]
  # print(values)
  # print(colnames(values))
  highest_val = apply(values[,-1], 1, max)
  highest_vals = c(highest_vals, highest_val)
  # print(colnames(values))
  # highest_spec= names(values == highest_val)
  highest_spec = names(values)[apply(values, 1, function(x) which(x == highest_val))][1]
  # highest_spec = lapply(apply(values, 2, function(x)which(x==highest_val)), names)
  # print(highest_spec)
  highest_specs = c(highest_specs, highest_spec)
}

# add to dataframe 
coweeta.species.nms.df.sm$highest_vals = highest_vals
coweeta.species.nms.df.sm$highest_specs = highest_specs

# graph everything in ggplot

ggplot() +
  
  # lines to show origin
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  
  # sites colored for dominant tree species + weighted for BA
  geom_point(data = coweeta.species.nms.df.sm, 
             aes(x = NMS2, y = NMS1, 
                 color = highest_specs,
                 size = highest_vals)) +
  
  # species name labels positioned according to weighted average vectors
  geom_text(data = coweeta.species.nms.wa.df.sm, 
            aes(x = NMS2, 
                y = NMS1, 
                label = row.names(coweeta.species.nms.wa.df.sm)),
            alpha = 0.8, size = 3) +
  
  # # environment vectors
  # geom_segment(data = coweeta.species.nms.vf.df,
  #              aes(x = 0, xend = NMS2.scaled, y = 0, yend = NMS1.scaled),
  #              arrow = arrow(length = unit(.5, "cm")), 
  #              colour = "black",
  #              alpha = 0.2,
  #              size = 2) +
  
  # # environment vector labels 
  # geom_text(data = coweeta.species.nms.vf.df, 
#           aes(x = NMS2.scaled, 
#               y = NMS1.scaled, 
#               label = row.names(coweeta.species.nms.vf.df)),
#           alpha = 0.8, 
#           size = 3) +

# other fomratting stuff
labs(color = "Species Code:", size = expression(paste('Basal Area (', "ft"^2, ')'))) +
  xlim(-.8,.8) +
  ylim(-.8,.8) +
  theme_gray() +
  theme(legend.position = "right",
        text = element_text(size = 12))

########################### ELEVATION GRADIENT ########################### 

# get species labels, positioned according to their weighted average vectors

# get just the important species from the weighted average dataset (these will be labelled)
coweeta.species.nms.wa.df.e = as.data.frame(coweeta.species.nms.wa) %>% 
  filter(row.names(.) %in% c("ACRU","QUAL","COFL","RHMA","ROPS"))

# get dominant trees + their BA for each plot (for nicer points on the graph)

# get just the important species from the relativized species data set
coweeta.species.rel.subset.e = coweeta.species.rel %>% 
  select(ACRU,QUAL,COFL,RHMA,ROPS)

# make nms result into a dataframe and bind the relativized dataset onto it
coweeta.species.nms.df.e = as.data.frame(coweeta.species.nms) %>%
  cbind(., coweeta.species.rel.subset.e)

# calculate the species with the highest BA in each plot 
highest_vals = c()
highest_specs = c()

for (i in 1:nrow(coweeta.species.nms.df.e)) { 
  values = coweeta.species.nms.df.e[i, -c(1,2)]
  # print(values)
  # print(colnames(values))
  highest_val = apply(values[,-1], 1, max)
  highest_vals = c(highest_vals, highest_val)
  # print(colnames(values))
  # highest_spec= names(values == highest_val)
  highest_spec = names(values)[apply(values, 1, function(x) which(x == highest_val))][1]
  # highest_spec = lapply(apply(values, 2, function(x)which(x==highest_val)), names)
  # print(highest_spec)
  highest_specs = c(highest_specs, highest_spec)
}

# add to dataframe 
coweeta.species.nms.df.e$highest_vals = highest_vals
coweeta.species.nms.df.e$highest_specs = highest_specs

# graph everything in ggplot

ggplot() +
  
  # lines to show origin
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  
  # sites colored for dominant tree species + weighted for BA
  geom_point(data = coweeta.species.nms.df.e, 
             aes(x = NMS2, y = NMS1, 
                 color = highest_specs,
                 size = highest_vals)) +
  
  # species name labels positioned according to weighted average vectors
  geom_text(data = coweeta.species.nms.wa.df.e, 
            aes(x = NMS2, 
                y = NMS1, 
                label = row.names(coweeta.species.nms.wa.df.e)),
            alpha = 0.8, size = 3) +
  
  # # environment vectors
  # geom_segment(data = coweeta.species.nms.vf.df,
  #              aes(x = 0, xend = NMS2.scaled, y = 0, yend = NMS1.scaled),
  #              arrow = arrow(length = unit(.5, "cm")), 
  #              colour = "black",
  #              alpha = 0.2,
  #              size = 2) +
  
  # # environment vector labels 
  # geom_text(data = coweeta.species.nms.vf.df, 
#           aes(x = NMS2.scaled, 
#               y = NMS1.scaled, 
#               label = row.names(coweeta.species.nms.vf.df)),
#           alpha = 0.8, 
#           size = 3) +

# other fomratting stuff
labs(color = "Species Code:", size = expression(paste('Basal Area (', "ft"^2, ')'))) +
  xlim(-.8,.8) +
  ylim(-.8,.8) +
  theme_gray() +
  theme(legend.position = "right",
        text = element_text(size = 12))

########################### ENVIRONMENT VECTORS ########################### 

coweeta.species.nms.df.full = as.data.frame(coweeta.species.nms) %>%
  cbind(., coweeta.species.rel)

ggplot() +
  
  # lines to show origin
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  
  # sites colored for dominant tree species + weighted for BA
  geom_point(data = coweeta.species.nms.df.full, 
             aes(x = NMS2, y = NMS1)) +
  
  # # species name labels positioned according to weighted average vectors
  # geom_text(data = coweeta.species.nms.wa.df, 
  #           aes(x = NMS2, 
  #               y = NMS1, 
  #               label = row.names(coweeta.species.nms.wa.df)),
  #           alpha = 0.8, size = 3) +
  
  # environment vectors
  geom_segment(data = coweeta.species.nms.vf.df,
               aes(x = 0, xend = NMS2.scaled, y = 0, yend = NMS1.scaled, colour = env.var),
               arrow = arrow(length = unit(.5, "cm")),
               alpha = 0.8,
               size = 2) +
  
  # environment vector labels
  geom_text(data = coweeta.species.nms.vf.df,
            aes(x = NMS2.scaled,
                y = NMS1.scaled,
                label = row.names(coweeta.species.nms.vf.df)),
            alpha = 1,
            size = 3) +
  
  # other fomratting stuff
  labs(color = "Environmental Variable:") +
  xlim(-.8,.8) +
  ylim(-.8,.8) +
  theme_gray() +
  theme(legend.position = "right",
        text = element_text(size = 12))

# example code

# geom_point(data = coweeta.species.nms.df, aes(x = NMS1, y = NMS2), 
#            color = "green", 
#            size = 6.0*coweeta.species.nms.df$"ACRU") +
# geom_point(data = coweeta.species.nms.df, aes(x = NMS1, y = NMS2), 
#            color = "red", 
#            size = 6.0*coweeta.species.nms.df$"PIRI") +
# scale_color_manual(values = inferno(15)[c(3, 8, 11)],
#                    name = "Aquatic System Type") +
# annotate(geom = "label", x = -1, y = 1.25, size = 10,
#          label = paste("Stress: ", round(nmds_results$stress, digits = 3))) +
