

# Script adapted from  Paniw et al. 2018 Ecol Lett - Appendix S2 - PCA analyses

# This script is divided into three parts:
# PART A : Create tree from phylogenies obtained from open sources
# PART B : Perform a phylogenetically-informed PCA on the life history quantities
# PART C : Create tree from phylogenies obtained from open sources

# Created: Feb 2021


################################################################################
############################  PART A ########################################  
################################################################################


#----------------------------------Clean memory---------------------------------
rm(list=ls(all=TRUE))

#---------------------------Setting working directory---------------------------
setwd("LOCATION")
getwd()


#-------------------------------Loading libraries-------------------------------
library(rotl)
library(taxize)
library(ape)
library(readxl)
library(phytools) 
library(ggplot2)  
library(ggrepel)  

#-----------Reading population quantities on which PCA is performed-------------

PopQs<-read_excel("Good_environment.xlsx")
# Reading in an excel file with in the first column the species names: make sure that 
# in your Excel file tthere is an underscore _ between the genus and species name (e.g. Cephalopholis_fulva)
# make sure that the names match those in the tree. You can check this by comparing
# smalltree$tip.label with PopQs$Species . This is also done below in the code.
# the other colums are the life history traits (DEB-IPM pars) 
# and derived life history traits (lambda, R0, generation time)

################################################################################
############################ PART A - BUILD THE TREE ###########################
################################################################################

#---------------------------Loading taxa to create tree-------------------------

a <- read.csv("LOCATION\\ott_names3.txt",sep = '\t') #Taxa that are used
a <- a[order(a$Scientific.name),] #set species names to alphabetic order

list <- rep(NA,dim(a)[1])

for (val in 1:dim(a)[1]) {
  list[val] <- a['OTT.name'][val,1]
}

taxa <- tnrs_match_names(names = c(list))

#---------------------------------Small tree------------------------------------

smalltree <- tol_induced_subtree(ott_ids = ott_id(taxa), label_format = "name")
plot(smalltree, cex = .8, label.offset = .1, no.margin = TRUE)

# set labels
smalltree$node.label<-as.character(1:length(smalltree$node.label))

# Drop all the tips in the tree that are not in the data
drop<-smalltree$tip.label[which(!smalltree$tip.label%in%PopQs$Species)]

#Should be TRUE
length(drop)+length(unique(PopQs$Species))==length(smalltree$tip.label)

# Trim tree of not needed species, if any at all
smalltree<-drop.tip(smalltree, drop)

# sorting required for PCA
all(sort(PopQs$Species)==sort(smalltree$tip.label))

#-----------------------------Editing the small tree----------------------------

tree=smalltree
tree$node.label<-as.character(1:length(tree$node.label))

#The tree must have branch lengths to be able to build the phylogenetic 
#variance-covariance matrix into the function phyl.pca for the phylogenetically-informed PCA
tree <- compute.brlen(tree)

# To use the tree in an MCMC analysis, it needs to be rooted
# By default the tree should be unrooted, but we can root it
# using the command root or root.tree
tree <- root(tree, tree$tip.label[1])

tree <- multi2di(tree, random=FALSE)  
is.rooted(tree)   #This is fo checking if the tree became rooted 

any(duplicated(tree$node.label))   
tree$node.label <- unique(tree$node.label)

#Another possible issue could be that the distance between taxa
#Could be 0, and this is not enabled for some analyses (including MCMCglmm)
#Thus we can solve this by adding an 0.001 to all tips
tree$edge.length <-  tree$edge.length + 0.001

#-------------------Transform the tree into an ultametric one-------------------

is.ultrametric(tree) 
tree <- chronos(tree, lambda=0, model="correlated") 

# continue
is.ultrametric(tree)
class(tree) <- "phylo"

any(duplicated(tree$node.label))
tree$node.label <- unique(tree$node.label)

#-----------------------Plot the final, ultametric tree-------------------------

plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)


################################################################################
#################### PART B PHYLOGENETICALLY - INFORMED PCA ####################
################################################################################


#-----------------------------Phylo PCA-----------------------------------------

row.names(PopQs)=PopQs$Species #This shows where the row names can be found

#Standardizing the data
pcaData1=log(PopQs[,c(2:8)])
pcaData=scale(pcaData1) # here non-standardised data are standardised.
# pcaData=PopQs[,c(2:4)] # use this if the data are standardised.

summary(pcaData)


#Checking the means and standard dev
colMeans(pcaData)
apply(pcaData, 2, sd)


#PCA
pca=phyl.pca(tree,pcaData,method="lambda",mode="corr") 

#Some results of the PCA
summary(pca)
pca$lambda 
pca$logL

#-------------------PCA is ready at this point----------------------------------

diag(pca$Eval)^2 # Apply Kaiser's criterion to determine how many PCA axes to keep

ncomp=2 # keep the first 2 axes

# Varimax correction on the axes
rawLoadings <- pca$L[,1:ncomp]

# find the variance maximizing rotation of loadings      
rotatedLoadings <- varimax(rawLoadings)$loadings

name1= summary(pca)$importance[[2,1]]   # fraction variance explained by first axis; this output is taken from summary(pca) 
name2= summary(pca)$importance[[2,2]]   # fraction variance explained by second axis; this output is taken from summary(pca)

# Create the inverse of the loading to calculate new scores: data multiplied by rotation matrix 
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- pcaData %*% invLoadings


x <- list() 
x$scores <- scores[,1:2]
colnames(x$scores)=c("PC1", "PC2")
x$scores[,1:2] <- (-1)*x$scores[,1:2]
x$loadings <- rotatedLoadings[,1:2]
colnames(x$loadings)=c("PC1", "PC2")
x$loadings[,1:2] <- (-1)*x$loadings[,1:2]

#-----------------Save the PCA scores for regression analyses-------------------

write.csv(cbind(x$scores,PopQs$Species),"PCAscores.csv",row.names=F)

  
  


################################################################################
#################### PART C Plotting the result of the PCA #####################
################################################################################



### With point labels (species ID)
data1 <- data.frame(x$scores,point.lab=1:nrow(x$scores))

options(ggrepel.max.overlaps = Inf)


PCbiplot <- function(PC, x="PC1", y="PC2") {
  
  #position on the plot
  plot <- ggplot(data1, aes_string(x=x, y=y)) 
  
  #editing the data points
  plot <- plot + geom_point(aes(color = PopQs$Species), size=4.5, alpha=0.5)
  
  #adding lines to the plot
  plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
  
  plot <- plot + scale_size(range = c(0, 6))    
  
  datapc <- data.frame(varnames=colnames(PopQs[2:8]), PC$loadings)    
  
  mult <- min(
    (max(data1[,y]) - min(data1[,y])/(max(datapc[,y])-min(datapc[,y]))), 
    (max(data1[,x]) - min(data1[,x])/(max(datapc[,x])-min(datapc[,x])))  
  )
  
  datapc <- transform(datapc, 
                      v1 = .50 * mult* (get(x)),      #I dont understand this line 
                      v2 = .50 * mult* (get(y))       #I dont understand this line 
  )
  
  #Adding the numbers and the title
  plot <- plot + geom_text_repel(aes(label=point.lab)) + labs(title = "High feeding level") + theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"))
  
  #adding the 3 arrows
  plot <- plot + geom_segment(data=datapc[1:1,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[1:1,], aes(x=1.2,y=-0.08,label = "\U03BB", parse = TRUE), size=5)
  plot <- plot + geom_segment(data=datapc[2:2,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[2:2,], aes(x=1.18,y=0.2,label = "R[0]"), parse = TRUE, size = 5)
  plot <- plot + geom_segment(data=datapc[3:3,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[3:3,], aes(x=-1.25,y=-0.08,label = "GT"), size = 5)
  plot <- plot + geom_segment(data=datapc[4:4,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[4:4,], aes(x=0.1,y=1.25,label = "L[m]"), parse = TRUE, size = 5)
  plot <- plot + geom_segment(data=datapc[5:5,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[5:5,], aes(x=-0.2,y=1.25,label = "L[p]"), parse = TRUE, size = 5)
  plot <- plot + geom_segment(data=datapc[6:6,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[6:6,], aes(x=0.13,y=1.40,label = "L[b]"), parse = TRUE, size = 5)
  plot <- plot + geom_segment(data=datapc[7:7,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=0.8, color="black") + geom_text(data = datapc[7:7,], aes(x=1.23,y=-0.2,label = "r[B]"), parse = TRUE, size = 5)
  
  
  #Axis settings
  plot <- plot+ylab(paste("PCA 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PCA 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=14))+theme(axis.title = element_text(size=16))
  
  #Legend settings
  plot<- plot + labs(color="Species") + theme(legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=14),legend.key.size = unit(1, "lines"))
  plot<- plot + scale_color_hue(labels = c(a$Scientific.name), l=60)
  
  plot
}

PCbiplot(x)

