########Loading Required Libraries#######
library("phyloseq")
library("ggplot2")  #Used for Graphics
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library(ape)
library('geodist')
library(picante)
library(dplyr)
library(MASS)
library('adespatial')
library('glmmTMB')
library(adespatial);library(sp);library(spdep)
library('ade4')
library('geostats')
library('geosphere')

set.seed(10)
testmet<-'unifrac'

####### Read in Saved Phyloseq object which has been processed and rarified#####
data_rarified<-readRDS("Greens_Biogeography_Phydata")

#Remove Environmental ASVs
#This function will be used to remove the bacteria also found in your environment samples
prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs) 
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps) 
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

#Defining the samples that are environment
Envneg.cts = subset_samples(data_rarified, Sample_Type != "Salamander")
Envneg.cts<- prune_taxa(taxa_sums(Envneg.cts) > 1, Envneg.cts) 
#Defining the samples that are salamanders

Sal_samples = subset_samples(data_rarified, Sample_Type == "Salamander")
Sal_samples<- prune_taxa(taxa_sums(Sal_samples) > 0, Sal_samples) 

#Creating a new Phyloseq object with bacteria not found in the environment
data_rarified = prune_negatives(data_rarified,Envneg.cts,Sal_samples)


Salamanders<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

clim<-read.csv('bioclim.csv')
climhit<-c("wc2_1_30s_bio_18", "wc2_1_30s_bio_9",  "wc2_1_30s_bio_8",  "wc2_1_30s_bio_7",  "wc2_1_30s_bio_3" )
climf<-scale(clim[,which(colnames(clim) %in% climhit)])
rownames(climf)<-clim[,1]
colnames(climf)<-c('PWQ','TDQ','TWQ','TR','IS')

if (testmet=='euclidean'){
#metric<-'euclidean'
edist<-vegdist(t(otu_table(Salamanders)),method='euclidean',upper=TRUE,diag=TRUE, binary=FALSE)
edist<-as.matrix(edist)
}else if (testmet=='jaccard'){
#metric<-'jaccard'
edist<-vegdist(sign(t(otu_table(Salamanders))),method='jaccard',upper=TRUE,diag=TRUE, binary=FALSE)
edist<-as.matrix(edist)
}else if (testmet=='bray'){
#metric<-'bray'
edist<-vegdist((t(otu_table(Salamanders))),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
edist<-as.matrix(edist)
}else if (testmet=='unifrac'){
metric<-'unifrac'
edist <- as.matrix(UniFrac(Salamanders,weighted=FALSE))
}else{
metric<-'weighted unifrac'
edist <- as.matrix(UniFrac(Salamanders,weighted=TRUE))
}
sitelist<-unique(sample_data(Salamanders)$Study.Site)


adists<-matrix(0,length(sitelist),length(sitelist))
rownames(adists)<-sitelist
colnames(adists)<-sitelist

cladal<-c()
for (k in 1:length(sitelist)){
  hit1<-which(sample_data(Salamanders)$Study.Site==sitelist[k])
  cladal<-c(cladal,sample_data(Salamanders)$Clade[hit1[1]])
  for (j in 1:length(sitelist)){
    hit2<-which(sample_data(Salamanders)$Study.Site==sitelist[j])
    
    if (k==j){
      p1<-which(rownames(adists)==sitelist[k])
      adists[p1,p1]<-0#mean(edist[hit1,hit2])
    }
    else{
      p1<-which(rownames(adists)==sitelist[k])
      p2<-which(rownames(adists)==sitelist[j])
      adists[p1,p2]<-mean(edist[hit1,hit2])
      adists[p2,p1]<-mean(edist[hit1,hit2])
    }
  }
  
}



rdaout<-dbrda(adists ~ PWQ + TDQ +TWQ+TR+IS, data = data.frame(climf))
plotter<-scores(rdaout, display="sites", choices=c(1,2), scaling=1)
ordiplot(rdaout,cex=1.5,display = 'sites',cex.lab=1.5,scaling=1,xlim=c(-0.5,0.5))
text(rdaout, display="bp", col='black',cex=1)
points(plotter[which(cladal=='HNG'),1],plotter[which(cladal=='HNG'),2],pch=16,col='green',cex=1.5)
points(plotter[which(cladal=='BRE'),1],plotter[which(cladal=='BRE'),2],pch=16,col='blue',cex=1.5)
points(plotter[which(cladal=='North Appalachian'),1],plotter[which(cladal=='North Appalachian'),2],pch=16,col='yellow',cex=1.5)
points(plotter[which(cladal=='Southern Appalachian'),1],plotter[which(cladal=='Southern Appalachian'),2],pch=16,col='purple',cex=1.5)

anova.cca(rdaout, step = 1000)
anova.cca(rdaout, step = 1000,by='term')
anova.cca(rdaout, step = 1000,by='axis')

