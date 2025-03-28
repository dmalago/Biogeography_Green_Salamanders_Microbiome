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


set.seed(100)
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

wtest<-adonis2(adists~cladal)
posthoctest<-pairwise.adonis(adists,as.factor(cladal),perm=99999)

pcoa_adists<-pcoa(adists)$vectors[,1:2]
types<-c('North Appalachian','HNG','BRE','Southern Appalachian')
colvec<-c('yellow','green','blue','purple')
no_type<-cladal
no_type[cladal=='North Appalachian']<-1
no_type[cladal=='HNG']<-2
no_type[cladal=='BRE']<-3
no_type[cladal=='Southern Appalachian']<-4
no_type<-as.numeric(no_type)
type<-cladal

plot(pcoa_adists,type='n',cex.lab=1.5,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_adists,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_adists, display = "sites", pch=16, col = colvec[no_type])

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_adists<-c()    #Value of the centroid along the first NMDS axis
ymeans_adists<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_adists<-c(xmeans_adists,mean(pcoa_adists[pts,1]))
  ymeans_adists<-c(ymeans_adists,mean(pcoa_adists[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_adists[k],ymeans_adists[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_adists[k],ymeans_adists[k],col='black',pch=1,cex=2)
}

testdisp<-betadisper(as.dist(adists),cladal)
dispersion_anova<-anova(testdisp)



