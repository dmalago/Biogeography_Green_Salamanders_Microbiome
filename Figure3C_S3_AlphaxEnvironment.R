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

#Clade Color Palette
#scale_color_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))


SD<-sample_data(Salamanders)

faithpd<-pd(t(otu_table(Salamanders)), phy_tree(Salamanders), include.root=FALSE)
shannon<-diversity(t(otu_table(Salamanders)),index='shannon')

alpha_df<-data.frame(SD$Longitude,SD$Latitude,SD$Clade,SD$Subclade,SD$Study.Site,colSums(sign(otu_table(Salamanders))),faithpd[,1],shannon)
colnames(alpha_df)<-c('long','lat','clade','subclade','site','rich','faith','shannon')
alpha_df<-alpha_df[order(alpha_df$clade),]
alpha_df$site[alpha_df$site=='Cumberland Mountain Trail']<-'Cumberland Mountan Trail'
alpha_df$site[alpha_df$site=='Campsite at World\'s Edge']<-'Campsite near World\'s Edge'

clim<-read.csv('bioclim.csv')
climhit<-c("wc2_1_30s_bio_18", "wc2_1_30s_bio_9",  "wc2_1_30s_bio_8",  "wc2_1_30s_bio_7",  "wc2_1_30s_bio_3" )
climf<-scale(clim[,which(colnames(clim) %in% climhit)])
rownames(climf)<-clim[,1]

w1<-c()
w2<-c()
w3<-c()
w4<-c()
w5<-c()
for (k in 1:length(alpha_df$site)){
  temp<-which(rownames(climf)==alpha_df$site[k])
  w1<-c(w1,climf[temp,1])
  w2<-c(w2,climf[temp,2])
  w3<-c(w3,climf[temp,3])
  w4<-c(w4,climf[temp,4])
  w5<-c(w5,climf[temp,5])
}

#Medians
mfalpha<-c()
mralpha<-c()
msalpha<-c()
nalpha<-c()
calpha<-c()
for (k in 1:length(climf[,1])){
  templ<-which(alpha_df$site==rownames(climf)[k])
  mfalpha<-c(mfalpha,median(alpha_df$faith[templ]))
  mralpha<-c(mralpha,median(alpha_df$rich[templ]))
  msalpha<-c(msalpha,median(alpha_df$shannon[templ]))
  nalpha<-c(nalpha,alpha_df$site[templ[1]])
  calpha<-c(calpha,alpha_df$clade[templ[1]])
}

mdf<-data.frame(climf,mfalpha,mralpha,msalpha,nalpha,calpha)
colnames(mdf)<-c('w1','w2','w3','w4','w5','faith','rich','shannon','site','clade')
mdf<-mdf[-1,]

richmodel<-glm(rich~w1+w2+w3+w4+w5+clade,data=mdf)
bestrichmodel<-stepAIC(richmodel,direction='both')
print(bestrichmodel)
finalrichmodel<-glm(rich~w2+w3+w5+clade,data=mdf)
summary(finalrichmodel)
r<-plot_model(finalrichmodel,title='',axis.labels=c('S. Appalachian','N. Appalachian','HNG','Isothermality','Temperature Wettest Quarter','Temperature Driest Quarter'),show.values=TRUE,vline.color='grey')
r<-r+theme(text = element_text(size = 17))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
plot(r)

shannonmodel<-glm(shannon~w1+w2+w3+w4+w5+clade,data=mdf)
bestshannonmodel<-stepAIC(shannonmodel,direction='both')
print(bestshannonmodel)
finalshannonmodel<-glm(rich~w2,data=mdf)
summary(finalshannonmodel)
s<-plot_model(finalshannonmodel,title='',axis.labels=c('Temperature \nDriest Quarter'),show.values=TRUE,vline.color='grey')
s<-s+theme(text = element_text(size = 17))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
plot(s)


faithmodel<-glm(faith~w1+w2+w3+w4+w5+clade,data=mdf)
bestfaithmodel<-stepAIC(faithmodel,direction='both')
print(bestfaithmodel)
finalfaithmodel<-glm(faith~w3+w4+clade,data=mdf)
summary(finalfaithmodel)
f<-plot_model(finalfaithmodel,title='',axis.labels=c('S. Appalachian','N. Appalachian','HNG','Temperature Range','Temperature Wettest Quarter'),show.values=TRUE,vline.color='grey')
f<-f+theme(text = element_text(size = 20))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
plot(f)


