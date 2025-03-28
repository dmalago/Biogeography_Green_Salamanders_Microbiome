#############Biogeography Greens Mantel Tests Code ###########

setwd("~/Desktop/Biogeography Greens/Github_Upload")



########Loading Required Libraries#######

library(lme4)
library("phyloseq")
library("ggplot2")  #Used for Graphics
library("VennDiagram") #Used for Venn Diagram
library(venneuler) #Used for Venn Diagram
library(gplots) #Used for Graphics
library(limma)
library(vegan)
library("metagMisc")
library("qiime2R")
library(devtools)
library(car)
library(metagenomeSeq)
library(tidyr)
library(dplyr)
library('stringr')
library('microshades')
library('speedyseq')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
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
library(lmerTest)
library("car")
library(cAIC4)
library("corrplot")
library(patchwork)
library("devtools")




set.seed(1)
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




data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
TAX<-tax_table(data_rarified)
#Clade Color Palette
#scale_color_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))


####Mantel Tests ####


clim<-read.csv('bioclim.csv')
climhit<-c("wc2_1_30s_bio_18", "wc2_1_30s_bio_9",  "wc2_1_30s_bio_8",  "wc2_1_30s_bio_7",  "wc2_1_30s_bio_3" )
climf<-scale(clim[,which(colnames(clim) %in% climhit)])
rownames(climf)<-clim[,1]
colnames(climf)<-c('PWQ','TDQ','TWQ','TR','IS')
rownames(climf)[which(rownames(climf)=='Cumberland Mountan Trail')]<-'Cumberland Mountain Trail'
rownames(climf)[which(rownames(climf)=='Campsite near World\'s Edge')]<-'Campsite at World\'s Edge'


PWQ<-c()
TDQ<-c()
TWQ<-c()
TR<-c()
IS<-c()

for (k in 1:length(sample_data(data_rarified)$Study.Site)){
  temp<-which(rownames(climf)==sample_data(data_rarified)$Study.Site[k])
  PWQ<-c(PWQ,climf[temp,1])
  TDQ<-c(TDQ,climf[temp,2])
  TWQ<-c(TWQ,climf[temp,3])
  TR<-c(TR,climf[temp,4])
  IS<-c(IS,climf[temp,5])
}

envmat<-data.frame(PWQ,TDQ,TWQ,TR,IS)
rownames(envmat)<-sample_names(data_rarified)


dist_matrix<-as.matrix(vegdist(envmat,method='euclidean',diag=TRUE,upper=TRUE))
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)

sitelist<-unique(sample_data(data_rarified)$Study.Site)
bigdiff<-0
listbig<-c()
first<-1
for (k in 1:length(sitelist)){
  temp<-which(sample_data(data_rarified)$Study.Site==sitelist[k])
  bigdiff<-max(bigdiff,max(dist_matrix[temp,temp]))
  listbig<-c(listbig,as.vector(dist_matrix[temp,temp]))
  if (bigdiff>1000 && first == 1){
    print(sitelist[k])
    print(temp)
    first<-0
  }
}


#Find  jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(data_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(data_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(data_rarified,weighted=FALSE)
wunifrac<-UniFrac(data_rarified,weighted=TRUE)


#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)

#The mantel tests looks to see whether there is a correlation between the two matrices (in our case, this is the matrix of spatial distances and the matrix of differences/distances between communities in 'community space')
#Again, you have options for the metric of correlation you use. I've chosen Pearson correlation

mantel_test_jaccard  <- mantel(jaccard_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_jaccard

mantel_test_bray  <- mantel(bray_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_bray

mantel_test_unifrac  <- mantel(unifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_unifrac

mantel_test_wunifrac  <- mantel(wunifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_wunifrac






#Set Distance classes for correlogram analyses

#You can also tell the program to go to all the distances, even if you have small numbers of replicates at those distances
temp<-mantel.correlog(jaccard_matrix,dist_matrix,nperm=9999)
mantel_correlog_jaccard<-mantel.correlog(jaccard_matrix,dist_matrix,nperm=9999, break.pts=c(0,0.01,temp$break.pts[2:length(temp$break.pts)]))
plot(mantel_correlog_jaccard, alpha=0.05)

#You can also tell the program to go to all the distances, even if you have small numbers of replicates at those distances
temp<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999)
mantel_correlog_bray<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999, break.pts=c(0,0.01,temp$break.pts[2:length(temp$break.pts)]))
plot(mantel_correlog_bray, alpha=0.05)

#You can also tell the program to go to all the distances, even if you have small numbers of replicates at those distances
temp<-mantel.correlog(unifrac_matrix,dist_matrix,nperm=9999)
mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,dist_matrix,nperm=9999, break.pts=c(0,0.01,temp$break.pts[2:length(temp$break.pts)]))
plot(mantel_correlog_unifrac, alpha=0.05,cex.lab=2)

#You can also tell the program to go to all the distances, even if you have small numbers of replicates at those distances
temp<-mantel.correlog(wunifrac_matrix,dist_matrix,nperm=9999)
mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,dist_matrix,nperm=9999, break.pts=c(0,0.01,temp$break.pts[2:length(temp$break.pts)]))
plot(mantel_correlog_wunifrac, alpha=0.05)


uni_df<-as.data.frame(mantel_correlog_unifrac$mantel.res)
sig<-c()
for (k in 1:length(mantel_correlog_unifrac$mantel.res[,5])){
  if (is.na(mantel_correlog_unifrac$mantel.res[k,5])){
    sig<-c(sig,100)
  }else{
  if (mantel_correlog_unifrac$mantel.res[k,5]<0.05){
    sig<-c(sig,1)
  }else{
    sig<-c(sig,0)
  }
  }
}

uni_df<-data.frame(uni_df,sig)
uni_df<-uni_df[-which(uni_df$sig==100),]
uni_df$sig<-as.factor(uni_df$sig)
color_palette <- c("white", "black")

######Unifrac Mantel Geographic plot#####
pdf(file = "Unifrac_Mantel_Geographic.pdf", width = 10, height = 10)
ggplot(uni_df, aes(x = class.index, y = Mantel.cor, fill = sig, group=1)) +
  geom_line()+
  geom_point(size = 7, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = "Environmental Distance Class", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=20))+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


dev.off()



