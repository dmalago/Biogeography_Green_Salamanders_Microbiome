#############Biogeography Greens Antifungal Code ###########

setwd("~/Desktop/Biogeography Greens/Antifungal_Analyses")


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




#Clade Color Palette
#scale_color_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))





Antifungal_data<-read.csv("Antifungal_data.csv")
##### Alpha Diversity Analyses ####

Antifungal_data<-Antifungal_data[which(Antifungal_data$Sample_Type!= "Crevice"),]


#This is just so you can see the overall pattern across all sites. No recognizable pattern I can discern...
a<-ggplot(Antifungal_data, aes(x = reorder(Clade, AntiFungal_Richness), y = AntiFungal_Richness, fill = Clade)) + 
  geom_violin(width = .5) + scale_fill_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))+
  stat_summary(fun="mean", color = "black", size= 2.5)+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs( 
    title = "",  
    x = "Population", 
    y = "Antifungal ASV Richness" 
  ) +
  geom_point(position = position_jitter(seed = 1, width = 0.15, ), color = "grey") +
  theme(text = element_text(size=30))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(text = element_text(size=25))+theme(legend.position = "none")


Richness_clade_kruskal_wallis<-kruskal.test(AntiFungal_Richness~Clade, data = Antifungal_data)
Richness_clade_kruskal_wallis



b<-ggplot(Antifungal_data, aes(x = reorder(Clade, Propor_TotalAntiFungal), y = Propor_TotalAntiFungal, fill = Clade)) + 
  geom_violin(width = .5) + scale_fill_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))+
  stat_summary(fun="mean", color = "black", size= 1.25)+
  labs( 
    title = "",  
    x = "Population", 
    y = "Proportion of Microbiota Chytrid-inhibitory" 
  ) + 
  geom_point(position = position_jitter(seed = 1, width = 0.15, ), color = "grey") +
  theme_classic()+theme_linedraw()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(text = element_text(size=25))
Propor_clade_kruskal_wallis<-kruskal.test(Propor_TotalAntiFungal~Clade, data = Antifungal_data)
Propor_clade_kruskal_wallis


pdf(file = "Antifungal_Clades.pdf", width = 30, height = 10)
a+b
dev.off()






Antifungal_data2<-subset(Antifungal_data, Bd !='?')

ggplot(Antifungal_data2, aes(x = Propor_TotalAntiFungal, y = as.numeric(Bd))) + 
  geom_point() +
  geom_smooth(method="glm", method.args=list(family="binomial"))+
  labs( 
    title = "",  
    x = "", 
    y = "" 
  )



Propor_clade_kruskal_wallis<-kruskal.test(Propor_TotalAntiFungal~Clade, data = Antifungal_data2)
Propor_clade_kruskal_wallis
Antifungal_data2<-subset(Antifungal_data2, Bd ==1)
Antifungal_data2<-subset(Antifungal_data2, Zoospore_Load !='#VALUE!')


#Antifungal_data2<-subset(Antifungal_data2, TotalAntiFungal !=1922)
ggplot(Antifungal_data2, aes(x=AntiFungal_Richness,y=as.numeric(Zoospores), group =1)) + 
  geom_point() +
  geom_line()+
  geom_smooth(method = "glm")+
  labs( 
    title = "",  
    x = "", 
    y = "" 
  ) 

Antifungal_data2$Zoospores<-as.numeric((Antifungal_data2$Zoospore_Load))
library(ggplot2)
library(ggpubr)  # Assuming gghistogram is from the ggpubr package

gghistogram(Antifungal_data2, 
            x = "Zoospores",  
            fill = "blue", 
            add = "mean", 
            add_density = TRUE) +
  theme(text = element_text(size = 25),  # Base text size
        axis.title = element_text(size = 25),  # Axis titles
        axis.text = element_text(size = 25),  # Axis text
        legend.title = element_text(size = 25),  # Legend title
        legend.text = element_text(size = 25))  # Legend text

mean(Antifungal_data2$Zoospores)

cor.test(Antifungal_data2$AntiFungal_Richness, as.numeric(Antifungal_data2$Zoospores))

CrevicevSalr<-ggplot(Antifungal_data, aes(Sample_Type,AntiFungal_Richness)) + 
  geom_violin(aes(fill =Sample_Type)) + 
  labs( 
    title = "",  
    x = "Location", 
    y = "ASV Richness" 
  ) + stat_summary(fun="mean", color = "black", size= 1.25)+
  geom_point(position = position_jitter(seed = 1, width = 0.25, ), color = "grey") +
  theme(text = element_text(size=30))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(text = element_text(size=25))




#####Do Antifungal communities differ across Clades?####

#           Perform Jaccard PERMANOVA on Antifungal ASVs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(vegan)
library(PERMANOVA)

data_rarified<-readRDS("Greens_Biogeography_Phydata")
library(phyloseq)
data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")

sample_names_physeq <- sample_names(data_rarified)
Antifungal_table = qiime2R::read_qza("Antifungal_Matches_99/clustered_table.qza")

Antifungal_table = as.data.frame(Antifungal_table$data)
# Convert the OTU table to a data frame (if it's not already)
otu_table_df <- as.data.frame(Antifungal_table)

# Filter the OTU table to keep only the columns that match the sample names
filtered_otu_table <- otu_table_df[, colnames(otu_table_df) %in% sample_names_physeq]


jaccard<-vegdist(t(filtered_otu_table),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)



# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(data_rarified))

# Ensure the OTU table has the correct sample names as columns
otu_samples <- colnames(filtered_otu_table)

# Subset the sample data to include only the samples present in the OTU table
matching_samples <- sample_data_df[rownames(sample_data_df) %in% otu_samples, ]

# Extract the clade information (assuming the clade information is in a column named 'clade')
clade_info <- matching_samples[, "Clade", drop = FALSE]


type<-clade_info$Clade


#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(Antifungal_table),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=1000)



#Here I am doing another PERMANOVA comparing antifungals with Unifrac distances
#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata 
data<-qza_to_phyloseq(
  features="Antifungal_Matches_99/clustered_table.qza",
  tree="Antifungal_Matches_99/rooted-tree.qza", metadata = "metadata.tsv")

#Plot the phylogenetic tree 
plot.phylo(phy_tree(data))


data<-subset_samples(data, Sample_Type == "Salamander")



set.seed(100)
testmet<-'unifrac'

####### Read in Saved Phyloseq object which has been processed and rarified#####


Salamanders<- data

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

plot(pcoa_adists,type='n',cex.lab=1.25,xlab='PCoA1',ylab='PCoA2')
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

#Faiths PD######




faithpd<-pd(t(otu_table(Salamanders)), phy_tree(Salamanders), include.root=TRUE)
faiths<-faithpd[,1]

Clade<-sample_data(Salamanders)$Clade

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(Clade,faiths)


a<-ggplot(diversity_df, aes(x = reorder(Clade, faiths), y = faiths, fill = Clade)) + 
  geom_violin(width = .5) + scale_fill_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))+
  stat_summary(fun="mean", color = "black", size= 2.5)+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs( 
    title = "",  
    x = "Population", 
    y = "Chytrid-inhibitory Faith's PD" 
  ) +
  geom_point(position = position_jitter(seed = 1, width = 0.15, ), color = "grey") +
  theme(text = element_text(size=30))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(text = element_text(size=25))+theme(legend.position = "none")


Richness_clade_kruskal_wallis<-kruskal.test(AntiFungal_Richness~Clade, data = Antifungal_data)
Richness_clade_kruskal_wallis

pdf(file = "Antifungal_Clades.pdf", width = 30, height = 10)
a+b
dev.off()

