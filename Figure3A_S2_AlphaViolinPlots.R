#############Biogeography Greens Figure X Code ###########

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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#####Calculate a variety of different ALPHA diversity metrics####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Find the richness for each sample
richness<-colSums(sign(otu_table(data_rarified)))

#Find the shannon diversity for each sample
shannon<-diversity(t(otu_table(data_rarified)),index = 'shannon')

#Find Faith's pd for each sample
faithpd<-pd(t(otu_table(data_rarified)), phy_tree(data_rarified), include.root=TRUE)
faiths<-faithpd[,1]

Clade<-sample_data(data_rarified)$Clade

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(Clade,richness,shannon,faiths)
#Switch the ordering of the types (to the order you want them presented in your plots)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Richness
k_richness<-kruskal.test(richness ~ Clade, data = diversity_df)
k_richness

#shannon Entropy
k_shannon<-kruskal.test(shannon ~ Clade, data = diversity_df)
k_shannon

#Faith's PD
k_faiths<-kruskal.test(faiths ~ Clade, data = diversity_df)
k_faiths

# Define a common theme for all plots
common_theme <- theme_classic() +
  theme(text = element_text(size=20),  # Sets a good size for readability
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = "none",  # Hides the legend
        plot.title = element_text(size=0, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18))

# Define colors
fill_colors <- c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple")

# Plot for ASV Richness
a <- ggplot(diversity_df, aes(x = reorder(Clade, richness), y = richness, fill = Clade)) + 
  geom_violin(width = .75) +
  stat_summary(fun="mean", color = "black", size= 1.25) +
  labs(title = "ASV Richness", x = "Clade", y = "ASV Richness") + 
  geom_point(position = position_jitter(seed = 1, width = 0.25), color = "grey") +
  scale_fill_manual(values = fill_colors) + scale_x_discrete(labels = c("North Appalachian" = "N. Appalachian","Southern Appalachian" = "S. Appalachian",'BRE' = 'BRE','HNG' = 'HNG'))+
  common_theme

# Plot for Shannon Diversity
b <- ggplot(diversity_df, aes(x = reorder(Clade, shannon), y = shannon, fill = Clade)) + 
  geom_violin(width = .75) +
  stat_summary(fun="mean", color = "black", size= 1.25) +
  labs(title = "Shannon Diversity", x = "Clade", y = "Shannon Diversity") + 
  geom_point(position = position_jitter(seed = 1, width = 0.25), color = "grey") +
  scale_fill_manual(values = fill_colors) + scale_x_discrete(labels = c("North Appalachian" = "N. Appalachian","Southern Appalachian" = "S. Appalachian",'BRE' = 'BRE','HNG' = 'HNG'))+
  common_theme

# Plot for Faith's PD
c <- ggplot(diversity_df, aes(x = reorder(Clade, faiths), y = faiths, fill = Clade)) + 
  geom_violin(width = .75) +
  stat_summary(fun="mean", color = "black", size= 1.25) +
  labs(title = "Faith's PD", x = "Clade", y = "Faith's PD") + 
  geom_point(position = position_jitter(seed = 1, width = 0.25), color = "grey") +
  scale_fill_manual(values = fill_colors) + scale_x_discrete(labels = c("North Appalachian" = "N. Appalachian","Southern Appalachian" = "S. Appalachian",'BRE' = 'BRE','HNG' = 'HNG'))+
  common_theme

# Print plot c
print(c)

# Save the plots to a PDF
pdf(file = "combined_alpha_violin_plots.pdf", width = 15, height = 15)
print(a + b + c + plot_layout(ncol = 1))
dev.off()




