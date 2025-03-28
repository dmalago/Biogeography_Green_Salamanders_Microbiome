#############Biogeography Greens Indicator Species and ANCOM Code ###########

setwd("~/Desktop/Biogeography Greens/")


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


#Modified function for creating the custom legend
custom_legend2column <- function (mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", x = "Sample",
                                  y = "Abundance", legend_key_size = 0.4, legend_text_size = 6)
{
  if (is.null(mdf[[group_level]])) {
    stop("mdf 'group_level' does not exist")
  }
  
  if (is.null(mdf[[subgroup_level]])) {
    stop("mdf 'subgroup_level' does not exist")
  }
  
  if (is.null(cdf$hex)) {
    stop("cdf 'hex' does not exist")
  }
  
  col_name_group <- paste0("Top_", group_level)
  col_name_subgroup <- paste0("Top_", subgroup_level)
  
  group_level_names <- unique(cdf[[col_name_group]])
  
  for (i in 1:length(group_level_names))
  {
    if( i == 1)
    {
      complete_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
      tracer<-14#length(unique(mdf[which(mdf[,48]==group_level_names[i]),49]))+0
    }
    else if (i ==2){
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-15#length(unique(mdf[which(mdf[,48]==group_level_names[i]),52]))-1
      tracer<-tracer+new_height
      print(c(complete_height,new_height))
      
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
      
    }
    else if (i==3){
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-12#length(unique(mdf[which(mdf[,48]==group_level_names[i]),52]))-1
      tracer<-tracer+new_height
      print(c(complete_height,new_height))
      
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
      
    }
    else if (i==4)
    {
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-28#length(unique(mdf[which(mdf[,48]==group_level_names[i]),52]))+1
      tracer<-tracer+new_height
      print(c(complete_height,new_height))
      print(complete_legend)
      print(new_legend)
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
    }
    else{
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-19#length(unique(mdf[which(mdf[,48]==group_level_names[i]),52]))+1
      tracer<-tracer+new_height
      print(c(complete_height,new_height))
      print(complete_legend)
      print(new_legend)
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
      
    }
  }
  plot(complete_legend)
  complete_legend
}

individual_legend2 <- function (mdf,
                                cdf,
                                group_name,
                                col_name_group = "Top_Phylum",
                                col_name_subgroup = "Top_Genus",
                                x = "Sample",
                                y = "Abundance",
                                legend_key_size = 0.4,
                                legend_text_size = 6)
{
  select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
  select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
  
  select_plot <- ggplot(select_mdf,
                        aes_string(x = x, y = y, fill = col_name_subgroup, text = col_name_subgroup)) +
    geom_col( position="fill") +
    scale_fill_manual(name = group_name,
                      values = select_cdf$hex,
                      breaks = select_cdf[[col_name_subgroup]]) +
    theme(legend.justification = "left") +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.key.size = unit(legend_key_size, "lines"), text=element_text(size=legend_text_size))
  
  legend <- get_legend(select_plot)
}


####### Read in Saved Phyloseq object which has been processed and rarified#####
data_rarified<-readRDS("Greens_Biogeography_Phydata")

data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
TAX<-tax_table(data_rarified)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
colnames(TAX)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
TAX<-gsub('p__', '', TAX)
TAX<-gsub('c__', '', TAX)
TAX<-gsub('o__', '', TAX)
TAX<-gsub('f__', '', TAX)
TAX<-gsub('g__', '', TAX)
TAX<-gsub('s__', '', TAX)


for (k in 1:length(TAX[,6])){
  if (is.na(TAX[k,6])){
    if (!(is.na(TAX[k,5]))){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (!(is.na(TAX[k,4]))){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (!(is.na(TAX[k,3]))){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}
for (k in 1:length(TAX[,6])){
  if (TAX[k,6]=='NA'){
    if (TAX[k,5]!='NA'){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (TAX[k,4]!='NA'){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (TAX[k,3]!='NA'){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}

data_rarified<-phyloseq(otu_table(data_rarified),TAX,sample_data(data_rarified),phy_tree(data_rarified))

#Clade Color Palette
#scale_color_manual( values =  c("BRE"= "blue", "HNG" = "green",  "North Appalachian" = "yellow", "Southern Appalachian" ="purple"))

#####Taxa Barplots#####
#%%%%%%%%%%%%%%%%%%%%%%   MAKE FANCY BAR GRAPHS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum
summed_phyla<-rowMeans(data.frame(otu_table(tax_glom(data_rarified,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list<-data.frame(tax_table(tax_glom(data_rarified,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list<-phyla_list[order(-summed_phyla)]

#Prepare the OTU table in the correct format for the plotting program (this must be done with level 6 data including genera!)
mdf_prep <- prep_mdf(data_rarified)
#Tell the function to plot the five most abundant phyla (you could pick something different... )
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c("Proteobacteria",   "Actinobacteria", "Acidobacteria", "Bacteroidetes","Planctomycetes"),cvd=TRUE)
#Extract the OTU table and color choices (again, putting things in the right format for the function)
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf
cdf_GP <- color_reassign(cdf_GP,group_assignment = c("Proteobacteria",   "Actinobacteria", "Acidobacteria", "Bacteroidetes","Planctomycetes"),color_assignment = c("micro_cvd_green", "micro_cvd_orange","micro_cvd_blue","micro_cvd_purple","micro_cvd_turquoise"))
#Define the legend
GP_legend <-custom_legend(mdf_GP, cdf_GP)

new_cdf_GP <- color_reassign(cdf_GP,
                             group_assignment = c("Bacteroidetes", "Planctomycetes","Acidobacteria","Actinobacteria","Proteobacteria"),
                             color_assignment = c("micro_purple", "micro_cvd_turquoise","micro_blue","micro_cvd_orange","micro_cvd_green"))


new_groups <- extend_group(mdf_GP, new_cdf_GP, "Phylum", "Genus", "Proteobacteria", existing_palette = "micro_green", new_palette = "micro_green", n_add = 5)
new_groups2 <- extend_group(new_groups$mdf, new_groups$cdf, "Phylum", "Genus", "Actinobacteria", existing_palette = "micro_cvd_orange", new_palette = "micro_orange", n_add = 5)


#Define your new legend with the expanded groups (FUNCTIONS TO GENERATE THIS CUSTOM LEGEND ARE AT THE BOTTOM OF THIS R FILE)
GP_legend_new <-custom_legend2column(new_groups2$mdf, new_groups2$cdf,legend_key_size=0.5,legend_text_size = 10)


#Define the plot that you will be making
plot <- plot_microshades(new_groups2$mdf, new_groups2$cdf)
#Define all of the formatting aspects of the plot
plot_diff <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 5)) +
  facet_grid(~fct_relevel(Clade,'HNG','Southern Appalachian','North Appalachian','BRE','Crevice'), scale="free_x", space = "free_x") +
  theme(axis.text.x = element_text(size= 5)) +
  theme(plot.margin = margin(6,20,6,6))+ylab("Relative Abundance")
#Make the plot

#Make the plot
plot_grid(plot_diff, GP_legend_new,  rel_widths = c(1, .2))


#pdf(file = "Taxa_Barplot_redo.pdf",   
#    width = 14, # The width of the plot in inches
#    height = 7) # The height of the plot in inches

#plot_grid(plot_diff, GP_legend,  rel_widths = c(1, .1))
dev.off()


######Subset Data #####

####This will divide your phyloseq object into 2, one with all salamander samples and one with environment
Salamanders<-subset_samples(data_rarified, Sample_Type == "Salamander")
Environment<-subset_samples(data_rarified, Sample_Type != "Salamander")
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the most abundant phyla on salamanders and crevices

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summed_phyla_salamanders<-rowMeans(data.frame(otu_table(tax_glom(Salamanders,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_salamanders<-data.frame(tax_table(tax_glom(Salamanders,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list_salamanders<-phyla_list_salamanders[order(-summed_phyla_salamanders)][1:5]
sorted_phyla_percentages_salamanders<-summed_phyla_salamanders[order(-summed_phyla_salamanders)][1:5]*100/20000
data.frame(sorted_phyla_list_salamanders,sorted_phyla_percentages_salamanders)

summed_phyla_crevices<-rowMeans(data.frame(otu_table(tax_glom(Crevicesdata,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list_crevices<-phyla_list_crevices[order(-summed_phyla_crevices)][1:5]
sorted_phyla_percentages_crevices<-summed_phyla_crevices[order(-summed_phyla_crevices)][1:5]*100/20000
data.frame(sorted_phyla_list_crevices,sorted_phyla_percentages_crevices)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the most abundant genera on salamanders and crevices

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summed_genus_salamanders<-rowMeans(data.frame(otu_table(tax_glom(Salamanders,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_salamanders<-data.frame(tax_table(tax_glom(Salamanders,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_genus_list_salamanders<-genus_list_salamanders[order(-summed_genus_salamanders)][1:5]
sorted_genus_percentages_salamanders<-summed_genus_salamanders[order(-summed_genus_salamanders)][1:5]*100/20000
data.frame(sorted_genus_list_salamanders,sorted_genus_percentages_salamanders)

summed_genus_crevices<-rowMeans(data.frame(otu_table(tax_glom(Crevicesdata,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Genus")))[,6]
family_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Genus")))[,5]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_genus_list_crevices<-genus_list_crevices[order(-summed_genus_crevices)][1:5]
sorted_family_list_crevices<-family_list_crevices[order(-summed_genus_crevices)][1:5]
sorted_genus_percentages_crevices<-summed_genus_crevices[order(-summed_genus_crevices)][1:5]*100/20000
data.frame(sorted_genus_list_crevices,sorted_genus_percentages_crevices)

