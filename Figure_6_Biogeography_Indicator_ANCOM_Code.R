#############Biogeography Greens Indicator Species and ANCOM Code ###########

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


######Indicator for Clades#####
set.seed(1)

library('indicspecies')
library('data.table')

#type_clade<-sample_data(data_rarified)$Clade
indval_clade = multipatt(t(otu_table(data_rarified)), type_clade, control = how(nperm=5000), duleg=TRUE)
summary(indval_clade,indvalcomp=TRUE)
#saveRDS(indval_clade, file = "indval_clade_final.rds")
indval_clade<-readRDS(file = "indval_clade_final.rds")
indisp.sign<-as.data.table(indval_clade$sign, keep.rownames=TRUE)
abundant_taxa<-rownames(otu_table(data_rarified))[which(apply(otu_table(data_rarified),MARGIN=1,FUN=max)/colSums(otu_table(data_rarified))[1]>0.01)]

indisp.sign<-indisp.sign[which(indisp.sign$rn %in% abundant_taxa)]
#add adjusted p-value
indisp.sign[ ,p.value.bh := p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
df_clade<-data.frame(indisp.sign[p.value.bh<=0.05, ])
temptax <- as.data.frame(tax_table(data_rarified)) #get your taxa
library(tibble)
temptax <- tibble::rownames_to_column(temptax, "rn")

mytax <- rbind(temptax[temptax$rn == "seq1069", ], 
               temptax[temptax$rn == "seq706", ], 
               temptax[temptax$rn == "seq360", ], 
               temptax[temptax$rn == "seq1433", ], 
               temptax[temptax$rn == "seq850", ], 
               temptax[temptax$rn == "seq1735", ], 
               temptax[temptax$rn == "seq962", ], 
               temptax[temptax$rn == "seq552", ], 
               temptax[temptax$rn == "seq689", ], 
               temptax[temptax$rn == "seq1725", ],
               temptax[temptax$rn == "seq265", ],
               temptax[temptax$rn == "seq472", ],
               temptax[temptax$rn == "seq387", ],
               temptax[temptax$rn == "seq300", ],
               temptax[temptax$rn == "seq1095", ],
               temptax[temptax$rn == "seq1217", ],
               temptax[temptax$rn == "seq1443", ],
               temptax[temptax$rn == "seq176", ],
               temptax[temptax$rn == "seq464", ],
               temptax[temptax$rn == "seq1667", ])

taxa_info<-mytax

colnames(taxa_info)[1] <- "rn"




df_clade<-inner_join(mytax,df_clade, by = "rn")




#write.csv(df_clade,'Clade_ASV_indicators.csv') 
#I added a clade column based on the index.

d<-read.csv('Clade_ASV_indicators.csv')


pdf(file = "Clade_Indicator.pdf", width = 15, height = 15)

ggplot(d, aes(Clade, Genus_No_NA)) +
  geom_point(aes(color = p.value.bh, size = stat)) +
  scale_color_gradient(low = "red", high = "gold") +
  scale_x_discrete(position = "top") +
  scale_size(range = c(2, 10)) +
  theme_bw() +theme(text = element_text(size=20))+
  theme(legend.position = "right")
dev.off()





#####ANCOM for Clades####


#ANCOMBC2 will not run with ASVs as is,it requires something taxonomic, so I am going to change the species names of all my ASVs to the actual ASVs.

df<-c()
df<-as.data.frame(tax_table(data_rarified))
df$Species<-c()

library(tibble)
df <- tibble::rownames_to_column(df, "Species")
rownames(df)<-df$Species
df2 <- df[, c("Kingdom", "Phylum" , "Class",   "Order",   "Family",  "Genus", "Species" )]
data_rarified@tax_table<-tax_table(as.matrix(df2))
rownames(tax_table(data_rarified))<-df2$Species

unload("phyloseq")
library(ANCOMBC)


library(dplyr)

# Get all unique clade pairs
clade_pairs <- combn(unique(sample_data(data_rarified)$Clade), 2, simplify = FALSE)

# Run ANCOMBC2 for each pair
pairwise_results <- lapply(clade_pairs, function(pair) {
  # Subset phyloseq object for the current pair of clades
  subset_data <- prune_samples(
    sample_data(data_rarified)$Clade %in% pair, 
    data_rarified
  )
  
  # Ensure levels are dropped for the subset
  sample_data(subset_data)$Clade <- factor(
    sample_data(subset_data)$Clade,
    levels = pair
  )
  
  # Run ANCOMBC2 on the subset data
  ancombc2(
    subset_data, group = "Clade", p_adj_method = "holm",
    fix_formula = "Clade", tax_level = "Species",
    verbose = TRUE, pairwise = FALSE, struc_zero = FALSE, pseudo_sens = TRUE
  )
})

# Optionally, name the results for clarity
names(pairwise_results) <- sapply(clade_pairs, paste, collapse = " vs ")

head(pairwise_results)

#saveRDS(pairwise_results,file = "Clade_ANCOM.rds")

pairwise_results<-readRDS(file = "Clade_ANCOM.rds")

str(pairwise_results)

# Load necessary library
library(phyloseq)

# Extract the taxonomy table and ASV IDs from the phyloseq object
tax_table <- as.data.frame(tax_table(data_rarified))
asv_ids <- rownames(otu_table(data_rarified))

# Ensure the row names of the taxonomy table are the ASV IDs
rownames(tax_table) <- asv_ids

# Extract genus-level information; check your taxonomy table for the correct column name
# Use gsub to remove the 'g__' prefix from the genus names
genus_names <- gsub("^g__", "", tax_table$Genus)

# Create the mapping from ASV IDs to genus names
genus_mapping <- setNames(genus_names, asv_ids)

# Handle missing or ambiguous assignments: Replace NA or empty genus names
# You can adjust this to use other taxonomic levels if genus is not resolved
genus_mapping[is.na(genus_mapping) | genus_mapping == ""] <- "Unclassified"

# Optionally, fill in unresolved genus names with the highest available resolved taxonomic rank
unresolved <- is.na(genus_mapping) | genus_mapping == "Unclassified"
genus_mapping[unresolved] <- apply(tax_table[unresolved, , drop = FALSE], 1, function(x) {
  non_na <- na.omit(x)
  if (length(non_na) > 0) tail(non_na, 1) else "Unclassified"
})

# Print the first few entries of the mapping to verify
head(genus_mapping)





update_res_row_names_conditional <- function(pairwise_results, genus_mapping) {
  # Loop over each comparison in the pairwise_results list
  for (comparison in names(pairwise_results)) {
    # Check if 'res' exists in the list
    if ("res" %in% names(pairwise_results[[comparison]])) {
      # Get the current 'res' data frame
      res_df <- pairwise_results[[comparison]]$res
      
      # Find all 'diff_' column names
      diff_col_names <- grep("^diff_", names(res_df), value = TRUE)
      
      # Check if exactly one 'diff_' column is present or select the correct one based on additional logic
      if (length(diff_col_names) == 1) {
        diff_col_name <- diff_col_names[1]  # straightforward case: one diff column
      } else if (length(diff_col_names) > 1) {
        # More complex logic needed if multiple diff columns are present
        # Example: select based on additional criteria or simply pick the first one
        diff_col_name <- diff_col_names[1]  # or any logic to select the appropriate column
      } else {
        warning(paste("No 'diff_' column found in", comparison))
        next  # Skip this iteration if no appropriate column is found
      }
      
      # Filter to update only rows where diff_ is TRUE
      rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
      
      # Update the row names using the genus mapping for filtered rows
      new_row_names <- genus_mapping[rownames(res_df)[rows_to_update]]
      # Handle missing mappings
      is_na <- is.na(new_row_names)
      if (any(is_na)) {
        warning(paste("Missing genus mappings found in", comparison, "for some rows."))
        new_row_names[is_na] <- rownames(res_df)[rows_to_update][is_na]  # Keep original ASV IDs if missing
      }
      
      # Update only the specific row names
      rownames(res_df)[rows_to_update] <- new_row_names
      
      # Put the updated data frame back into the list
      pairwise_results[[comparison]]$res <- res_df
    }
  }
  return(pairwise_results)
}

# Apply the function to update the pairwise_results object
pairwise_results <- update_res_row_names_conditional(pairwise_results, genus_mapping)


# Print structure or a part of the updated pairwise_results to verify changes
str(pairwise_results)



# Extract each component of pairwise_results into separate variables
North_Appalachian_vs_BRE <- pairwise_results[["North Appalachian vs BRE"]]
update_res_row_names_for_North_Appalachian_vs_BRE <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' - assumed to be 'diff_CladeBRE' as per structure
  diff_col_name <- "diff_CladeBRE"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c( new_row_names[i], "ASV",rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Print all rows where 'diff_CladeBRE' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function
North_Appalachian_vs_BRE <- update_res_row_names_for_North_Appalachian_vs_BRE(North_Appalachian_vs_BRE, genus_mapping)






North_Appalachian_vs_Southern_Appalachian <- pairwise_results[["North Appalachian vs Southern Appalachian"]]
update_res_row_names_for_North_Appalachian_vs_Southern_Appalachian <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' - tailored for this specific comparison
  diff_col_name <- "diff_CladeSouthern Appalachian"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c( new_row_names[i], "ASV",rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Print all rows where 'diff_CladeSouthern Appalachian' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function to the specific comparison
North_Appalachian_vs_Southern_Appalachian <- update_res_row_names_for_North_Appalachian_vs_Southern_Appalachian(North_Appalachian_vs_Southern_Appalachian, genus_mapping)



North_Appalachian_vs_HNG <- pairwise_results[["North Appalachian vs HNG"]]
update_res_row_names_for_North_Appalachian_vs_HNG <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' specific to HNG - 'diff_CladeHNG'
  diff_col_name <- "diff_CladeHNG"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_CladeHNG' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c( new_row_names[i],"ASV", rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Print all rows where 'diff_CladeHNG' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function to the specific comparison
North_Appalachian_vs_HNG <- update_res_row_names_for_North_Appalachian_vs_HNG(North_Appalachian_vs_HNG, genus_mapping)

BRE_vs_Southern_Appalachian <- pairwise_results[["BRE vs Southern Appalachian"]]
update_res_row_names_for_BRE_vs_Southern_Appalachian <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' specific to Southern Appalachian - 'diff_CladeSouthern Appalachian'
  diff_col_name <- "diff_CladeSouthern Appalachian"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_CladeSouthern Appalachian' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c( new_row_names[i], "ASV",rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Print all rows where 'diff_CladeSouthern Appalachian' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function to the specific comparison
BRE_vs_Southern_Appalachian <- update_res_row_names_for_BRE_vs_Southern_Appalachian(BRE_vs_Southern_Appalachian, genus_mapping)


BRE_vs_HNG <- pairwise_results[["BRE vs HNG"]]
# Assuming BRE_vs_HNG is correctly extracted from pairwise_results
BRE_vs_HNG <- pairwise_results[["BRE vs HNG"]]

update_res_row_names_for_BRE_vs_HNG <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' specific to HNG - 'diff_CladeHNG'
  diff_col_name <- "diff_CladeHNG"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_CladeHNG' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c(new_row_names[i], "ASV", rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Optionally print all rows where 'diff_CladeHNG' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function to the specific comparison
BRE_vs_HNG <- update_res_row_names_for_BRE_vs_HNG(BRE_vs_HNG, genus_mapping)




# Retrieve the dataset from a list or other structure where it's stored
Southern_Appalachian_vs_HNG <- pairwise_results[["Southern Appalachian vs HNG"]]

# Define the function to update row names based on the genus mapping
update_res_row_names_for_Southern_Appalachian_vs_HNG <- function(comparison_data, genus_mapping) {
  # Access the 'res' data frame
  res_df <- comparison_data$res
  
  # Identify the column for 'diff_' specific to HNG - 'diff_CladeHNG'
  diff_col_name <- "diff_CladeHNG"
  
  if (diff_col_name %in% names(res_df)) {
    # Find rows where 'diff_CladeHNG' is TRUE
    rows_to_update <- which(res_df[[diff_col_name]] == TRUE)
    
    # Update the row names using the genus mapping for these rows
    original_row_names <- rownames(res_df)[rows_to_update]
    new_row_names <- genus_mapping[original_row_names]
    is_na <- is.na(new_row_names)
    if (any(is_na)) {
      warning("Missing genus mappings found for some rows.")
      new_row_names[is_na] <- original_row_names[is_na]  # Keep original ASV IDs if missing
    }
    
    # Construct new row names and remove any "NA" from the constructed name
    new_row_names <- sapply(seq_along(new_row_names), function(i) {
      parts <- c(new_row_names[i], "ASV", rows_to_update[i])
      name <- paste(parts[parts != "NA"], collapse = "_")  # Join parts, excluding "NA"
      return(name)
    })
    
    # Update only the specific row names
    rownames(res_df)[rows_to_update] <- new_row_names
    
    # Optionally print all rows where 'diff_CladeHNG' is TRUE with updated row names
    print(res_df[rows_to_update, ], row.names = TRUE)
  } else {
    warning("The specified 'diff_' column does not exist in the dataframe.")
  }
  
  # Update the data frame in the original dataset
  comparison_data$res <- res_df
  return(comparison_data)
}

# Apply the function to the specific comparison
Southern_Appalachian_vs_HNG <- update_res_row_names_for_Southern_Appalachian_vs_HNG(Southern_Appalachian_vs_HNG, genus_mapping)






##### ANCOMBC2 Pairwise Figure Code####




library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)  # For combining plots


#New Function



library(ggplot2)
library(dplyr)
library(tibble)

generate_heatmap <- function(data, clade_name, diff_column, lfc_column, p_column) {
  # Convert row names to a column before processing
  df_fig_res1 <- data %>%
    tibble::rownames_to_column(var = "taxon") %>%
    dplyr::filter(!!sym(diff_column) == TRUE) %>%  # Filter significant taxa
    dplyr::mutate(
      lfc = round(!!sym(lfc_column), 2)  # Round log fold change
    )
  
  df_fig_res2 <- data %>%
    tibble::rownames_to_column(var = "taxon") %>%
    dplyr::filter(!!sym(diff_column) == TRUE) %>%
    dplyr::mutate(
      color = ifelse(!!sym(p_column) < 0.05, "aquamarine3", "black")  # Determine significance color
    )
  
  # Merge the two data frames
  df_fig_res <- df_fig_res1 %>%
    dplyr::left_join(df_fig_res2, by = "taxon") %>%
    dplyr::arrange(desc(lfc))  # Sort by the magnitude of log fold change
  
  # Set factor levels for 'taxon' to enforce plotting order
  df_fig_res$taxon <- factor(df_fig_res$taxon, levels = unique(df_fig_res$taxon))
  
  # Fixed plotting limits
  fixed_lo <- -2.5  # Fixed lower limit for LFC
  fixed_up <- 2.5   # Fixed upper limit for LFC
  fixed_mid <- 0  # Fixed midpoint, usually zero for LFC
  
  # Generate the heatmap
  heatmap <- df_fig_res %>%
    ggplot(aes(x = clade_name, y = taxon, fill = lfc)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "yellow",
      na.value = "white", midpoint = fixed_mid, limit = c(fixed_lo, fixed_up),
      name = "Log Fold Change"
    ) +
    geom_text(aes(label = lfc, color = color), size = 4) +
    scale_color_identity(guide = FALSE) +
    labs(x = "Clade Comparison", y = "Taxon", title = paste(clade_name)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels if needed
    )
  
  return(heatmap)
}



# Function to remove the legend from a ggplot object
remove_legend <- function(heatmap_plot) {
  # Modify the heatmap to remove the legend
  updated_heatmap <- heatmap_plot + theme(legend.position = "none")
  return(updated_heatmap)
}




# Generate heatmaps for each pairwise comparison
fig_res1 <- generate_heatmap(North_Appalachian_vs_BRE$res, 
                             "North Appalachian vs BRE", 
                             "diff_CladeBRE", 
                             "lfc_CladeBRE", 
                             "p_CladeBRE")

fig_res2 <- generate_heatmap(North_Appalachian_vs_Southern_Appalachian$res, 
                             "North Appalachian vs Southern Appalachian", 
                             "diff_CladeSouthern Appalachian", 
                             "lfc_CladeSouthern Appalachian", 
                             "p_CladeSouthern Appalachian")


fig_res3 <- generate_heatmap(North_Appalachian_vs_HNG$res, 
                             "North Appalachian vs HNG", 
                             "diff_CladeHNG", 
                             "lfc_CladeHNG", 
                             "p_CladeHNG")

fig_res4 <- generate_heatmap(BRE_vs_Southern_Appalachian$res, 
                             "BRE vs Southern Appalachian", 
                             "diff_CladeSouthern Appalachian", 
                             "lfc_CladeSouthern Appalachian", 
                             "p_CladeSouthern Appalachian")


fig_res_BRE_HNG <- generate_heatmap(BRE_vs_HNG$res, 
                                    "BRE vs HNG", 
                                    "diff_CladeHNG", 
                                    "lfc_CladeHNG", 
                                    "p_CladeHNG")

fig_res_SA_HNG <- generate_heatmap(Southern_Appalachian_vs_HNG$res, 
                                   "Southern Appalachian vs HNG", 
                                   "diff_CladeHNG", 
                                   "lfc_CladeHNG", 
                                   "p_CladeHNG")

# Combine all plots using patchwork
final_plot <- fig_res1  + (remove_legend(fig_res3)+remove_legend(fig_res4)) / (remove_legend(fig_res_BRE_HNG)+remove_legend(fig_res_SA_HNG))
print(final_plot)



pdf(file = "Clade_ANCOM_ASV.pdf", width = 18, height = 14)
final_plot
dev.off()





