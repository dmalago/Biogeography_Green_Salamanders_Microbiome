
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

#Pull out metadata
SD<-sample_data(Salamanders)

#Calculate Diversity metrics
faithpd<-pd(t(otu_table(Salamanders)), phy_tree(Salamanders), include.root=FALSE)
shannon<-diversity(t(otu_table(Salamanders)),index='shannon')
richness<-colSums(sign(otu_table(Salamanders)))

#Make a dataframe with diversity and metadata
alpha_df<-data.frame(SD$Longitude,SD$Latitude,SD$Clade,SD$Subclade,SD$Study.Site,richness,faithpd[,1],shannon)
colnames(alpha_df)<-c('long','lat','clade','subclade','site','rich','faith','shannon')

#Reorder the dataframe according to clade
alpha_df<-alpha_df[order(alpha_df$clade),]

##################################### Test for spatial autocorrelation using Moran's I############################################
set.seed(1)
newlat<-c()
newlong<-c()
for (k in 1:length(alpha_df$lat)){
  newlat<-c(newlat,alpha_df$lat[k]+sample(1E2,1)/2E4)
  newlong<-c(newlong,alpha_df$long[k]+sample(1E2,1)/2E4)
}

alpha_df<-data.frame(alpha_df,newlat,newlong)
colnames(alpha_df)<-c('long','lat','clade','subclade','site','rich','faith','shannon','newlat','newlong')

latlong<-as.matrix(data.frame(newlong,newlat))

nb <- chooseCN(coordinates(latlong), type = 5, d1 = 0, d2 = 10000, plot.nb = FALSE)
distnb <- nbdists(nb, latlong)
fdist <- lapply(distnb, function(x) 1 / x^1)
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)  

set.seed(1)
moran.mc(alpha_df$rich, lw,alternative='two.sided',nsim=10000)   #spatial autocorrelation
moranNP.randtest(alpha_df$rich, lw,alter='two-sided')           #positive spatial autocorrelation

moran.mc(alpha_df$shannon, lw,alternative='two.sided',nsim=10000)   #spatial autocorrelation
moranNP.randtest(alpha_df$shannon, lw,alter='two-sided')           #positive spatial autocorrelation

moran.mc(alpha_df$faith, lw,nsim=10000,alternative='two.sided') #spatial autocorrelation in Faith's PD
moranNP.randtest(alpha_df$faith, lw,alter = 'two-sided')        

par(mar=c(5,5,5,5))
lag<-lag.listw(lw, alpha_df$faith)
plot(lag~alpha_df$faith, pch=16,xlab='Faith\'s PD',ylab='Spatially Lagged Faith\'s PD',cex.lab=1.5)
M1 <- lm(lag ~ alpha_df$fait)
abline(M1, col="blue")


within_dist<-c()

for (k in 1:102){
  listme<-which(alpha_df$site==alpha_df$site[k])
  listme<-listme[listme!=k]
  for (j in 1:length(listme)){
    if (k<listme[j]){
      listme[j]<-listme[j]-1
    }
  }
  if (max(distnb[[k]][listme])>0.01){print(k)}
  within_dist<-c(within_dist,distnb[[k]][listme])
}
max(within_dist)


nb <- chooseCN(coordinates(latlong), type = 5, d1 = 1.1*max(within_dist), d2 = 10000, plot.nb = FALSE)
distnb <- nbdists(nb, latlong)
fdist <- lapply(distnb, function(x) 1 / x^1)
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)  

set.seed(1)
moran.mc(alpha_df$rich, lw,alternative='two.sided',nsim=10000)   #spatial autocorrelation
moranNP.randtest(alpha_df$rich, lw,alter='two-sided')           #positive spatial autocorrelation

moran.mc(alpha_df$shannon, lw,alternative='two.sided',nsim=10000)   #spatial autocorrelation
moranNP.randtest(alpha_df$shannon, lw,alter='two-sided')           #positive spatial autocorrelation

moran.mc(alpha_df$faith, lw,nsim=10000,alternative='two.sided') #spatial autocorrelation in Faith's PD
moranNP.randtest(alpha_df$faith, lw,alter = 'two-sided')        


