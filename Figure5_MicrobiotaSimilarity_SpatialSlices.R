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


set.seed(1)
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

lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitude')

#Adding a tiny amount of jitter around salamanders found at same site
j = 1E-5
num_rows <- nrow(lonlat)  # Get the number of rows directly

lonlat2 <- lonlat %>%
  mutate(
    longitude = longitude + runif(n = num_rows, min = -j, max = j),
    latitude = latitude + runif(n = num_rows, min = -j, max = j)
  )

dist_matrix<-geodist(lonlat2,measure="geodesic")
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)

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

env_matrix<-as.matrix(vegdist(envmat+10,metric='euclidean',diag=0,upper=TRUE))



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

ddists<-matrix(0,length(sitelist),length(sitelist))
rownames(ddists)<-sitelist
colnames(ddists)<-sitelist

cladal<-c()
for (k in 1:length(sitelist)){
  hit1<-which(sample_data(Salamanders)$Study.Site==sitelist[k])
  cladal<-c(cladal,sample_data(Salamanders)$Clade[hit1[1]])
  for (j in 1:length(sitelist)){
    hit2<-which(sample_data(Salamanders)$Study.Site==sitelist[j])
    
    if (k==j){
      p1<-which(rownames(adists)==sitelist[k])
      adists[p1,p1]<-0#mean(edist[hit1,hit2])
      ddists[p1,p1]<-0
    }
    else{
      p1<-which(rownames(adists)==sitelist[k])
      p2<-which(rownames(adists)==sitelist[j])
      adists[p1,p2]<-mean(edist[hit1,hit2])
      adists[p2,p1]<-mean(edist[hit1,hit2])
      ddists[p1,p2]<-dist_matrix[hit1[1],hit2[1]]
      ddists[p2,p1]<-dist_matrix[hit1[1],hit2[1]]
    }
  }
  
}

cladalcount<-cladal
cladalcount[cladalcount=='HNG']<-1
cladalcount[cladalcount=='BRE']<-2
cladalcount[cladalcount=='North Appalachian']<-3
cladalcount[cladalcount=='Southern Appalachian']<-4
cladalcount<-as.numeric(cladalcount)

binner<-list()
fer<-c()
for (m in 1:20){
  tlists<-c()
for (k in 1:length(cladal)){
  for (j in 1:length(cladal)){
    low<-(m-1)*30000
    high<-m*30000
    if (ddists[k,j]>low && ddists[k,j]<high){
      if (cladalcount[k]<cladalcount[j]){
      tlists<-c(tlists,cladalcount[k]+10*cladalcount[j])
      }else{
        tlists<-c(tlists,cladalcount[j]+10*cladalcount[k])
      }
    }
  }
}
  fer<-c(fer,length(tlists))
  binner<-list.append(binner,tlists)
}


aa<-c()
bb<-c()
cc<-c()
ff<-c()
ab<-c()
ac<-c()
ad<-c()
bc<-c()
bd<-c()
cd<-c()
for (k in 1:length(binner)){
  aa<-c(aa,length(which(binner[[k]]==11)))
  bb<-c(bb,length(which(binner[[k]]==22)))
  cc<-c(cc,length(which(binner[[k]]==33)))
  ff<-c(ff,length(which(binner[[k]]==44)))
  ab<-c(ab,length(which(binner[[k]]==21)))
  ac<-c(ac,length(which(binner[[k]]==31)))
  ad<-c(ad,length(which(binner[[k]]==41)))
  bc<-c(bc,length(which(binner[[k]]==32)))
  bd<-c(bd,length(which(binner[[k]]==42)))
  cd<-c(cd,length(which(binner[[k]]==43)))
  
}

mm<-rbind(aa,bb,cc,ff,ab,ac,ad,bc,bd,cd)
m1<-c()
m2<-c()
m3<-c()
for (k in 1:20){
m1<-c(m1,rep((k)*30,10))
m2<-c(m2,c('HNG-HNG','BRE-BRE','NA-NA','SA-SA','HNG-BRE','HNG-NA','HNG-SA','BRE-NA','BRE-SA','NA-SA'))
m3<-c(m3,t(mm)[k,])
}

mdf<-data.frame(m1,m2,m3)
colnames(mdf)<-c('distance','connect','count')
mdf$connect<-factor(mdf$connect,levels=c('HNG-HNG','BRE-BRE','NA-NA','SA-SA','HNG-BRE','HNG-NA','HNG-SA','BRE-NA','BRE-SA','NA-SA'))

dbar<-ggplot(mdf, aes(fill=connect, y=count, x=distance))
dbar<-dbar+ ylab('Number of Site Pairs')+xlab('Distance (km)')+theme(axis.title = element_text(size = 18))
dbar<-dbar+theme(axis.text = element_text(size = 12))
dbar<-dbar+theme(panel.border = element_rect(color = 'black',fill = NA,size = 1))
dbar<-dbar+geom_col_pattern(aes(fill=connect,pattern_fill=connect,width=20),color='black',pattern='stripe', pattern_angle=60, pattern_spacing=0.025,pattern_density=0.5,width=1,pattern_linetype = 0)
dbar<-dbar+scale_fill_manual(values=c("green", "blue","yellow",'purple','green','green','green','blue','blue','yellow'))
dbar<-dbar+scale_pattern_fill_manual(values=c("green", "blue","yellow",'purple','blue','yellow','purple','yellow','purple','purple'))
dbar<-dbar+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dbar<-dbar+guides(fill=guide_legend("Site Pairs"),pattern_fill=guide_legend("Site Pairs"))
plot(dbar)                                                                                                                                                                                                                                                                                                                                                                       


common_theme <- theme_classic() +
  theme(text = element_text(size=20),  # Sets a good size for readability
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = "none",  # Hides the legend
        plot.title = element_text(size=0, hjust = 0.5),
        axis.title = element_text(size=30),
        axis.text = element_text(size=22))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                                     NA-NA vs NA-HNG

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin1<-5
bin2<-6
which_sites<-c()
which_clades<-c()
site1<-c()
site2<-c()
clade1<-c()
clade2<-c()
for (k in 1:length(cladal)){
  hitter<-k
  for (j in 1:hitter){
    low<-bin1*30000
    high<-bin2*30000
    if (ddists[k,j]>low && ddists[k,j]<high){
      if (cladalcount[k]<cladalcount[j]){
        which_sites<-c(which_sites,paste0(rownames(ddists)[k],' & ',rownames(ddists)[j]))
        which_clades<-c(which_clades,paste0(cladal[k],' & ',cladal[j]))
        site1<-c(site1,rownames(ddists)[k])
        site2<-c(site2,rownames(ddists)[j])
        clade1<-c(clade1,cladal[k])
        clade2<-c(clade2,cladal[j])
      }else{
        which_sites<-c(which_sites,paste0(rownames(ddists)[j],' & ',rownames(ddists)[k]))
        which_clades<-c(which_clades,paste0(cladal[j],' & ',cladal[k]))
        site1<-c(site1,rownames(ddists)[j])
        site2<-c(site2,rownames(ddists)[k])
        clade1<-c(clade1,cladal[j])
        clade2<-c(clade2,cladal[k])
      }
    }
  }
}

pair11<-which(which_clades=='North Appalachian & North Appalachian')

pair21<-which(which_clades=='HNG & North Appalachian')
pair22<-which(which_clades=='North Appalachian & HNG')

sharedpairs<-intersect(unique(c(site1[pair11],site2[pair11])),unique(site2[pair21]))

NAgood<-sharedpairs
HNGgood<-unique(site1[pair11][which(site2[pair11] %in% sharedpairs)])
BREgood<-unique(site1[pair21][which(site2[pair21] %in% sharedpairs)])

sharedpairs<-intersect(unique(c(site1[pair11],site2[pair11])),unique(site2[pair21]))

NAgood<-unique(c(site2[pair11][which(site1[pair11] %in% sharedpairs)],site1[pair11][which(site2[pair11] %in% sharedpairs)]))
HNGgood<-unique(site1[pair21][which(site2[pair21] %in% sharedpairs)])

NANAcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% NAgood)]
NAHNGcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]

f1<-rowMeans(NANAcomparison)
f2<-rowMeans(NAHNGcomparison)

testme<-wilcox.test(f1,f2,paired=TRUE)
median(f1)
median(f2)

envNANAcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% NAgood)]
envNAHNGcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]

envf1<-rowMeans(envNANAcomparison)
envf2<-rowMeans(envNAHNGcomparison)

testme<-wilcox.test(envf1,envf2,paired=TRUE)
median(envf1)
median(envf2)

dNANAcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% NAgood)]
dNAHNGcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% sharedpairs),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]

df1<-rowMeans(dNANAcomparison)
df2<-rowMeans(dNAHNGcomparison)
dtestme<-wilcox.test(df1,df2,paired=TRUE)
median(df1)
median(df2)


Pdf<-data.frame(c(rep('NA-NA',length(f1)),rep('HNG-NA',length(f2))),c(f1,f2),c(envf1,envf2),c(df1/1000,df2/1000))
colnames(Pdf)<-c('types','cdist','edist','ddist')
#Switch the ordering of the types (to the order you want them presented in your plots)
Pdf$types<-factor(Pdf$types,levels=c('NA-NA','HNG-NA'))


##########################            Violin Plots          #########################################################################################################################################

#Define violin plot
plotDistance <- ggplot(Pdf, aes(x=types, y=ddist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Spatial Distance (km)')+xlab('')
#Choose the size of font for the axes titles and labels
plotDistance<-plotDistance+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotDistance<-plotDistance+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotDistance<-plotDistance+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('yellow', 'green'))+scale_pattern_fill2_manual(values = c('yellow', 'green'))
#Add boxplots inside the violins
plotDistance<-plotDistance+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-dtestme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('***')
plotDistance<-plotDistance+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=170,size=5)
plotDistance<-plotDistance+common_theme
plotDistance


#Define violin plot
plotEnvironment <- ggplot(Pdf, aes(x=types, y=edist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Environmental Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotEnvironment<-plotEnvironment+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotEnvironment<-plotEnvironment+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotEnvironment<-plotEnvironment+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('yellow', 'green'))+scale_pattern_fill2_manual(values = c('yellow', 'green'))
#Add boxplots inside the violins
plotEnvironment<-plotEnvironment+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-envtestme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('*')
plotEnvironment<-plotEnvironment+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.0875,size=5)
plotEnvironment<-plotEnvironment+common_theme
plotEnvironment


#Define violin plot
plotUnifrac <- ggplot(Pdf, aes(x=types, y=cdist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('UniFrac Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotUnifrac<-plotUnifrac+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotUnifrac<-plotUnifrac+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotUnifrac<-plotUnifrac+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('yellow', 'green'))+scale_pattern_fill2_manual(values = c('yellow', 'green'))
#Add boxplots inside the violins
plotUnifrac<-plotUnifrac+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-testme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('***')
plotUnifrac<-plotUnifrac+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.95,size=5)
plotUnifrac<-plotUnifrac+common_theme
plotUnifrac






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                                     NA-BRE vs NA-HNG (Shorter Distance Bin)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin1<-6
bin2<-7
which_sites<-c()
which_clades<-c()
site1<-c()
site2<-c()
clade1<-c()
clade2<-c()
for (k in 1:length(cladal)){
  hitter<-k
  for (j in 1:hitter){
    low<-bin1*30000
    high<-bin2*30000
    if (ddists[k,j]>low && ddists[k,j]<high){
      if (cladalcount[k]<cladalcount[j]){
      which_sites<-c(which_sites,paste0(rownames(ddists)[k],' & ',rownames(ddists)[j]))
      which_clades<-c(which_clades,paste0(cladal[k],' & ',cladal[j]))
      site1<-c(site1,rownames(ddists)[k])
      site2<-c(site2,rownames(ddists)[j])
      clade1<-c(clade1,cladal[k])
      clade2<-c(clade2,cladal[j])
      }else{
        which_sites<-c(which_sites,paste0(rownames(ddists)[j],' & ',rownames(ddists)[k]))
        which_clades<-c(which_clades,paste0(cladal[j],' & ',cladal[k]))
        site1<-c(site1,rownames(ddists)[j])
        site2<-c(site2,rownames(ddists)[k])
        clade1<-c(clade1,cladal[j])
        clade2<-c(clade2,cladal[k])
      }
    }
  }
}

pair11<-which(which_clades=='HNG & North Appalachian')
pair12<-which(which_clades=='North Appalachian & HNG')

pair21<-which(which_clades=='BRE & North Appalachian')
pair22<-which(which_clades=='North Appalachian & BRE')

sharedpairs<-intersect(unique(site2[pair11]),unique(site2[pair21]))

NAgood<-sharedpairs
HNGgood<-unique(site1[pair11][which(site2[pair11] %in% sharedpairs)])
BREgood<-unique(site1[pair21][which(site2[pair21] %in% sharedpairs)])

uni<-as.matrix(UniFrac(Salamanders))

NAHNGcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
NABREcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

f1<-rowMeans(NAHNGcomparison)
f2<-rowMeans(NABREcomparison)

testme<-wilcox.test(f1,f2,paired=TRUE)
median(f1)
median(f2)



envNAHNGcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
envNABREcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

envf1<-rowMeans(envNAHNGcomparison)
envf2<-rowMeans(envNABREcomparison)
envtestme<-wilcox.test(envf1,envf2,paired=TRUE)
median(envf1)
median(envf2)


dNAHNGcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
dNABREcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

df1<-rowMeans(dNAHNGcomparison)
df2<-rowMeans(dNABREcomparison)
dtestme<-wilcox.test(df1,df2,paired=TRUE)
median(df1)
median(df2)


Pdf<-data.frame(c(rep('HNG-NA',length(f1)),rep('BRE-NA',length(f2))),c(f1,f2),c(envf1,envf2),c(df1/1000,df2/1000))
colnames(Pdf)<-c('types','cdist','edist','ddist')
#Switch the ordering of the types (to the order you want them presented in your plots)
Pdf$types<-factor(Pdf$types,levels=c('BRE-NA','HNG-NA'))


##########################            Violin Plots          #########################################################################################################################################

#Define violin plot
plotDistance <- ggplot(Pdf, aes(x=types, y=ddist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Spatial Distance (km)')+xlab('')
#Choose the size of font for the axes titles and labels
plotDistance<-plotDistance+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotDistance<-plotDistance+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotDistance<-plotDistance+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotDistance<-plotDistance+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
plotDistance<-plotDistance+common_theme
plotDistance


#Define violin plot
plotEnvironment <- ggplot(Pdf, aes(x=types, y=edist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Environmental Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotEnvironment<-plotEnvironment+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotEnvironment<-plotEnvironment+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotEnvironment<-plotEnvironment+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotEnvironment<-plotEnvironment+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-envtestme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('**')
plotEnvironment<-plotEnvironment+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.078,size=5)
plotEnvironment<-plotEnvironment+common_theme
plotEnvironment


#Define violin plot
plotUnifrac <- ggplot(Pdf, aes(x=types, y=cdist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('UniFrac Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotUnifrac<-plotUnifrac+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotUnifrac<-plotUnifrac+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotUnifrac<-plotUnifrac+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotUnifrac<-plotUnifrac+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-testme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('**')
plotUnifrac<-plotUnifrac+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.76,size=5)
plotUnifrac<-plotUnifrac+common_theme
plotUnifrac




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                                     NA-BRE vs NA-HNG  (Larger Distance Bins)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bin1<-7
bin2<-8
which_sites<-c()
which_clades<-c()
site1<-c()
site2<-c()
clade1<-c()
clade2<-c()
for (k in 1:length(cladal)){
  hitter<-k
  for (j in 1:hitter){
    low<-bin1*30000
    high<-bin2*30000
    if (ddists[k,j]>low && ddists[k,j]<high){
      if (cladalcount[k]<cladalcount[j]){
        which_sites<-c(which_sites,paste0(rownames(ddists)[k],' & ',rownames(ddists)[j]))
        which_clades<-c(which_clades,paste0(cladal[k],' & ',cladal[j]))
        site1<-c(site1,rownames(ddists)[k])
        site2<-c(site2,rownames(ddists)[j])
        clade1<-c(clade1,cladal[k])
        clade2<-c(clade2,cladal[j])
      }else{
        which_sites<-c(which_sites,paste0(rownames(ddists)[j],' & ',rownames(ddists)[k]))
        which_clades<-c(which_clades,paste0(cladal[j],' & ',cladal[k]))
        site1<-c(site1,rownames(ddists)[j])
        site2<-c(site2,rownames(ddists)[k])
        clade1<-c(clade1,cladal[j])
        clade2<-c(clade2,cladal[k])
      }
    }
  }
}

pair11<-which(which_clades=='HNG & North Appalachian')
pair12<-which(which_clades=='North Appalachian & HNG')

pair21<-which(which_clades=='BRE & North Appalachian')
pair22<-which(which_clades=='North Appalachian & BRE')

sharedpairs<-intersect(unique(site2[pair11]),unique(site2[pair21]))

NAgood<-sharedpairs
HNGgood<-unique(site1[pair11][which(site2[pair11] %in% sharedpairs)])
BREgood<-unique(site1[pair21][which(site2[pair21] %in% sharedpairs)])

uni<-as.matrix(UniFrac(Salamanders))

NAHNGcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
NABREcomparison<-uni[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

f1<-rowMeans(NAHNGcomparison)
f2<-rowMeans(NABREcomparison)

testme<-wilcox.test(f1,f2,paired=TRUE)
median(f1)
median(f2)



envNAHNGcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
envNABREcomparison<-env_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

envf1<-rowMeans(envNAHNGcomparison)
envf2<-rowMeans(envNABREcomparison)
envtestme<-wilcox.test(envf1,envf2,paired=TRUE)
median(envf1)
median(envf2)


dNAHNGcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% HNGgood)]
dNABREcomparison<-dist_matrix[which(sample_data(Salamanders)$Study.Site %in% NAgood),which(sample_data(Salamanders)$Study.Site %in% BREgood)]

df1<-rowMeans(dNAHNGcomparison)
df2<-rowMeans(dNABREcomparison)
dtestme<-wilcox.test(df1,df2,paired=TRUE)
median(df1)
median(df2)


Pdf<-data.frame(c(rep('HNG-NA',length(f1)),rep('BRE-NA',length(f2))),c(f1,f2),c(envf1,envf2),c(df1/1000,df2/1000))
colnames(Pdf)<-c('types','cdist','edist','ddist')
#Switch the ordering of the types (to the order you want them presented in your plots)
Pdf$types<-factor(Pdf$types,levels=c('BRE-NA','HNG-NA'))


##########################            Violin Plots          #########################################################################################################################################


#Define violin plot
plotDistance <- ggplot(Pdf, aes(x=types, y=ddist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Spatial Distance (km)')+xlab('')
#Choose the size of font for the axes titles and labels
plotDistance<-plotDistance+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotDistance<-plotDistance+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotDistance<-plotDistance+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotDistance<-plotDistance+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-dtestme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('**')
plotDistance<-plotDistance+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=216.5,size=5)
plotDistance<-plotDistance+common_theme
plotDistance

#Define violin plot
plotEnvironment <- ggplot(Pdf, aes(x=types, y=edist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('Environmental Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotEnvironment<-plotEnvironment+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotEnvironment<-plotEnvironment+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotEnvironment<-plotEnvironment+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotEnvironment<-plotEnvironment+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-envtestme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('*')
plotEnvironment<-plotEnvironment+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.078,size=5)
plotEnvironment<-plotEnvironment+common_theme
plotEnvironment



#Define violin plot
plotUnifrac <- ggplot(Pdf, aes(x=types, y=cdist)) +geom_violin_pattern(aes(fill=types,pattern_fill=types,pattern_fill2=types,pattern=types),pattern_spacing=0.05,pattern_density=0.5,width=1) + theme(axis.title.x = element_blank())+ylab('UniFrac Distance')+xlab('')
#Choose the size of font for the axes titles and labels
plotUnifrac<-plotUnifrac+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
plotUnifrac<-plotUnifrac+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
plotUnifrac<-plotUnifrac+scale_pattern_manual(values=c('stripe','stripe'))+scale_fill_manual(values=c('yellow', 'yellow'))+scale_pattern_fill_manual(values = c('blue', 'green'))+scale_pattern_fill2_manual(values = c('blue', 'green'))
#Add boxplots inside the violins
plotUnifrac<-plotUnifrac+geom_boxplot(aes(fill=types),width=0.1) + theme(axis.title.x = element_blank())
group1<-1
group2<-2
p.adj<-testme$p.value
stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
pstar<-c('*')
plotUnifrac<-plotUnifrac+stat_pvalue_manual(stat.test_richness,label="{pstar}",y.position=0.76,size=5)
plotUnifrac<-plotUnifrac+common_theme
plotUnifrac







