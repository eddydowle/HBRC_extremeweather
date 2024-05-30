#bring in files
library('tidyverse')
library('readxl')
library('devtools')
library('tidyverse')
library('viridis')
library('indicspecies') 
library('vegan')
library('glue')
library('reshape2')
library('phyloseq')
library('iNEXT.3D') 
library('ggplot2') 
library('microbiome')


setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
records<-read.csv('HBRC_records_Eddy.csv')
head(records)
otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)

#split by primer set
otu_table_splitout <- split(otu_table, otu_table$PrimerSet)
#so there is 20 different primer sets
#no idea what each of them is
otu_table_splitout


#just work through one primer set for now and then come back 

#just choosing a random not super diverse but present in most samples to play with
head(otu_table_splitout$MZ)

out_table_MZ<- otu_table_splitout$MZ %>% select(-PrimerSet) %>% group_by(UID, TaxID)

#switching to CI which is insect
out_table_MZ<- otu_table_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)

#count unique samples (no idea why dplyr thinks it grouped)
out_table_MZ %>% ungroup() %>% select(UID) %>% n_distinct()

#there is a TaxID that is NA sometimes in some assays - I have no idea what that means:
out_table_MZ[!complete.cases(out_table_MZ), ]

#sample count per dataset
#BE - 1710 - 18S general eukaryote
#BU - 1723 - 18S general eukaryote
#BX - 723 - 18S general eukaryote
#CI - 1718 - COI insects, inverts, some fish
#DG - 652
#EA - 809
#GV - 728 - vascular plants ITS1
#LG - 987 - Mifish
#LV - 642
#MZ - 1708 - Vascular plants rbcl
#RV - 1684 - 12S vertebrate ecoprimers 
#TP - 1717 - vascular plants
#UM - 1723 -microbe
#WG - 727 0 venerid clams 16S
#WV - 1713 - 16S vertebrate primers (picks up several worms too)
#XG - 861
#YG - 298
#ZC - 642
#ZP - 586
#ZV - 80

#bring in meta data
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)


#samples present and select from meta data
meta_data_clean <- meta_data %>% filter(UID %in% out_table_MZ$UID)

#what is the stuff missing (e.g. all pre or all post)

meta_data %>% filter(!UID %in% out_table_MZ$UID)

#mix pre and post from various sites no idea - maybe failed samples?

#long format OTU table
#rows samples, columns tax, values counts
library(reshape2)

out_table_MZ_wide<-dcast(out_table_MZ, TaxID~UID,value.var = "Count", fun.aggregate = sum)


#so 322 unique taxid's
#have just merged on taxid here (with aggregate sum)

#generate taxa table from TaxIDs
#from records file
taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()


length(unique(out_table_MZ$TaxID))
#321 but one is a NA

#might need to filter by taxid's in otu table
taxa_table_clean <- taxa_table %>% filter(TaxID %in% out_table_MZ$TaxID)


#generate meta file
#critical one
#Location_date
#Location_precyclone , Location_postcyclone
head(meta_data)
meta_data<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))

meta_data_clean <- meta_data %>% filter(UID %in% out_table_MZ$UID)

#okay so my three correctly shaped files
#out_table_MZ_wide
#taxa_table_clean
#meta_data_clean

#going to have to remove that weird taxid == 'na'
#move taxid to row names
out_table_MZ_wide<-out_table_MZ_wide[complete.cases(out_table_MZ_wide), ]
row.names(out_table_MZ_wide) <- out_table_MZ_wide$TaxID
out_table_MZ_wide[1] <- NULL

#just double checking everything has a value
#out_table_MZ_wide2<-out_table_MZ_wide%>%filter(rowSums(across(where(is.numeric)))!=0)

#taxa table
#remove na
#move taxid to row.names
taxa_table_clean<-taxa_table_clean[rowSums(is.na(taxa_table_clean)) != ncol(taxa_table_clean), ]
row.names(taxa_table_clean) <- taxa_table_clean$TaxID
taxa_table_clean[1] <- NULL

#metadata UID to row.name
row.names(meta_data_clean) <- meta_data_clean$UID
meta_data_clean[1] <- NULL

#adding in root so everything as something (not just NA)
taxa_table_clean$root<-'Root'
taxa_table_clean <- taxa_table_clean %>%
  select(root, everything())

#setting up for phyloseq
tax_mat<-as.matrix(taxa_table_clean)
otu <- otu_table(out_table_MZ_wide, taxa_are_rows = TRUE) 
taxa <- tax_table(tax_mat)
sample <- sample_data(meta_data_clean)

####################################
#really overwhelming
#just cuting down to a couple of sites to make things manageable for now
#if you want to keep everything in just skip this and go down to the phyloseq section
####################################

#filtering out a couple of random sites
meta_data_cut <- meta_data %>% filter( grepl('Esk|Poukawa', HBRC_Site_Name))
out_table_MZ_cut <- out_table_MZ %>% filter(UID %in% meta_data_cut$UID)
#reverse in case there is one missing the other way for this marker
meta_data_cut_clean <- meta_data_cut %>% filter(UID %in% out_table_MZ$UID)

#cast wide
out_table_MZ_cut_wide<-dcast(out_table_MZ_cut, TaxID~UID,value.var = "Count", fun.aggregate = sum)

#going to have to remove that weird taxid == 'na'
#and fix row.names again
out_table_MZ_cut_wide<-out_table_MZ_cut_wide[complete.cases(out_table_MZ_cut_wide), ]
row.names(out_table_MZ_cut_wide) <- out_table_MZ_cut_wide$TaxID
out_table_MZ_cut_wide[1] <- NULL

#cleen up taxa table
taxa_table_cut_clean <- taxa_table %>% filter(TaxID %in% out_table_MZ_cut$TaxID)

taxa_table_cut_clean<-taxa_table_cut_clean[rowSums(is.na(taxa_table_cut_clean)) != ncol(taxa_table_cut_clean), ]
row.names(taxa_table_cut_clean) <- taxa_table_cut_clean$TaxID
taxa_table_cut_clean[1] <- NULL

row.names(meta_data_cut_clean) <- meta_data_cut_clean$UID
meta_data_cut_clean[1] <- NULL

taxa_table_cut_clean$root<-'Root'
taxa_table_cut_clean <- taxa_table_cut_clean %>%
  select(root, everything())

tax_mat<-as.matrix(taxa_table_cut_clean)
otu <- otu_table(out_table_MZ_cut_wide, taxa_are_rows = TRUE) 
taxa <- tax_table(tax_mat)
sample <- sample_data(meta_data_cut_clean)

###############################################
###phyloseq
###############################################

#some of this wont work on all sites and some wont work on cutdown just fucking about
FisheDNA<-phyloseq(otu, taxa, sample)
FisheDNA
sample_names(FisheDNA)[1:5]
rank_names(FisheDNA)
otu_table(FisheDNA)[1:5, 1:5]
tax_table(FisheDNA)[1:5, 1:4]
taxa_names(FisheDNA)[1:10]

FisheDNA.5 <- filter_taxa(FisheDNA, function(x){sum(x > 0) > 1}, prune=TRUE) #removing singletons
FisheDNA.5 #drop from 321 to 250 taxa

FisheDNA.5.fish = subset_taxa(FisheDNA.5, Class=="Actinopteri" |Class=="Cladistia" |Class=='Hyperoartia')
#7 fish taxa
FisheDNA.5 <-FisheDNA.5.fish

Fish.eDNA.pa <- microbiome::transform(FisheDNA.5, 'pa') #pa dataset
p = plot_bar(FisheDNA.5,  fill="Class")
p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Read abundance across sample number ")  #makes a more aesthetic plot

p1 = plot_bar(FisheDNA.5, "HBRC_Site_Name", fill="Class")
p1 + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Read abundance across HBRC_Sites") 
#
p1 = plot_bar(FisheDNA.5, "Cyclone", fill="Species")
p1 + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack") +
  ggtitle("Read abundance across cyclone") 

#Lets make some relative abundance plots
#First step is to transform to relative abundance
FisheDNA.6 = transform_sample_counts(Fish.eDNA.pa, function(x) x / sum(x) )
FisheDNA.6 

data_fish <- psmelt(FisheDNA.6)
relabundance_plot <- ggplot(data=data_fish, aes(x=Sample, y=Abundance, fill=Class)) + facet_grid(~HBRC_Site_Name, scales = "free")
relabundance_plot + geom_bar(aes(), stat="identity", position="fill") 

#Lets get these objects out of phyloseq and use them for further analysis

#First, the raw OTU table
Fish.OTU <- as(otu_table(FisheDNA.5), "matrix")
Fish_OTU_df = as.data.frame((Fish.OTU))
Fish_OTU_df= t(Fish_OTU_df)

#now presence-absence
Fish.pa.OTU <- as(otu_table(Fish.eDNA.pa), "matrix")
Fish.pa_OTU_df = as.data.frame((Fish.pa.OTU))
Fish.pa_OTU_df= t(Fish.pa_OTU_df)

Fish_Sample = as(sample_data(Fish.eDNA.pa), "matrix")
Fish_Sample_df =  as.data.frame(Fish_Sample)

Fish_taxa <- tax_table(Fish.eDNA.pa)
Fish_taxa_df <- as.data.frame(Fish_taxa)

Fish_PA_OTU_Sample = cbind(Fish_Sample_df, Fish.pa_OTU_df)

########################
## RAREFACTION CURVES ##
########################

## identify the lowers number of reads for samples and generate rarefaction curves
raremax_df <- min(rowSums(Fish_OTU_df))
#rarecurve(Fish_OTU_df, step = 100, sample = 1, col = 'blue', cex = 0.5)
#?rarecurve
rarecurve(Fish_OTU_df, step = 100, sample = 1, col = 'blue', cex = 0.5,label=FALSE)

######################
## SPECIES RICHNESS ##
######################
## calculate number of taxa detected per sample and group per sampling location 

rich<-estimate_richness(Fish.eDNA.pa, measures=c("Observed"))

rich<-cbind(rich, Fish_Sample)


## test assumptions of statistical test, first normal distribution, next homoscedasticity
histogram(~ Observed | location_cyclone, data = rich, layout = c(2,2))

shapiro.test(rich$Observed)

bartlett.test(Observed ~location_cyclone , data = rich)
#non-parametric - need to use Kruskal-Wallis test
kruskal.test(Observed ~ location_cyclone, data = rich)
boxplot(Observed ~ location_cyclone, data = rich, ylab = 'Species Richness', xlab = 'Location')


ggplot(rich, aes(x=location_cyclone, y=Observed, colour = location_cyclone)) +
  geom_boxplot() +
  ylab("Observed zotu richness")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#pairwise test - lets look at the difference in richness
pairwise.wilcox.test(rich$Observed, rich$location_cyclone,
                     p.adjust.method = "BH")

#iNEXT
library(iNEXT.3D)
categories <- unique(sample_data(Fish.eDNA.pa)$location_cyclone)
split_physeq_list <- list()
for (category in categories) {
  sub_physeq.100 <- subset_samples(Fish.eDNA.pa, location_cyclone == category)
  split_physeq_list[[category]] <- otu_table(sub_physeq.100)
}  
matrix_list <- lapply(split_physeq_list, function(x) {
  otu_table <- as(x, "matrix")
  return(otu_table)
})
matrix_list<-NULL
matrix_list <- list(data = list())
for (category in categories) {
  otu_table <- as(otu_table(split_physeq_list[[category]]), "matrix")
  matrix_list[["data"]][[category]] <- otu_table
}

out.raw <- iNEXT3D(data = matrix_list$data, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_raw', nboot = 50)
ggiNEXT3D(out.raw, type = 1, facet.var = 'Assemblage') + facet_wrap(~Assemblage, nrow = 3)
ggiNEXT3D(out.raw, type = 1, facet.var = "Order.q")
ggiNEXT3D(out.raw, type = 2, facet.var = "Order.q", color.var = "Assemblage")

#NMDS Plot


OTU_Dist <- vegdist(Fish.pa_OTU_df, method="jaccard") #create distance matrix

ord <- metaMDS(OTU_Dist,
               distance = "jaccard", #dissimilarity method
               k = 2, #number of dimensions
               maxit = 999, #max number of iterations 
               trymax = 500, #max number of random starts - may need to play around with this if you have many zeros in your dataset
               wascores = TRUE) #is a method of calculating species scores, default is TRUE

ord$stress #stress value

plot(ord)
# Shepards test/goodness of fit
goodness(ord) # Produces a results of test statistics for goodness of fit for each point
stressplot(ord) # Produces a Shepards diagram

both.scores = cbind(scores(ord), Fish_Sample_df) 

ggplot(both.scores) +
  geom_point( aes(x=NMDS1,y=NMDS2, colour = location_cyclone),size=2) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, colour= location_cyclone))+ 
  coord_equal() +
  ggtitle("NMDS") +
  theme_classic() 

#MDS through phyoseq
FisheDNA.ord <- ordinate(Fish.eDNA.pa, "NMDS", "jaccard") 
FisheDNA.ord
NMDSphylo = plot_ordination(Fish.eDNA.pa, FisheDNA.ord, type="samples", color="location_cyclone")
NMDSphylo + stat_ellipse(aes(x=NMDS1, y=NMDS2, colour= location_cyclone))+
  geom_point(size=2) +  coord_equal() +
  ggtitle("NMDS phyloseq") +
  theme_classic() 

#PCoA

ordination_mds <- wcmdscale(OTU_Dist, eig = TRUE)

pcoa_df <- data.frame(ordination_mds$points)
colnames(pcoa_df) <- c("PCo1", "PCo2")
pcoa_df$location_cyclone <- factor(Fish_Sample_df$location_cyclone) #add group of interest, mine was Location data
percent_explained <- 100 * ordination_mds$eig /sum(ordination_mds$eig)
pretty_pe <- round(percent_explained[1:2], digits = 1)
pretty_pe

labs <- c(glue("PCo1 ({pretty_pe[1]}%)"),
         glue("PCo2 ({pretty_pe[2]}%)"))


ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = location_cyclone)) + 
  geom_point(size = 2) +
  stat_ellipse(aes(x = PCo1, y = PCo2, color = location_cyclone)) +
  ggtitle("PCOA") +
  theme_classic() +
  labs(x=labs[1], y=labs[2])

ordination_eigen <- ordination_mds$eig
ordination_eigenvalue <- ordination_eigen/sum(ordination_eigen) 
ordination_eigen_frame <- data.frame(Inertia = ordination_eigenvalue*100, Axes = c(1:95)) #here we can see the axes we will be using and the eigen values
head(ordination_eigen_frame)

eigenplot <- ggplot(data = ordination_eigen_frame, aes(x = factor(Axes), y = Inertia)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
  theme_classic() +
  xlab("Axes") +
  ylab("Inertia %") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

eigenplot
ordination_mds$GOF #goodness of fit

#pcoa phyloseq
FisheDNA.PCoA = ordinate(Fish.eDNA.pa, "PCoA", "jaccard", weighted=TRUE)

PCoAphylo = plot_ordination(Fish.eDNA.pa, FisheDNA.PCoA, type="samples", color="location_cyclone")

PCoAphylo + stat_ellipse(aes(x=Axis.1, y=Axis.2, colour= location_cyclone))+
  geom_point(size=2) +  coord_equal() +
  ggtitle("PCoA phyloseq") +
  theme_classic() 

#fit enviro variables
meta_env <- meta_data_cut %>% select(-UID,-CollectionDate,-HBRC_Site_ID,-location_date)
en = envfit(ordination_mds, meta_env, permutations = 999, na.rm = TRUE)
en
{ordiplot(ordination_mds, display = 'sites')
  plot(en, p.max = 0.05, col = "red")}

#PERMANOVA

adonis2(OTU_Dist ~ location_cyclone, method = "jaccard", data = Fish_Sample_df, permutations = 9999)

#indication species
indval <- multipatt(Fish.pa_OTU_df, Fish_Sample_df$location_cyclone, 
                    control = how(nperm=999)) 

summary(indval)
tax_table(Fish.eDNA.pa)

#Group Esk Rvr at Berry Rd_Post+Esk Rvr at Waipunga Br_Post#+Poukawa Strm at Stock Rd_Post+Poukawa Strm Te Mahanga #Rd_Post  #sps.  1 
#stat p.value  
#2836 0.601   0.039 *
taxa_table_clean %>% filter(row.names(taxa_table_clean) %in% "2836")

#Group Esk Rvr at Berry Rd_Post+Poukawa Strm at Stock #Rd_Post+Poukawa Strm Te Mahanga Rd_Post  #sps.  1 
#stat p.value  
#1512276 0.582   0.018 *
taxa_table_clean %>% filter(row.names(taxa_table_clean) %in% "1512276")

#alright definitly going to have to filter some of this out. Could try:
#1- only analysing insecta
#2- only analysing true freshwater groups

########################
#JUST FRESHWATER GROUPS#
########################
#I've created a list of species I consider to have some evidence of freshwater life stages (most insects have freshwater and terrestrial stages so that could stuff things up)

freshwater_invertsfish<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')

freshwater_invertsfish %>% distinct(Rank) 
family<-freshwater_invertsfish %>% filter(Rank == 'Family')
genus<-freshwater_invertsfish %>% filter(Rank == 'Genus')
phylum<-freshwater_invertsfish %>% filter(Rank == 'Phylum')
order<-freshwater_invertsfish %>% filter(Rank == 'Order')
class<-freshwater_invertsfish %>% filter(Rank == 'Class')

records_eukaryote_freshwater<-records %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)

#>100K records for stuff that is eukaryote animal and posibly freshwater

otu_table_freshwater <- records_eukaryote_freshwater %>% select(UID, PrimerSet, Count, TaxID)

#split by primer set
otu_table_freshwater_splitout <- split(otu_table_freshwater, otu_table_freshwater$PrimerSet)

otu_table_freshwater_splitout_CI<- otu_table_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)

#BE - 1710 - 18S general eukaryote
#BU - 1723 - 18S general eukaryote
#BX - 723 - 18S general eukaryote
#CI - 1718 - COI insects, inverts, some fish
#DG - 652
#EA - 809
#GV - 728 - vascular plants ITS1
#LG - 987 - Mifish
#LV - 642
#MZ - 1708 - Vascular plants rbcl
#RV - 1684 - 12S vertebrate ecoprimers 
#TP - 1717 - vascular plants
#UM - 1723 -microbe
#WG - 727 0 venerid clams 16S
#WV - 1713 - 16S vertebrate primers (picks up several worms too)
#XG - 861
#YG - 298
#ZC - 642
#ZP - 586
#ZV - 80

#if we take the top 8 and filter out anything that doesnt have representation across all the groups

#bring in meta data
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)


#samples present and select from meta data
meta_data_clean <- meta_data %>% filter(UID %in% otu_table_freshwater_splitout_CI$UID)


#long format OTU table
#rows samples, columns tax, values counts
library(reshape2)
otu_table_freshwater_splitout_CI_wide<-dcast(otu_table_freshwater_splitout_CI, TaxID~UID,value.var = "Count", fun.aggregate = sum)

#so 1337 unique taxid's
#have just merged on taxid here (with aggregate sum)

#generate taxa table from TaxIDs
#from records file
taxa_table_eukaryote_freshwater<-records_eukaryote_freshwater %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()


length(unique(otu_table_freshwater_splitout_CI$TaxID))
#1337 but one is a NA

#might need to filter by taxid's in otu table
taxa_table_eukaryote_freshwater_clean <- taxa_table_eukaryote_freshwater %>% filter(TaxID %in% otu_table_freshwater_splitout_CI$TaxID)


#generate meta file
#critical one
#Location_date
#Location_precyclone , Location_postcyclone
head(meta_data)
meta_data<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))

meta_data_clean <- meta_data %>% filter(UID %in% otu_table_freshwater_splitout_CI$UID)

#okay so my three correctly shaped files
#out_table_MZ_wide
#taxa_table_clean
#meta_data_clean

#going to have to remove that weird taxid == 'na'
otu_table_freshwater_splitout_CI_wide<-otu_table_freshwater_splitout_CI_wide[complete.cases(otu_table_freshwater_splitout_CI_wide), ]
row.names(otu_table_freshwater_splitout_CI_wide) <- otu_table_freshwater_splitout_CI_wide$TaxID
otu_table_freshwater_splitout_CI_wide[1] <- NULL

#just double checking everything has a value
otu_table_freshwater_splitout_CI_wide2<-otu_table_freshwater_splitout_CI_wide%>%filter(rowSums(across(where(is.numeric)))!=0)

taxa_table_eukaryote_freshwater_clean<-taxa_table_eukaryote_freshwater_clean[rowSums(is.na(taxa_table_eukaryote_freshwater_clean)) != ncol(taxa_table_eukaryote_freshwater_clean), ]
row.names(taxa_table_eukaryote_freshwater_clean) <- taxa_table_eukaryote_freshwater_clean$TaxID
taxa_table_eukaryote_freshwater_clean[1] <- NULL

row.names(meta_data_clean) <- meta_data_clean$UID
meta_data_clean[1] <- NULL

taxa_table_eukaryote_freshwater_clean$root<-'Root'
taxa_table_eukaryote_freshwater_clean <- taxa_table_eukaryote_freshwater_clean %>%
  select(root, everything())

tax_mat<-as.matrix(taxa_table_eukaryote_freshwater_clean)
otu <- otu_table(otu_table_freshwater_splitout_CI_wide, taxa_are_rows = TRUE) 
taxa <- tax_table(tax_mat)
sample <- sample_data(meta_data_clean)

#should go into phyloseq now

###############
#create a map#
##############

library(leaflet)
library(htmltools)
library(tidyverse)
library(RColorBrewer)

#load the metadata file
meta_data_map<-read.table('Wilderlab_meta_map.txt',sep='\t',quote='',header =T)


leaflet(meta_data_map) %>% addTiles() %>%
  addMarkers(~Longitude_hbrc, ~Latitude_hbrc, popup = ~paste(HBRC_Site_Name, paste('N Uniq Sample Dates:', n_uniq_samplingdates,sep =' '), sep = "<br>"))


pal <- colorFactor(
  palette = 'Dark2',
  domain = meta_data_map$n_uniq_samplingdates
)

leaflet(meta_data_map) %>% 
  addTiles() %>%
  addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                   popup = ~paste(HBRC_Site_Name),
                    color = ~ pal(n_uniq_samplingdates),
                   stroke = FALSE, fillOpacity = 1) %>% 
  addLegend("bottomright", pal = pal, values = ~n_uniq_samplingdates,
            title = "# unique sampling dates",
            opacity = 1) 