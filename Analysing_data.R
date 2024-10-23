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
test<-otu_table %>% filter(PrimerSet=='CI')
test2<-otu_table %>% filter(PrimerSet=='CI')%>% select(-PrimerSet) %>% group_by(UID, TaxID)


#split by primer set
otu_table_splitout <- split(otu_table, otu_table$PrimerSet)
#so there is 20 different primer sets
#no idea what each of them is
#otu_table_splitout$

#otu_table_splitout %>% (filter %in%'CI')
#just work through one primer set for now and then come back 

#just choosing a random not super diverse but present in most samples to play with
head(otu_table_splitout$MZ)

out_table_MZ<- otu_table_splitout$MZ %>% select(-PrimerSet) %>% group_by(UID, TaxID)
out_table_MZ<- otu_table_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)

#switching to CI which is insect
out_table_MZ<- otu_table_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)
#switching to RV which has issues in app
out_table_MZ<- otu_table_splitout$RV %>% select(-PrimerSet) %>% group_by(UID, TaxID)

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
meta_data_cut <- meta_data %>% filter( grepl('Esk Rvr at Berry Rd', HBRC_Site_Name))
#weird sites acting out
meta_data_cut <- meta_data %>% filter( grepl('Wharerangi US Ahuriri Estuary', HBRC_Site_Name))
meta_data_cut <- meta_data %>% filter( grepl('Awanui Strm at Flume', HBRC_Site_Name))
meta_data_cut <- meta_data %>% filter( grepl('Aropaoanui at Aropaoanui', HBRC_Site_Name))


#othero one is Esk Rvr at Waipunga Br
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
sample_data(FisheDNA)
filter<-as.data.frame(sample_data(FisheDNA)) %>% group_by(CollectionDate) %>% summarise(count=n()) %>% filter(count >=5) 
test<-ps_filter(FisheDNA,CollectionDate != "16/11/23")
sample_names(FisheDNA)
otu_table(FisheDNA)
FisheDNA
subset_samples(FisheDNA, CollectionDate == "16/11/23")
#this:
subset_samples(FisheDNA, CollectionDate %in% filter$CollectionDate)


data("enterotype", package = "phyloseq")
enterotype
sample_data(enterotype)
ps1 <- ps_filter(enterotype, SeqTech != "Sanger")
ps1

test<-meta_data_clean %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name=='Esk Rvr at Berry Rd') %>% pull(sample_id)
testphyloseq<-prune_samples(test,FisheDNA)
#order samples by date
#modify for which metadata using:
meta_data_clean
test<-meta_data_cut_clean[order(as.Date(meta_data_cut_clean$CollectionDate, format="%d/%m/%Y")),]
test<-meta_data_clean[order(as.Date(meta_data_clean$CollectionDate, format="%d/%m/%Y")),]
test_names<-row.names(test)
library(microViz)
?ps_reorder
FisheDNA<-ps_reorder(FisheDNA, test_names)
#looks reordered now

FisheDNA
sample_names(FisheDNA)[1:5]
rank_names(FisheDNA)
otu_table(FisheDNA)[1:5, 1:5]
tax_table(FisheDNA)[1:5, 1:4]
taxa_names(FisheDNA)[1:10]

FisheDNA.5 <- filter_taxa(FisheDNA, function(x){sum(x > 0) > 1}, prune=TRUE) #removing singletons
FisheDNA.5 #drop from 321 to 250 taxa

FisheDNA.5.fish = subset_taxa(FisheDNA.5, Class=="Actinopteri" |Class=="Cladistia" |Class=='Hyperoartia')
otu_table(FisheDNA.5.fish)
#7 fish taxa
#FisheDNA.5 <-FisheDNA.5.fish

Fish.eDNA.pa <- microbiome::transform(FisheDNA.5, 'pa') #pa dataset
library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
?plot_bar
p = plot_bar(FisheDNA.5,  fill="Class")
p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Read abundance across sample number ") +
  geom_vline(aes(xintercept = 646.5), 
             linetype = "dashed", colour = "red",size = 1)+
  theme(legend.key.size = unit(0.05, 'cm'))+
 theme(legend.position="none")#+
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)

  p$data
  
  ptest = plot_bar(FisheDNA.5, "CollectionDate", fill="Class") + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") + ggtitle("Read Abundance") + theme_bw()+  theme(legend.position="none")
  ptest
  test_data<-ptest$data
  ptest$plot_env
  
  ggplot(test_data, aes(fill=Class, y=Abundance, x=CollectionDate)) + 
    geom_bar(position="stack", stat="identity")
  
  #labeling only the top ten classes
  top_10class<-test_data %>% group_by(Class) %>% drop_na(Class) %>%  summarise(Abundance = sum(Abundance,na.rm=TRUE)) %>% arrange(desc(Abundance)) %>% slice_max(Abundance, n = 10,with_ties=FALSE)
      data_fish_top10 <- test_data %>% 
       mutate(Class_top10 = ifelse (!Class %in% top_10class$Class, "Other", Class))
  #make a colour pallet of 11 colours (grey at the end for other)
    col_brew<-c(brewer.pal(n = 10, name = "Paired"),"#808080")
  #relevel other to end
   data_fish_top10$Class_top10 <- forcats::fct_relevel(data_fish_top10$Class_top10, "Other", after = Inf)
  
   ggplot(data_fish_top10, aes(fill=Class_top10, y=Abundance, x=CollectionDate)) + 
     geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=col_brew)+ guides(fill=guide_legend(title="Top 10 Classes"))+ theme(legend.title = element_text(size = 8),legend.text = element_text(size = 6)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Read Abundance")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   
#presence absence plot
  
  p = plot_bar(Fish.eDNA.pa,  fill="Class")
  p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
    ggtitle("Presence/Absence across sites") +
    geom_vline(aes(xintercept = 646.5), 
               linetype = "dashed", colour = "red",size = 1)+
    theme(legend.key.size = unit(0.05, 'cm'))#+
    theme(legend.position="none")

#adding a vertical line to indicate cyclone date 
  #makes a more aesthetic plot
#such a clunky way to do colours, not sure why I cant set it in main 
#colours will only work for subset of freshwater (14 classes)

Fish.eDNA.pa <- microbiome::transform(FisheDNA.5, 'pa') #pa dataset
?plot_bar
p = plot_bar(FisheDNA.5,  fill="Class")
p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Read abundance across sample number ")  #makes a more aesthetic plot

p1 = plot_bar(FisheDNA.5, "location_date", fill="Class")
#do presence/absence plot
#p1 = plot_bar(Fish.eDNA.pa, "location_date", fill="Class")

unique(p1$data$location_date)
#relevel to order chart
#p1$data$location_date <- factor(p$data$location_date, levels = c('Esk Rvr at Berry Rd_24/02/22',"Esk Rvr at Berry Rd_20/04/23","Esk Rvr at Berry Rd_30/08/23","Esk Rvr at Berry Rd_16/11/23","Esk Rvr at Berry Rd_24/01/24"))

p1$data$location_date <- factor(p1$data$location_date, levels = c('Esk Rvr at Berry Rd_24/02/22',"Esk Rvr at Berry Rd_20/04/23","Esk Rvr at Berry Rd_30/08/23","Esk Rvr at Berry Rd_16/11/23","Esk Rvr at Berry Rd_24/01/24"))
p1 + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Presence/absence across HBRC_Sites") +
#  ggtitle("Read abundance across HBRC_Sites")+
  theme(legend.key.size = unit(0.03, 'cm'))+theme(axis.text=element_text(size=4),axis.title=element_text(size=10,face="bold"))

#
p1 = plot_bar(FisheDNA.5, "Cyclone", fill="Species")
p1 + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack") +
  ggtitle("Read abundance across cyclone") 

#Lets make some relative abundance plots
#First step is to transform to relative abundance
FisheDNA.6 = transform_sample_counts(Fish.eDNA.pa, function(x) x / sum(x) )
FisheDNA.6 
#or from raw data: probably more sensible but shit looking
FisheDNA.6 = transform_sample_counts(FisheDNA.5, function(x) x / sum(x) )
FisheDNA.6 

data_fish <- psmelt(FisheDNA.6)
#data_fish$CollectionDate <- factor(data_fish$CollectionDate, levels = c('24/02/22',"20/04/23","30/08/23","16/11/23","24/01/24"))
data_fish<-data_fish[order(as.Date(data_fish$CollectionDate, format="%d/%m/%Y")),]
relabundance_plot <- ggplot(data=data_fish, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(~Year, scales = "free")
relabundance_plot + geom_bar(aes(), stat="identity", position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #theme(legend.position="none")
#+
  theme(legend.key.size = unit(0.03, 'cm'))


#find top 10 classes in dataset
top_10class<-data_fish %>% group_by(Class) %>% drop_na(Class) %>%  summarise(Abundance = sum(Abundance,na.rm=TRUE)) %>% arrange(desc(Abundance)) %>% slice_max(Abundance, n = 10)
top_10class$Class
#then everything that is not in top_10class gets mutated into 'other' category
#https://stackoverflow.com/questions/71595962/conditional-mutate-by-matching-strings-or-characters

data_fish_top10 <- data_fish %>% 
  mutate(Class_top10 = ifelse (!Class %in% top_10class$Class, "Other", Class))
relabundance_plot <- ggplot(data=data_fish_top10, aes(x=Sample, y=Abundance, fill=Class_top10)) + facet_grid(~Year, scales = "free")
relabundance_plot + geom_bar(aes(), stat="identity", position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

col_brew<-c(brewer.pal(n = 10, name = "Paired"),"#808080")
col_brew

data_fish_top10$Class_top10 <- forcats::fct_relevel(data_fish_top10$Class_top10, "Other", after = Inf)

relabundance_plot <- ggplot(data=data_fish_top10, aes(x=Sample, y=Abundance, fill=Class_top10)) + facet_grid(~Year, scales = "free")
relabundance_plot + geom_bar(aes(), stat="identity", position="fill") +theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=col_brew)+ guides(fill=guide_legend(title="Top 10 Classes"))+ theme(legend.title = element_text(size = 8),legend.text = element_text(size = 6)) +theme(legend.key.height=unit(0.3, "cm"))

?slice_max
pbar = plot_bar(FisheDNA.5, "CollectionDate", fill="Class")
#do presence/absence plot
pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"))

pbar + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Presence/absence") +
  #  ggtitle("Read abundance across HBRC_Sites")+
  theme(legend.key.size = unit(0.03, 'cm'))+theme(axis.text=element_text(size=4),axis.title=element_text(size=10,face="bold"))
pbar

#desired_order <- names(sort(meta_data_cut_clean$CollectionDate, decreasing=TRUE))


p1 = plot_bar(FisheDNA.5, "location_date", fill="Class")
unique(p1$data$location_date)
p1$data$location_date <- factor(p1$data$location_date, levels = c('Esk Rvr at Berry Rd_24/02/22',"Esk Rvr at Berry Rd_20/04/23","Esk Rvr at Berry Rd_30/08/23","Esk Rvr at Berry Rd_16/11/23","Esk Rvr at Berry Rd_24/01/24"))
p1 + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  ggtitle("Presence/absence across HBRC_Sites") +
  #  ggtitle("Read abundance across HBRC_Sites")+
  theme(legend.key.size = unit(0.03, 'cm'))+theme(axis.text=element_text(size=4),axis.title=element_text(size=10,face="bold"))



sort(as.POSIXct(unique(pbar$data$CollectionDate)))
format(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y"), "20%y-%m-%d")
sort(format(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y"), "20%y-%m-%d"))
dplyr::arrange(data, as.POSIXct(pbar$data$CollectionDate))
pbar$data$Sample <- factor(pbar$data$Sample, levels = desired_order)
pbar

#turning everything under say 5% to other and to cut down on the legend
plot_5<-transform_sample_counts(FisheDNA.5, function(x) 100 * x/sum(x))
tax_table(plot_5)
#FisheDNA.6 = transform_sample_counts(FisheDNA.5, function(x) x / sum(x) )
glom <- tax_glom(plot_5, taxrank = 'Class')
glom # should list taxa as phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Class <- as.character(data_glom$Class) #convert to character
data_glom$Class[data_glom$Abundance < 0.25] <- "<25% abundance"
Count = length(unique(data_glom$Class)); Count
unique(data_glom$Class) #add into line below

spatial_plot <- ggplot(data_glom, aes(x=Sample, y=Abundance, fill=Class)) + facet_grid(~location_date, scales = "free")

spatial_plot + geom_bar(aes(), stat="identity", position="stack") + labs (title = "Relative abudance", x="Stage of Decomposition", y= "Relative Abundance", ) + theme(legend.position="bottom", plot.title = element_text(hjust=0.5)) + guides(fill=guide_legend(nrow=4))

names(sort(taxa_sums(FisheDNA.5), TRUE)[1:10])


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
histogram(~ Observed | location_date, data = rich, layout = c(2,2))

shapiro.test(rich$Observed)

bartlett.test(Observed ~location_date , data = rich)
#non-parametric - need to use Kruskal-Wallis test
kruskal.test(Observed ~ location_date, data = rich)
boxplot(Observed ~ location_date, data = rich, ylab = 'Species Richness', xlab = 'Location')

rich$location_date

rich$location_date <- factor(rich$location_date, levels = c('Esk Rvr at Berry Rd_24/02/22',"Esk Rvr at Berry Rd_20/04/23","Esk Rvr at Berry Rd_30/08/23","Esk Rvr at Berry Rd_16/11/23","Esk Rvr at Berry Rd_24/01/24"))

ggplot(rich, aes(x=location_date, y=Observed, colour = location_date)) +
  geom_boxplot() +
  ylab("Observed zotu richness")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1))+theme(axis.text=element_text(size=4),axis.title=element_text(size=10,face="bold"))

#pairwise test - lets look at the difference in richness
pairwise.wilcox.test(rich$Observed, rich$location_date,
                     p.adjust.method = "BH")

#iNEXT
library(iNEXT.3D)
categories <- unique(sample_data(Fish.eDNA.pa)$location_date)
split_physeq_list <- list()
for (category in categories) {
  sub_physeq.100 <- subset_samples(Fish.eDNA.pa, location_date == category)
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

innext_plot<-ggiNEXT3D(out.raw, type = 1, facet.var = 'Assemblage') + facet_wrap(~Assemblage, nrow = 3)
innext_plot
innext_plot$theme
?innext_plot

plotdf<-fortify.iNEXT(out.raw, type=1)


innext_plot$data$Assemblage <- factor(innext_plot$data$Assemblage, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))
out.raw$TDInfo$Assemblage<- factor(out.raw$TDInfo$Assemblage, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))
factor(out.raw$TDInfo$Assemblage)

innext_plot<-ggiNEXT3D(out.raw, type = 1, facet.var = "Order.q")
innext_plot+ theme(plot.title = element_text(size=10),text = element_text(size = 10))+ theme(legend.position="none")
innext_plot$layers

#this does not solve all the issues for some reason I cant change the factor of the shape but can the colour so shifting the colour over to ordered and just turning of the shape legend. 
innext_plot$data$shape <- factor(innext_plot$data$shape, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))
#have tried the same on Assemblage and shape columns in data table but no dice, not sure why its not working. Also tried changing factor of inext.3d object out.raw and didnt work either
innext_plot+ theme(plot.title = element_text(size=8),text = element_text(size = 8))+theme(legend.position = "right") +theme(legend.text=element_text(size=6))+guides(shape = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +theme(legend.key.size = unit(0.5,"line")) 
#cant seem to change the legend size without stuffing up the rest of the legend

+ guides(colour = guide_legend(override.aes = list(size=2)))

+ guides(colour = guide_legend(override.aes = list(alpha=1)))

+ guides(linetype = guide_legend(override.aes = list(size = 0.5)))

+guides(col = guide_legend(override.aes = list(size = 0.5)))
  
  theme(legend.key.size = unit(0.5,"line"))



innext_plot$data$Assemblage <- factor(innext_plot$data$Assemblage, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))

innext_plot$data$shape <- factor(innext_plot$data$shape, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))

innext_plot$data$col <- factor(innext_plot$data$col, levels = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))

+geom_line(size=1)

+ scale_fill_discrete(limits = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))+ scale_fill_discrete(limits = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))



innext_plot+ theme(plot.title = element_text(size=10),text = element_text(size = 10),axis.text.x=element_blank())

innext_plot+ theme(legend.title=element_text(size=6))+theme(axis.text=element_text(size=6))
innext_plot+ scale_fill_discrete(limits = c('Esk Rvr at Berry Rd_24/02/22','Esk Rvr at Berry Rd_20/04/23','Esk Rvr at Berry Rd_30/08/23','Esk Rvr at Berry Rd_16/11/23','Esk Rvr at Berry Rd_24/01/24'))

unique(innext_plot$data$Assemblage)

ggiNEXT3D(out.raw, type = 2, facet.var = "Order.q", color.var = "Assemblage")

ggiNEXT3D(out.raw, type = 2, facet.var = "Order.q", color.var = "Assemblage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.key.size = unit(0.03, 'cm'))+theme(axis.text=element_text(size=4),axis.title=element_text(size=4,face="bold"),text=element_text(size=6))+  theme(legend.key.size = unit(0.05, 'cm'))


?ggiNEXT3D
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
  geom_point( aes(x=NMDS1,y=NMDS2, colour = location_date),size=2) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, colour= location_date))+ 
  coord_equal() +
  ggtitle("NMDS") +
  theme_classic() 

ggplot(both.scores) +
  geom_point( aes(x=NMDS1,y=NMDS2, colour = location_date),size=2) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, colour= location_date))+ 
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
#pernova in app
phyloseq_object_prune.pa <- microbiome::transform(FisheDNA, 'pa')
phyloseq_object_prune.pa.OTU <- as(otu_table(phyloseq_object_prune.pa), "matrix")
phyloseq_object_prune.pa_OTU_df = as.data.frame(phyloseq_object_prune.pa.OTU)
phyloseq_object_prune.pa_OTU_df= t(phyloseq_object_prune.pa_OTU_df)
phyloseq_object_prune.pa_Sample = as(sample_data(phyloseq_object_prune.pa), "matrix")
phyloseq_object_prune.pa_Sample_df =  as.data.frame(phyloseq_object_prune.pa_Sample)
OTU_Dist <- vegdist(phyloseq_object_prune.pa_OTU_df, method="jaccard")
mod_beta<-adonis2(OTU_Dist ~ location_date, method = "jaccard", data = phyloseq_object_prune.pa_Sample_df, permutations = 9999)
print(mod_beta)
print(OTU_Dist)

stargazer(mod_beta)

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

otu_table_freshwater_splitout_CI<- otu_table_freshwater_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)

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
unique(otu_table_freshwater_splitout_CI$TaxID)
otu_table_freshwater_splitout_CI_wide<-dcast(otu_table_freshwater_splitout_CI, TaxID~UID,value.var = "Count", fun.aggregate = sum)

#so 1337 unique taxid's
#have just merged on taxid here (with aggregate sum)

#generate taxa table from TaxIDs
#from records file
taxa_table_eukaryote_freshwater<-records_eukaryote_freshwater %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
#so this has 499 unique taxids but otu_table_freshwater_splitout_CI_wide has 1337
test<-row.names(otu_table_freshwater_splitout_CI_wide)
test
records_eukaryote_freshwater %>% filter(TaxID %in% test) %>% select(TaxID) %>% distinct() %>% nrow()

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

nrow(taxa_table_eukaryote_freshwater_clean)
tax_mat<-as.matrix(taxa_table_eukaryote_freshwater_clean)
otu <- otu_table(otu_table_freshwater_splitout_CI_wide, taxa_are_rows = TRUE) 
taxa <- tax_table(tax_mat)
sample <- sample_data(meta_data_clean)

#add the seasonality to the map meta_data
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)



test<-meta_data %>% mutate(Date = as.Date(CollectionDate, format= "%d/%m/%y")) %>% mutate(season=case_when(Date>='2019-07-01' & Date<='2020-06-30' ~'2019-2020',Date>='2020-07-01' & Date<='2021-06-30' ~'2020-2021',Date>='2021-07-01' & Date<='2022-06-30' ~'2021-2022', Date>='2022-07-01' & Date<='2023-06-30' ~'2022-2023', Date>='2023-07-01' & Date<='2024-06-30' ~'2023-2024'))

unique(test$season)
#another column pre and post cyclone
test %>% select(CollectionDate,Cyclone,season) %>% unique()
#Cyclone column is correct in pre/post
#all 2022-2023 are post cyclone so to compare first year impact to second year impact can do: 2022-2023 vs 2023-2024
write.table(test,'Wilderlab_meta_out_season.txt',sep='\t',row.names=F)
test<-read.table('Wilderlab_meta_out_season.txt',sep='\t',header =T)



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
meta_data_map<-read.table('OneDrive_1_7-05-2024/Wilderlab_meta_map.txt',sep='\t',quote='',header =T)

leaflet(meta_data_map) %>% addTiles() %>%
  addMarkers(~Longitude_hbrc, ~Latitude_hbrc, popup = ~paste(HBRC_Site_Name, paste('N Uniq Sample Dates:', n_uniq_samplingdates,sep =' '), sep = "<br>"))


pal <- colorFactor(
  palette = 'Blues',
  domain = meta_data_map$n_uniq_samplingdates
)

my_palette <- brewer.pal(name="YlOrBr",n=9)[4:9]

pal <- colorFactor(
  palette = my_palette,
  domain = meta_data_map$n_uniq_samplingdates
)

leaflet(meta_data_map) %>% 
  addTiles() %>%
  addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                   popup = ~paste(HBRC_Site_Name),
                   color='black',
                    fillColor = ~ pal(n_uniq_samplingdates),
                   stroke = FALSE, fillOpacity = 1,radius =6) %>% 
  addProviderTiles("PenStreetMap.France") %>%
  addLegend("bottomright", pal = pal, values = ~n_uniq_samplingdates,
            title = NULL,
            opacity = 1) 
leafletOptions(resolutions = 1200)
?leafletOptions
?leaflet
pal <- colorFactor(
  palette = my_palette,
  domain = meta_data_map$n_uniq_samplingdates
)

leaflet(meta_data_map) %>% 
  addTiles() %>%
  addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                   popup = ~paste(HBRC_Site_Name),
                   color='black',
                   fillColor = ~ my_palette,
                   stroke = FALSE, fillOpacity = 1,radius =6) %>% 
  addLegend("bottomright", colors = my_palette, values = ~n_uniq_samplingdates,
            title = "# unique sampling dates",
            opacity = 1) 
length(unique(meta_data_map$n_uniq_samplingdates))
#toanga species map#
#just from everything?
#just plot out presence before and after
head(records_eukaryote_freshwater)

#short finned eel
#Anguilla australis
#long finned eel
#Anguilla dieffenbachii
#freshwater cray genus: 
#Paranephrops 
test<-left_join(records_eukaryote_freshwater,meta_data,by='UID')

#year2022<-test %>% filter(Species=='Anguilla dieffenbachii' & Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
#year2022<-test %>% filter(Genus=='Nesameletus' & Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
year2022<-test %>% filter(Family=='Nesameletidae' & Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()

#year2023<-test %>% filter(Species=='Anguilla dieffenbachii' & Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()
#year2023<-test %>% filter(Genus=='Nesameletus' & Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()
year2023<-test %>% filter(Family=='Nesameletidae' & Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()

#year2024<-test %>% filter(Species=='Anguilla dieffenbachii' & Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
#year2024<-test %>% filter(Genus=='Nesameletus' & Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
year2024<-test %>% filter(Family=='Nesameletidae' & Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
#sites sampled each year

year2022_allsites<-test %>% filter(Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
year2023_allsites<-test %>% filter( Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()
year2024_allsites<-test %>% filter(Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()

year2022_allsites<-year2022_allsites %>% mutate(present = case_when(HBRC_Site_Name %in% year2022$HBRC_Site_Name~"present",TRUE ~ 'Absent'))

year2023_allsites<-year2023_allsites %>% mutate(present = case_when(HBRC_Site_Name %in% year2023$HBRC_Site_Name~"present",TRUE ~ 'Absent'))

year2024_allsites<-year2024_allsites %>% mutate(present = case_when(HBRC_Site_Name %in% year2024$HBRC_Site_Name~"present",TRUE ~ 'Absent'))

year2023_allsites<-year2023_allsites %>% mutate(present_2022 = case_when(HBRC_Site_Name %in% year2022$HBRC_Site_Name~"present",TRUE ~ 'Absent'))

year2023_allsites<-year2023_allsites %>% mutate(present_overtime = case_when(
  (present == 'Absent' & present_2022 == 'Absent') ~ 'Absent',
  (present == 'present' & present_2022 == 'present') ~ 'Present_2022&2023',
  (present == 'present' & present_2022 == 'Absent') ~ 'Gain 2023',
  (present == 'Absent' & present_2022 == 'present') ~ 'Loss 2023', TRUE ~ NA_character_
))

meta_data_map<-read.table('Wilderlab_meta_map.txt',sep='\t',quote='',header =T)

year2023_allsites<-left_join(year2023_allsites,meta_data_map,by='HBRC_Site_Name')

pal <- colorFactor(
  palette = 'Dark2',
  domain = year2023_allsites$present_overtime
)

leaflet(year2023_allsites) %>% 
  addTiles() %>%
  addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                   popup = ~paste(HBRC_Site_Name),
                   color = ~ pal(present_overtime),
                   stroke = FALSE, fillOpacity = 1) %>% 
  addLegend("bottomright", pal = pal, values = ~present_overtime,
            title = "Presence/absence",
            opacity = 1) 

#######
#fixing site list for app
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
actualsites<-read.csv('Site_list_HBRC_Sep2024_justin.csv',header=T)
unique(meta_data$HBRC_Site_Name)
sort(unique(meta_data$HBRC_Site_Name))
#should only be 93
sites<-actualsites$HBRC_Site_Name
sites
meta_data_justsites<-meta_data %>% filter(HBRC_Site_Name %in% sites)
write.table(meta_data_justsites,'Wilderlab_meta_out_clean.txt',row.names=F,sep='\t',quote=F)
