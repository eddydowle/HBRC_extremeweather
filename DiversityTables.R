#Eddy Nov 2024
#code to build diversity tables to then work out diversity change between sites

#first build diversity tables

#for each assay create a giant phyloseq table
#find unique sample names in phyloseq table
#make loop
#for each sample name find the unique sampling dates
#for each unique sampling date generate a diversity estimate
#create a output that is season | diversity score
#build it into a table of 
#Season diversiiy location 
#season2020-2021 1.4 loc_x
#season2021-2022 1.6 loc_x
#season2022-2023 1.8 loc_x
#season2023-2024 1.9

#merge into a larger table that is rows season column diversity score per site

#then should be able to grab a row of diversity scores per year for the comparison

#observed diversity
#shannon
#simpson

#only do this on insect and fish sets and filter by freshwater

#use raw table rather than pa for shannon and simpson

library('tidyverse')
library('readxl')
library('devtools')
library('tidyverse')
library('viridis')
library('vegan')
library('glue')
library('reshape2')
library('phyloseq')
library('microbiome')


#set up phyloseq table as previously done for wilderlabs

setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
records<-read.csv('HBRC_records_Eddy.csv')
head(records)
otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
otu_table_splitout <- split(otu_table, otu_table$PrimerSet)

#will need to do this for 3 freshwater animal sets CI, RV and WV
#out_table_MZ<- otu_table_splitout$CI %>% select(-PrimerSet) %>% group_by(UID, TaxID)
#out_table_MZ<- otu_table_splitout$RV %>% select(-PrimerSet) %>% group_by(UID, TaxID)
out_table_MZ<- otu_table_splitout$WV %>% select(-PrimerSet) %>% group_by(UID, TaxID)

#freshwater table
freshwater_taxa<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')

#bring in meta data
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
#samples present and select from meta data
meta_data_clean <- meta_data %>% filter(UID %in% out_table_MZ$UID)

#cast otu table wide
out_table_MZ_wide<-dcast(out_table_MZ, TaxID~UID,value.var = "Count", fun.aggregate = sum)

#generate taxa table from TaxIDs
#from records file
taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
#filter freshwater
family<-freshwater_taxa %>% filter(Rank == 'Family')
genus<-freshwater_taxa %>% filter(Rank == 'Genus')
phylum<-freshwater_taxa %>% filter(Rank == 'Phylum')
order<-freshwater_taxa %>% filter(Rank == 'Order')
class<-freshwater_taxa %>% filter(Rank == 'Class')
taxa_table<-taxa_table %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)
#filter by taxid's in otu table
taxa_table_clean <- taxa_table %>% filter(TaxID %in% out_table_MZ$TaxID)

#clean metadata
meta_data<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))
meta_data_clean <- meta_data %>% filter(UID %in% out_table_MZ$UID)

#move taxid to row names
out_table_MZ_wide<-out_table_MZ_wide[complete.cases(out_table_MZ_wide), ]
row.names(out_table_MZ_wide) <- out_table_MZ_wide$TaxID
out_table_MZ_wide[1] <- NULL
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
FisheDNA<-phyloseq(otu, taxa, sample)
FisheDNA

#generate sample list for loops need site names a dataframe with collection dates
samp<-as.data.frame(sample_data(FisheDNA)) #includes dates
sample_list<-unique(samp$HBRC_Site_Name) #just site names for loop


#here is the loop to make the file and the reasoning is below
outtab<-NULL

for (item in sample_list){
  print(item)
  to_prune<-meta_data_clean %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==item) %>% pull(sample_id)
  print(to_prune)
  FisheDNA_cut<-prune_samples(to_prune,FisheDNA)
  sample_data(FisheDNA_cut)
  collection_cut<-samp %>% as.matrix %>% as.data.frame() %>% filter(HBRC_Site_Name==item) %>% select(CollectionDate) %>% unique() %>% pull(CollectionDate)
  print(collection_cut)
  for (subdate in collection_cut){
    to_prune_within<-meta_data_clean %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==item & CollectionDate==subdate) %>% pull(sample_id)
    print(to_prune_within)
    FisheDNA_cut<-prune_samples(to_prune_within,FisheDNA)
    sample_data(FisheDNA_cut)
    phyloseq_object_prune.pa <- microbiome::transform(FisheDNA_cut, 'identity') #raw
    #      phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa')#pa dataset
    rich<-estimate_richness(phyloseq_object_prune.pa, measures=c("Shannon", "Simpson", "Observed"))
    print(rich)
    date_selected<-as.Date(subdate,format="%d/%m/%y")
    print(date_selected)
    if (date_selected>='2019-11-01' & date_selected<='2020-05-31'){
      print('2019-2020_summer')
      season<-'2019-2020_summer'
    }
    if (date_selected>='2019-06-01' & date_selected<='2019-10-31'){
      print('2019_winter')
      season<-'2019_winter'
    }
    if (date_selected>='2020-11-01' & date_selected<='2021-05-31'){
      print('2020-2021_summer')
      season<-'2020-2021_summer'
    }
    if (date_selected>='2020-06-01' & date_selected<='2020-10-31'){
      print('2020_winter')
      season<-'2020_winter'
    }
    if (date_selected>='2021-11-01' & date_selected<='2022-05-31'){
      print('2021-2022_summer')
      season<-'2021-2022_summer'
    }
    if (date_selected>='2021-06-01' & date_selected<='2021-10-31'){
      print('2021_winter')
      season<-'2021_winter'
    }
    if (date_selected>='2022-11-01' & date_selected<='2023-05-31'){
      print('2022-2023_summer')
      season<-'2022-2023_summer'
    }
    if (date_selected>='2022-06-01' & date_selected<='2022-10-31'){
      print('2022_winter')
      season<-'2022_winter'
    }
    if (date_selected>='2023-11-01' & date_selected<='2024-05-31'){
      print('2023-2024_summer')
      season<-'2023-2024_summer'
    }
    if (date_selected>='2023-06-01' & date_selected<='2023-10-31'){
      print('2023_winter')
      season<-'2023_winter'
    }
    rich_mod<-  colMeans(rich) %>% as.data.frame() %>% rename_with(~season) %>% rownames_to_column("Richness_test") %>%  pivot_longer(cols=c(-Richness_test),names_to="Season") %>%  pivot_wider(names_from=c(Richness_test)) %>% mutate(location=item) %>% mutate(date_col=subdate)
    print(rich_mod)
    outtab<-bind_rows(outtab,rich_mod)
  }
}

#write.table(outtab,'Diversity_estimates_persite_assay_CI_rawotu.txt',row.names=F,quote=F,sep='\t')
#write.table(outtab,'Diversity_estimates_persite_assay_RV_rawotu.txt',row.names=F,quote=F,sep='\t')
write.table(outtab,'Diversity_estimates_persite_assay_WV_rawotu.txt',row.names=F,quote=F,sep='\t')


#









###################################
#eddy's working to get to the loop#
###################################
#just starting with generating the output for one site
to_prune<-meta_data_clean %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name=='Tukituki Rvr at Red Br') %>% pull(sample_id)
#prune to a site
FisheDNA_cut<-prune_samples(to_prune,FisheDNA)
sample_data(FisheDNA_cut)
#prune to a site and date
collection_cut<-samp %>% as.matrix %>% as.data.frame() %>% filter(HBRC_Site_Name=='Tukituki Rvr at Red Br') %>% select(CollectionDate) %>% unique() %>% pull(CollectionDate)
collection_cut
to_prune_within<-meta_data_clean %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name=='Tukituki Rvr at Red Br' & CollectionDate=="12/01/22") %>% pull(sample_id)
FisheDNA_cut<-prune_samples(to_prune_within,FisheDNA)
sample_data(FisheDNA_cut)
#this is actually unneccessary seen as Im not transforming here but if we want to switch it to presence absence later on it makes it easier to leave in
phyloseq_object_prune.pa <- microbiome::transform(FisheDNA_cut, 'identity') #raw
otu_table(phyloseq_object_prune.pa)
#      phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa')#pa dataset
rich<-estimate_richness(phyloseq_object_prune.pa, measures=c("Shannon", "Simpson", "InvSimpson"))
rich

colMeans(rich)
date_selected<-"12/01/22"
date_selected<-as.Date(date_selected,format="%d/%m/%y")

#need to drop the winter sampling Im saying summer is 1-nov-31-May 
#this is kinda late but many of the post cyclone samples didnt happen until April May 2023
#but the really weird ones are in Aug etc so hopefully its okay

if (date_selected>='2019-11-01' & date_selected<='2020-05-31'){
  print('2019-2020_summer')
  season<-'2019-2020_summer'
}
if (date_selected>='2019-06-01' & date_selected<='2019-10-31'){
  print('2019_winter')
  season<-'2019_winter'
}
if (date_selected>='2020-11-01' & date_selected<='2021-05-31'){
  print('2020-2021_summer')
  season<-'2020-2021_summer'
}
if (date_selected>='2020-06-01' & date_selected<='2020-10-31'){
  print('2020_winter')
  season<-'2020_winter'
}
if (date_selected>='2021-11-01' & date_selected<='2022-05-31'){
  print('2021-2022_summer')
  season<-'2021-2022_summer'
}
if (date_selected>='2021-06-01' & date_selected<='2021-10-31'){
  print('2021_winter')
  season<-'2021_winter'
}
if (date_selected>='2022-11-01' & date_selected<='2022-05-31'){
  print('2022-2023_summer')
  season<-'2022-2023_summer'
}
if (date_selected>='2022-06-01' & date_selected<='2022-10-31'){
  print('2022_winter')
  season<-'2022_winter'
}
if (date_selected>='2023-11-01' & date_selected<='2024-05-31'){
  print('2023-2024_summer')
  season<-'2023-2024_summer'
}
if (date_selected>='2023-06-01' & date_selected<='2023-10-31'){
  print('2023_winter')
  season<-'2023_winter'
}

#output into a table
rich_mod<-  colMeans(rich) %>% as.data.frame() %>% rename_with(~season) %>% rownames_to_column("Richness_test") %>%  pivot_longer(cols=c(-Richness_test),names_to="Season") %>%  pivot_wider(names_from=c(Richness_test)) %>% mutate(location='Tukituki Rvr at Red Br')
rich_mod