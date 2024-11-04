setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
records<-read.csv('HBRC_records_Eddy.csv')
#just cutting down to make it lighter for testing
#records<-read.csv('HBRC_records_Eddy_random50k.csv')
#meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
meta_data<-read.table('Wilderlab_meta_out_clean.txt',sep='\t',quote='',header =T)
otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
#write.csv(otu_table,'HBRC_records_Eddy_random50k_otuapp.csv',row.names=F)
#write.csv(taxa_table,'HBRC_records_Eddy_random50k_taxaapp.csv',row.names=F)
taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
#otu_table<-read.csv('HBRC_records_Eddy_random50k_otuapp.csv')
#taxa_table<-read.csv('HBRC_records_Eddy_random50k_taxaapp.csv')
freshwater_taxa<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')

meta_data_2<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))
meta_data_cl <-meta_data_2 %>% filter(UID %in% otu_table_subset()$UID)
row.names(meta_data_cl) <- meta_data_cl$UID
meta_data_cl[1] <- NULL

