#eddy dowle
#march 2024

#trying to align the sites across from the file Shaun gave me and the one from the RC

setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/")
library('tidyverse')
library('readxl')

#hbrc table
hbrc_table<-read_excel('eDNAsitenames_HBRC_table.xlsx', sheet = "Sheet2")
print(hbrc_table$E,digits=12)
#have converted E/N to lat long
#have used NZ transverse mercator projection to NZ geodetic datum 2000 (version 20180701)
#wasnt 100% on input but it seems close
#have added a 'b' to sample sites that are collected via a passive filter. 

#there is a few duplicate rows Im not sure why but just remove all rows that are complete dups
hbrc_table_nodups <- hbrc_table[!duplicated(hbrc_table),]

#wilder lab
#wilder_table<-read.csv('HBRC_samples_2024_wilder.csv',header=T, row.names=NULL, sep=',',quote='')
wilder_table<-read_excel('HBRC_samples_2024_wilder.xlsx', sheet = "HBRC_samples_2024_wilder")
#wilder_table<-read_excel('Wilder_simple.xlsx', sheet = "Sheet1")

#just do a quick map to see where we are:
library(leaflet)
library(htmltools)
library(tidyverse)
library(RColorBrewer)

#load the metadata file
#this must include the longitude and latitude data + whatever data you wish to display

#leaflet(hbrc_table_nodups) %>% addTiles() %>%
#  addMarkers(~Longitude, ~Latitude, popup = ~paste(Wilderlab_Site_Names,HBRC_Site_Name, HBRC_Site_ID, sep = "<br>"))
#a few are a bit of the river, not sure if that is leaflet or a GPS problem (either I have used the wrong NE version going in or spot was not taken in right place)
#for example 3731 appears to be slightly off

#okay so relabel wilderlab as either the site or 

wilder_meta<-left_join(wilder_table, hbrc_table_nodups, by = c("ClientSampleID" = "Wilderlab_Site_Names"))

#which sites dont have a match n the hbrc table:
no_match<-wilder_meta %>% filter(is.na(HBRC_Site_Name))
#sweet now everything has a match between the two

#Okay dokey
#we want:
#HB site ID, HB site name, dates sampled, No replicates at each sample date, lat longs, 
test<-wilder_meta %>% group_by(HBRC_Site_Name) %>% summarise(total_replicates = n(), HBRC_Site_ID = toString(unique(HBRC_Site_ID)), samplingdates = toString(unique(CollectionDate)),  n_uniq_samplingdates = n_distinct(unique(CollectionDate)),Latitude_wilder = toString(unique(Latitude.x)), Longitude_wilder = toString(unique(Longitude.x)), Latitude_hbrc = toString(unique(Latitude.y)), Longitude_hbrc = toString(unique(Longitude.y)), cyclone = toString(unique(Cyclone)),.groups = 'drop') %>% drop_na(.,HBRC_Site_ID)


#okay so sites that we have pre and post cyclone data and appropriate metadata from HBRC

test_prepost<-test %>% filter(cyclone=='Pre, Post'|cyclone =='Post, Pre')
test_pre<-test %>% filter(cyclone=='Pre')
test_post<-test %>% filter(cyclone=='Post')

#what is the lowest number of replicates per sample point
lowest_sample_number<-wilder_meta %>% group_by(HBRC_Site_ID,CollectionDate) %>%  summarise(lowest_rep_n=n(),.groups='drop') %>% select(-CollectionDate)  %>% group_by(HBRC_Site_ID) %>%  filter(rank(lowest_rep_n, ties.method= "first") == 1) %>% drop_na(.,HBRC_Site_ID)

#lowest value for each site

test_prepost<-left_join(test_prepost,lowest_sample_number,by='HBRC_Site_ID')
test_pre<-left_join(test_pre,lowest_sample_number,by='HBRC_Site_ID')
test_post<-left_join(test_post,lowest_sample_number,by='HBRC_Site_ID')

#get rid of that weird fishing site in prepost
#test_prepost <- test_prepost %>% filter(HBRC_Site_ID!='Fishing')

#okay lets output this:
#write.table(test_prepost,"HBRC_sampling_prepost_metadata_reconciled.txt",row.names=F,quote=F,sep="\t")

#write.table(test_pre,"HBRC_sampling_pre_metadata_reconciled.txt",row.names=F,quote=F,sep="\t")

#write.table(test_post,"HBRC_sampling_post_metadata_reconciled.txt",row.names=F,quote=F,sep="\t")

#all right lets deal with those samples that we dont have meta data for. 
test_nometa<-wilder_meta %>% filter(is.na(HBRC_Site_ID) | HBRC_Site_ID=='Fishing') %>% select(JobID,ClientSampleID,CollectionDate,Latitude.x,Longitude.x,Cyclone,HBRC_Site_Name,HBRC_Site_ID,E,N,Latitude.y,Longitude.y) %>% unique()
unique(test$ClientSampleID)

#write.table(test_nometa,"HBRC_sampling_notreconciled.txt",row.names=F,quote=F,sep='\t')

#playing with the table of all sites that HBRC gave us
#turns out this table is not one we should consider dont know why
#eDNA_plan<-read_excel('eDNA_plan_OG_HBRC.xlsx', sheet = "Sheet2")

#test_allsites<-left_join(test,lowest_sample_number,by='HBRC_Site_ID')

#not sure whether best to join by sample name or site id
#try site name first as some have no site ID
#test_allsites_eDNAplan<-full_join(eDNA_plan,test_allsites,by='HBRC_sitename')
#write.table(test_allsites_eDNAplan,'eDNA_plan_merge_wilder.txt',row.names=F, quote=F,sep='\t')


#alright new table from wilderlab
#this includes 2024 samples sequenced and those from niwa
wilder_table<-read_csv('HBRC_samples_2024_HB_NIWA_MfE_v3.csv')
#deployment duration indicates a passive site I think so drop those:

wilder_syring<-wilder_table %>% filter(is.na(DeploymentDuration))

#lets drop the niwa samples for now as we are not sure what the deal is with them:
wilder_syring_HB_MfE<-wilder_table %>% filter(!AccountNumber=='NIW001')
#to not kick out NIWA
wilder_syring_HB_MfE<-wilder_table
#okay align with the one from HB

#hbrc table
hbrc_table<-read_excel('eDNAsitenames_HBRC_table.xlsx', sheet = "Sheet2")
print(hbrc_table$E,digits=12)
#have converted E/N to lat long
#have used NZ transverse mercator projection to NZ geodetic datum 2000 (version 20180701)
#wasnt 100% on input but it seems close
#have added a 'b' to sample sites that are collected via a passive filter. 

#there is a few duplicate rows Im not sure why but just remove all rows that are complete dups
hbrc_table_nodups <- hbrc_table[!duplicated(hbrc_table),]

wilder_meta<-left_join(wilder_syring_HB_MfE, hbrc_table_nodups, by = c("ClientSampleID" = "Wilderlab_Site_Names"))
#ffs there is more naming issues. 
no_match<-wilder_meta %>% filter(is.na(HBRC_Site_Name))
#Deep Stream at Old Taupo Coach Rd
#Poporangi Stream at Wakarara Rd
#Tarere Stream at Devils Elbow
#Triplex Creek at Triplex hut (NIWA)
#not sure what is happening to those sites just remove for now

wilder_meta_clean<-wilder_meta %>% filter(!is.na(HBRC_Site_Name))

#alright clean it up
#the new wilderlab table has 60 samples that dont have a date (not sure what is going on there)
test<-wilder_meta_clean %>% group_by(HBRC_Site_Name) %>% summarise(total_replicates = n(), HBRC_Site_ID = toString(unique(HBRC_Site_ID)), samplingdates = toString(unique(CollectionDate)),  n_uniq_samplingdates = n_distinct(unique(CollectionDate)),Latitude_wilder = toString(unique(Latitude.x)), Longitude_wilder = toString(unique(Longitude.x)), Latitude_hbrc = toString(unique(Latitude.y)), Longitude_hbrc = toString(unique(Longitude.y)), cyclone = toString(unique(Cyclone)),.groups = 'drop') %>% drop_na(.,HBRC_Site_ID)

#putting in UID for shaun to pull samples out by
test<-wilder_meta_clean %>% group_by(HBRC_Site_Name) %>% summarise(total_replicates = n(), HBRC_Site_ID = toString(unique(HBRC_Site_ID)), samplingdates = toString(unique(CollectionDate)),  n_uniq_samplingdates = n_distinct(unique(CollectionDate)),Latitude_wilder = toString(unique(Latitude.x)), Longitude_wilder = toString(unique(Longitude.x)), Latitude_hbrc = toString(unique(Latitude.y)), Longitude_hbrc = toString(unique(Longitude.y)), UIDs = toString(UID), UIDs_count = n_distinct(unique(UID)),cyclone = toString(unique(Cyclone)),.groups = 'drop') %>% drop_na(.,HBRC_Site_ID)
#double checking we have the same number of UIDs a we do total replicates:
identical(test[['total_replicates']],test[['UIDs_count']])
#true so should be able to send list of UIDs to shaun

#okay so sites that we have pre and post cyclone data and appropriate metadata from HBRC
#those nodates are mucking up the filters to switched to a grepl
#test_prepost<-test %>% filter(cyclone=='Pre, Post'|cyclone =='Post, Pre')
test_prepost<-test %>% filter(grepl ("Pre",cyclone) & grepl("Post", cyclone))
#test_pre<-test %>% filter(cyclone=='Pre')
test_pre<-test %>% filter(grepl ("Pre",cyclone) & !grepl("Post", cyclone))
#test_post<-test %>% filter(cyclone=='Post')
test_post<-test %>% filter(!grepl ("Pre",cyclone) & grepl("Post", cyclone))

#what is the lowest number of replicates per sample point
lowest_sample_number<-wilder_meta_clean %>% group_by(HBRC_Site_Name,CollectionDate) %>%  summarise(lowest_rep_n=n(),.groups='drop') %>% select(-CollectionDate)  %>% group_by(HBRC_Site_Name) %>%  filter(rank(lowest_rep_n, ties.method= "first") == 1) %>% drop_na(.,HBRC_Site_Name)

#lowest value for each site

test_prepost<-left_join(test_prepost,lowest_sample_number,by='HBRC_Site_Name')
test_pre<-left_join(test_pre,lowest_sample_number,by='HBRC_Site_Name')
test_post<-left_join(test_post,lowest_sample_number,by='HBRC_Site_Name')

#get rid of that weird fishing site in prepost
#test_prepost <- test_prepost %>% filter(HBRC_Site_ID!='Fishing')

#okay lets output this:
#write.table(test_prepost,"HBRC_sampling_prepost_metadata_reconciled_new_niwa_UID.txt",row.names=F,quote=F,sep="\t")

#write.table(test_pre,"HBRC_sampling_pre_metadata_reconciled_new_niwa_UID.txt",row.names=F,quote=F,sep="\t")

#write.table(test_post,"HBRC_sampling_post_metadata_reconciled_new_niwa_UID.txt",row.names=F,quote=F,sep="\t")

#create table for meta data analyses
#.x on lat long is wilder
wilder_meta_clean_out<-wilder_meta_clean %>% select(UID, CollectionDate,Latitude.x,Longitude.x,HBRC_Site_Name,HBRC_Site_ID,Cyclone)

write.table(wilder_meta_clean_out,'ASV_files/OneDrive_1_7-05-2024/Wilderlab_meta_out_v2.txt',row.names=F,quote=F,sep='\t')
