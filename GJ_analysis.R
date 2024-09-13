#using the data processed by GJ:
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


setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/")
samples<-read.csv('OneDrive_1_7-05-2024/HBRC_samples_Eddy.csv')
#records<-read.csv('OneDrive_1_7-05-2024/HBRC_records_Eddy.csv')
meta_data<-read.table('OneDrive_1_7-05-2024/Wilderlab_meta_out.txt',sep='\t',quote='',header =T)

#gj data from BE
gj_analysis<-read.table('GJ_analysis/BE_filter_first_look.txt',sep='\t',header=T,row.names=NULL)

#rename header so that it matches the meta_data:
names(gj_analysis)<-gsub(x = names(gj_analysis), pattern = ".*?_BE_", replacement = "") 

names(gj_analysis)<-gsub(x = names(gj_analysis), pattern = "_.*", replacement = "") 
duplicated(colnames(gj_analysis))
#taxonomy file:
taxonomy<-gj_analysis %>% select(OTU.ID,superkingdom,phylum,class,order,family,genus,species)


#gj data from BE
gj_analysis<-read.table('GJ_analysis/BE_filter_first_look.txt',sep='\t',header=T,row.names=NULL)

#otu 
#need to sum across duplicates
gj_analysis_t<-as.data.frame(t(gj_analysis))
colnames(gj_analysis_t) <- as.character(unlist(gj_analysis_t[2,]))

#make the rownames a columns and stringsub to UID number
gj_analysis_t$UID <- rownames(gj_analysis_t) #%>% gsub(x = names(gj_analysis), pattern = ".*?_BE_", replacement = "")

gj_analysis_t$UID<-gsub(x = gj_analysis_t$UID, pattern = ".*?_BE_", replacement = "")
gj_analysis_t$UID<-gsub(x = gj_analysis_t$UID, pattern = "_.*", replacement = "") 

#drop none otu columns


trans_gj_analysis<-t(gj_analysis)

otu<-gj_analysis %>% select(-superkingdom,-ID,-phylum,-class,-order,-family,-genus,-species,-pident,-qcov,-otherSpeciesID)
