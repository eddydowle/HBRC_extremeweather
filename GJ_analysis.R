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
records<-read.csv('OneDrive_1_7-05-2024/HBRC_records_Eddy.csv')
meta_data<-read.table('OneDrive_1_7-05-2024/Wilderlab_meta_out.txt',sep='\t',quote='',header =T)





