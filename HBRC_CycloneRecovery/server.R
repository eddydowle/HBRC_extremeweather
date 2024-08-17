#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library('shiny')
library('DT')
library('data.table')
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
library('leaflet')


setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
meta_data_map<-read.table('Wilderlab_meta_map.txt',sep='\t',quote='',header =T)

records<-read.csv('HBRC_records_Eddy.csv')
#bring in meta data
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
freshwater_invertsfish<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')
#freshwater_invertsfish %>% distinct(Rank) 
family<-freshwater_invertsfish %>% filter(Rank == 'Family')
genus<-freshwater_invertsfish %>% filter(Rank == 'Genus')
phylum<-freshwater_invertsfish %>% filter(Rank == 'Phylum')
order<-freshwater_invertsfish %>% filter(Rank == 'Order')
class<-freshwater_invertsfish %>% filter(Rank == 'Class')

records_eukaryote_freshwater<-records %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)

test<-left_join(records_eukaryote_freshwater,meta_data,by='UID')


function(input, output, session) {
  
 output$secondSelection<-renderUI({selectInput("taxa", "Select taxa", choices =sort((test %>% select(input$level)  %>% unique()) [,1]))})
  
  #probably should be more like summer 2021, summer 2022 or something
  test_filtered<- reactive({
    test_filtered_select<-test %>% filter(eval(parse(text = paste0(input$level))) == input$taxa)
    if (input$integer == 2021) {
      year2021_filtered<-test_filtered_select %>% filter(Year==2021 ) %>% select(HBRC_Site_Name) %>% distinct()
      year2021_all_sites<-test %>% filter(Year==2021 ) %>% select(HBRC_Site_Name) %>% distinct()
      df<-year2021_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2021_filtered$HBRC_Site_Name~"present_2021",TRUE ~ 'Absent'))
      }
         if (input$integer == 2022) {
           year2022_filtered<-test_filtered_select %>% filter(Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
           year2022_all_sites<-test %>% filter(Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2022_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2022_filtered$HBRC_Site_Name~"present_2022",TRUE ~ 'Absent'))
         }
         if (input$integer == 2023) {
           year2023_filtered<-test_filtered_select %>% filter(Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()
           year2023_all_sites<-test %>% filter(Year==2023 ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2023_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2023_filtered$HBRC_Site_Name~"present_2023",TRUE ~ 'Absent'))         }
         if (input$integer == 2024) {
           year2024_filtered<-test_filtered_select %>% filter(Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
           year2024_all_sites<-test %>% filter(Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2024_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2024_filtered$HBRC_Site_Name~"present_2024",TRUE ~ 'Absent'))         }
    if (input$integer == "Precyclone (2021 to March 2023) vs Post cyclone (March 2023 to 2024)") {
      pre_filtered<-test_filtered_select %>% filter(Year==2021|Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
    post_filtered<-test_filtered_select %>% filter(Year==2023|Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
    pre_all_sites<-test %>% filter(Year==2021|Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
    post_all_sites<-test %>% filter(Year==2023|Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
    pre_df<-pre_all_sites %>% mutate(present_pre = case_when(HBRC_Site_Name %in% pre_filtered$HBRC_Site_Name~"present_pre",TRUE ~ 'absent_pre'))
    post_df<-post_all_sites %>% mutate(present_post = case_when(HBRC_Site_Name %in% post_filtered$HBRC_Site_Name~"present_post",TRUE ~ 'absent_post'))
    df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
    df<-df %>% mutate(present = case_when(
      (present_pre == 'absent_pre' & present_post == 'absent_post') ~ 'Absent_pre&post',
      (present_pre == 'present_pre' & present_post == 'present_post') ~ 'Present_pre&post',
      (present_pre == 'present_pre' & present_post == 'absent_post') ~ 'Loss post',
      (present_pre == 'absent_pre' & present_post == 'present_post') ~ 'Gain post', TRUE ~ NA_character_
    ))
    #this results in some NA's as there is some sites (2 that are only sampled post)
    df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')
    }
    #"2022 vs 2023 (Post cyclone)"
    if (input$integer == "2022 vs 2023 (Post cyclone)") {
      filtered_2022<-test_filtered_select %>% filter(Year==2022 ) %>% select(HBRC_Site_Name) %>% distinct()
      filtered_2023<-test_filtered_select %>% filter(Year==2023) %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2022<-test %>% filter(Year==2022) %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2023<-test %>% filter(Year==2023) %>% select(HBRC_Site_Name) %>% distinct()
      pre_df<-all_sites_2022 %>% mutate(present_pre = case_when(HBRC_Site_Name %in% filtered_2022$HBRC_Site_Name~"present_2022",TRUE ~ 'absent_2022'))
      post_df<-all_sites_2023 %>% mutate(present_post = case_when(HBRC_Site_Name %in% filtered_2023$HBRC_Site_Name~"present_2023",TRUE ~ 'absent_2023'))
      df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
      df<-df %>% mutate(present = case_when(
        (present_pre == 'absent_2022' & present_post == 'absent_2023') ~ 'Absent_2022&2023',
        (present_pre == 'present_2022' & present_post == 'present_2023') ~ 'Present_2022&2023',
        (present_pre == 'present_2022' & present_post == 'absent_2023') ~ 'Loss 2023',
        (present_pre == 'absent_2022' & present_post == 'present_2023') ~ 'Gain 2023', TRUE ~ NA_character_
      ))
      #this results in some NA's as there is some sites (2 that are only sampled post)
      df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')%>% filter(!present=='NA')
    }
    if (input$integer == "2023 (Post cyclone) vs 2024") {
      filtered_2024<-test_filtered_select %>% filter(Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
      filtered_2023<-test_filtered_select %>% filter(Year==2023) %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2024<-test %>% filter(Year==2024) %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2023<-test %>% filter(Year==2023) %>% select(HBRC_Site_Name) %>% distinct()
      pre_df<-all_sites_2023 %>% mutate(present_pre = case_when(HBRC_Site_Name %in% filtered_2023$HBRC_Site_Name~"present_2023",TRUE ~ 'absent_2023'))
      post_df<-all_sites_2024 %>% mutate(present_post = case_when(HBRC_Site_Name %in% filtered_2024$HBRC_Site_Name~"present_2024",TRUE ~ 'absent_2024'))
      df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
      df<-df %>% mutate(present = case_when(
        (present_pre == 'absent_2023' & present_post == 'absent_2024') ~ 'Absent_2023&2024',
        (present_pre == 'present_2023' & present_post == 'present_2024') ~ 'Present_2023&2024',
        (present_pre == 'present_2023' & present_post == 'absent_2024') ~ 'Loss 2024',
        (present_pre == 'absent_2023' & present_post == 'present_2024') ~ 'Gain 2024', TRUE ~ NA_character_
      ))
      #this results in some NA's as there is some sites (2 that are only sampled post)
      df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br') %>% filter(!present=='NA')
    }
    
    
  return(left_join(df,meta_data_map,by='HBRC_Site_Name'))
  })
  
  test_filtered_col<- reactive({
  test_filtered_pal<-test_filtered()
  colorFactor(
  palette = 'Dark2',
  domain = test_filtered_pal$present)
})

  output$mymap <- renderLeaflet({ leaflet(test_filtered()) %>% 
    addTiles() %>%
    addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
            popup = ~paste(HBRC_Site_Name),
            color = ~test_filtered_col()(present),
            stroke = FALSE, fillOpacity = 1) %>% 
    addLegend("bottomright", pal = test_filtered_col(), values = ~present,
              title = "Presence/absence",
              opacity = 1) 
  })
  
    }
