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


#setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
#records<-read.csv('HBRC_records_Eddy.csv')
#bring in meta data
#meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
#freshwater_invertsfish<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')
#freshwater_invertsfish %>% distinct(Rank) 
#family<-freshwater_invertsfish %>% filter(Rank == 'Family')
#genus<-freshwater_invertsfish %>% filter(Rank == 'Genus')
#phylum<-freshwater_invertsfish %>% filter(Rank == 'Phylum')
#order<-freshwater_invertsfish %>% filter(Rank == 'Order')
#class<-freshwater_invertsfish %>% filter(Rank == 'Class')

#records_eukaryote_freshwater<-records %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)

#test<-left_join(records_eukaryote_freshwater,meta_data,by='UID')


# Define server logic required to draw a histogram
function(input, output, session) {
  
 # output$secondSelection<-renderUI({selectInput("taxa", "Select taxa", choices =sort((test %>% select(input$level)  %>% unique()) [,1]))})
  
  #probably should be more like summer 2021, summer 2022 or something
  test_filtered<- reactive({
#    test_filtered_select<-test %>% filter(.data[[input$level]] == input$secondSelection)
    test_filtered_select<-test %>% filter(Family == input$taxa)
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
           year2024_filtered<-test_filtered_select %>% filter(Year==2021 ) %>% select(HBRC_Site_Name) %>% distinct()
           year2024_all_sites<-test %>% filter(Year==2024 ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2024_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2024_filtered$HBRC_Site_Name~"present_2024",TRUE ~ 'Absent'))         }
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
