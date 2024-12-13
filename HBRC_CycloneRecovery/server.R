#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library('shinyfullscreen')
library('shiny')
#library('DT')
#library('data.table')
library('tidyverse')
#library('readxl')
library('devtools')
library('tidyverse')
#library('viridis')
#library('indicspecies') 
#library('vegan')
#library('glue')
#library('reshape2')
#library('phyloseq')
#library('iNEXT.3D') 
library('ggplot2') 
#library('microbiome')
library('leaflet')


#setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
meta_data_map<-read.table('Wilderlab_meta_map2.txt',sep='\t',quote='',header =T)

#records<-read.csv('HBRC_records_Eddy.csv')
#bring in meta data
#meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
#modified with additional column for season
#meta_data<-read.table('Wilderlab_meta_out_season.txt',sep='\t',header =T)

#freshwater_invertsfish<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')
#freshwater_invertsfish %>% distinct(Rank) 
#family<-freshwater_invertsfish %>% filter(Rank == 'Family')
#genus<-freshwater_invertsfish %>% filter(Rank == 'Genus')
#phylum<-freshwater_invertsfish %>% filter(Rank == 'Phylum')
#order<-freshwater_invertsfish %>% filter(Rank == 'Order')
#class<-freshwater_invertsfish %>% filter(Rank == 'Class')

#records_eukaryote_freshwater<-records %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)

#test<-left_join(records_eukaryote_freshwater,meta_data,by='UID')
#write.table(test,'HBRC_records_Eddy_freshwater_season.txt',sep='\t',quote=F,row.names = F)
test<-read.table('HBRC_records_Eddy_freshwater_season.txt',header=T,row.names=NULL,sep='\t',quote='')
#need to filter it for the 8 common assays
#CI: Mostly insects","BE: General Eukaryote","BU: General Eukaryote","MZ: Vascular Plants","RV: Vertebrates","TP: Vascular Plants","UM: Microbe","WV: Vertebrates"
#test2<-test %>% filter(PrimerSet=='CI'|PrimerSet=='BE'|PrimerSet=='BU'|PrimerSet=='MZ'|PrimerSet=='RV'|PrimerSet=='TP'|PrimerSet=='UM'|PrimerSet=='WV')
#write.table(test2,'HBRC_records_Eddy_freshwater_season.txt',sep='\t',quote=F,row.names = F)

function(input, output, session) {

#updateSelectizeInput(session, 'Select Taxa', choices = sort((test %>% select(input$level)  %>% unique()) [,1]), server = TRUE)

    
 output$secondSelection<-renderUI({
   if (input$filter_group=='All'){
     test_filtered<-test
   }
   if (input$filter_group=='Insecta'){
     test_filtered<-test %>% filter(Class=='Insecta')
   }
   if (input$filter_group=='Freshwater fish'){
     test_filtered<-test %>% filter(Class %in% c('Actinopteri','Hyperoartia'))
   }
#   selectInput("taxa", "Select taxa", choices =sort((test_filtered %>% select(input$level)  %>% unique()) [,1]))
  pickerInput("taxa", "Select taxa", choices =sort((test_filtered %>% select(input$level)  %>% unique()) [,1]),multiple=F)
    })

  #probably should be more like summer 2021, summer 2022 or something
  test_filtered<- reactive({
    test_filtered_select<-test %>% filter(eval(parse(text = paste0(input$level))) == input$taxa)
  #  print('here')
  #  dropping 2019-2020 as there is only one site
    if (input$integer == '2023-2024') {
      year2021_filtered<-test_filtered_select %>% filter(season=='2023-2024' ) %>% select(HBRC_Site_Name) %>% distinct()
      year2021_all_sites<-test %>% filter(season=='2023-2024' ) %>% select(HBRC_Site_Name) %>% distinct()
      df<-year2021_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2021_filtered$HBRC_Site_Name~"Present 2023-2024",TRUE ~ 'Absent'))
     # print(df)
      }
         if (input$integer == '2020-2021') {
           year2022_filtered<-test_filtered_select %>% filter(season=='2020-2021' ) %>% select(HBRC_Site_Name) %>% distinct()
           year2022_all_sites<-test %>% filter(season=='2020-2021' ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2022_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2022_filtered$HBRC_Site_Name~"Present 2020-2021",TRUE ~ 'Absent'))
         }
         if (input$integer == '2021-2022') {
           year2023_filtered<-test_filtered_select %>% filter(season=='2021-2022' ) %>% select(HBRC_Site_Name) %>% distinct()
           year2023_all_sites<-test %>% filter(season=='2021-2022' ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2023_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2023_filtered$HBRC_Site_Name~"Present 2021-2022",TRUE ~ 'Absent'))         }
         if (input$integer == '2022-2023') {
           year2024_filtered<-test_filtered_select %>% filter(season=='2022-2023' ) %>% select(HBRC_Site_Name) %>% distinct()
           year2024_all_sites<-test %>% filter(season=='2022-2023' ) %>% select(HBRC_Site_Name) %>% distinct()
           df<-year2024_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2024_filtered$HBRC_Site_Name~"Present 2022-2023",TRUE ~ 'Absent'))         
           }
    if (input$integer == "Pre-cyclone (2019 to March 2023) vs Post-cyclone (March 2023 to 2024)") {
      pre_filtered<-test_filtered_select %>% filter(Cyclone=='Pre' ) %>% select(HBRC_Site_Name) %>% distinct()
    post_filtered<-test_filtered_select %>% filter(Cyclone=='Post' ) %>% select(HBRC_Site_Name) %>% distinct()
    pre_all_sites<-test %>% filter(Cyclone=='Pre') %>% select(HBRC_Site_Name) %>% distinct()
    post_all_sites<-test %>% filter(Cyclone=='Post') %>% select(HBRC_Site_Name) %>% distinct()
    pre_df<-pre_all_sites %>% mutate(present_pre = case_when(HBRC_Site_Name %in% pre_filtered$HBRC_Site_Name~"present_pre",TRUE ~ 'absent_pre'))
    post_df<-post_all_sites %>% mutate(present_post = case_when(HBRC_Site_Name %in% post_filtered$HBRC_Site_Name~"present_post",TRUE ~ 'absent_post'))
    df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
    df<-df %>% mutate(present = case_when(
      (present_pre == 'absent_pre' & present_post == 'absent_post') ~ 'Absent pre & post',
      (present_pre == 'present_pre' & present_post == 'present_post') ~ 'Present pre & post',
      (present_pre == 'present_pre' & present_post == 'absent_post') ~ 'Loss post-cyclone',
      (present_pre == 'absent_pre' & present_post == 'present_post') ~ 'Gain post-cyclone', TRUE ~ NA_character_
    ))
    #this results in some NA's as there is some sites (2 that are only sampled post)
    df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')
    }
   # 2021-2022 Pre-cyclone vs 2023 (March-June) Post-cyclone
    if (input$integer == "2021-2022 Pre-cyclone vs 2023 (March-June) Post-cyclone") {
      filtered_2022<-test_filtered_select %>% filter(season=='2021-2022' ) %>% select(HBRC_Site_Name) %>% distinct()
      filtered_2023<-test_filtered_select %>% filter(season=='2022-2023' ) %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2022<-test %>% filter(season=='2021-2022') %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2023<-test %>% filter(season=='2022-2023') %>% select(HBRC_Site_Name) %>% distinct()
      pre_df<-all_sites_2022 %>% mutate(present_pre = case_when(HBRC_Site_Name %in% filtered_2022$HBRC_Site_Name~"present_pre-cyclone",TRUE ~ 'absent_pre-cyclone'))
      post_df<-all_sites_2023 %>% mutate(present_post = case_when(HBRC_Site_Name %in% filtered_2023$HBRC_Site_Name~"present_post-cyclone",TRUE ~ 'absent_post-cyclone'))
      df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
      df<-df %>% mutate(present = case_when(
        (present_pre == 'absent_pre-cyclone' & present_post == 'absent_post-cyclone') ~ 'Absent pre & post',
        (present_pre == 'present_pre-cyclone' & present_post == 'present_post-cyclone') ~ 'Present pre & post',
        (present_pre == 'present_pre-cyclone' & present_post == 'absent_post-cyclone') ~ 'Loss post-cyclone',
        (present_pre == 'absent_pre-cyclone' & present_post == 'present_post-cyclone') ~ 'Gain post-cyclone', TRUE ~ NA_character_
      ))
      #this results in some NA's as there is some sites (2 that are only sampled post)
      df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')%>% filter(!present=='NA')
    }
    #choices =c("2020-2021","2021-2022","2022-2023","2023-2024","Pre-cyclone (2019 to March 2023) vs Post-cyclone (March 2023 to 2024)","2021-2022 Pre-cyclone vs 2022-2023 Post-cyclone","2022-2023 (Post-cyclone only) vs 2023-2024"), selected='2023-2024'),
    #2023 (March-June) vs 2023-2024 (July 2023-April 2024)
    if (input$integer == "2023 (March-June) vs 2023-2024 (July 2023-April 2024)") {
      filtered_2024<-test_filtered_select %>% filter(season=='2022-2023' ) %>% select(HBRC_Site_Name) %>% distinct()
      filtered_2023<-test_filtered_select %>% filter(season=='2023-2024') %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2024<-test %>% filter(season=='2023-2024') %>% select(HBRC_Site_Name) %>% distinct()
      all_sites_2023<-test %>% filter(season=='2022-2023') %>% select(HBRC_Site_Name) %>% distinct()
      pre_df<-all_sites_2023 %>% mutate(present_pre = case_when(HBRC_Site_Name %in% filtered_2023$HBRC_Site_Name~"present_2023",TRUE ~ 'absent_2023'))
      post_df<-all_sites_2024 %>% mutate(present_post = case_when(HBRC_Site_Name %in% filtered_2024$HBRC_Site_Name~"present_2024",TRUE ~ 'absent_2024'))
      df<-full_join(pre_df,post_df,by='HBRC_Site_Name')
      df<-df %>% mutate(present = case_when(
        (present_pre == 'absent_2023' & present_post == 'absent_2024') ~ 'Absent 2023 & 2024',
        (present_pre == 'present_2023' & present_post == 'present_2024') ~ 'Present 2023 & 2024',
        (present_pre == 'present_2023' & present_post == 'absent_2024') ~ 'Loss 2024',
        (present_pre == 'absent_2023' & present_post == 'present_2024') ~ 'Gain 2024', TRUE ~ NA_character_
      ))
      #this results in some NA's as there is some sites (2 that are only sampled post)
      df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br') %>% filter(!present=='NA')
    }
    
   # print(df)
    #set colours for map
    df<-df %>% mutate(colour=case_when(present == 'Absent' ~ 'grey',present == 'Present 2023-2024' ~ 'steelblue', present=='Present 2020-2021' ~'steelblue', present=='Present 2021-2022'~ 'steelblue', present=='Present 2022-2023'~'steelblue',present=='Absent pre & post'~'grey',present=='Present pre & post' ~'steelblue',present=='Loss post-cyclone'~'tomato',present=='Gain post-cyclone'~"limegreen",present=='Absent 2023 & 2024'~'grey',present=='Present 2023 & 2024'~'steelblue',present=='Loss 2024'~'tomato',present=='Gain 2024'~'limegreen'))
  # print(head(df))
  return(left_join(df,meta_data_map,by='HBRC_Site_Name'))
  })
  
  test_filtered_col<- reactive({
  test_filtered_pal<-test_filtered()
  #colorFactor(
  #palette = 'Dark2',
  #domain = test_filtered_pal$present)
 # print(unique(test_filtered_pal$colour))
 # print(unique(as.factor(test_filtered_pal$present)))
  test_filtered_pal$present<-factor(test_filtered_pal$present, levels=unique(test_filtered_pal$present))
  colorFactor(palette=unique(test_filtered_pal$colour),
              domain=test_filtered_pal$present
  )
})

  output$mymap <- renderLeaflet({ 
    validate(
      need(input$taxa != "", "No taxa selected") # display custom message in need
    )
    if (input$map_background == "Satellite") {
      mp<-leaflet(test_filtered()) %>% 
      addProviderTiles(
        providers$Esri.WorldImagery,
        options = providerTileOptions(opacity = 0.90)
      ) %>%
    addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
            popup = ~paste(HBRC_Site_Name),
            color = ~test_filtered_col()(present),
            radius=7,
 #           color = colour,
            stroke = FALSE, fillOpacity = 1) %>% 
    addLegend("bottomright", pal = test_filtered_col(), values = ~present,
              title = "Presence/absence",
              opacity = 1)
    }
    if (input$map_background == "Stylised") {
      mp<-leaflet(test_filtered()) %>% 
        addTiles() %>%
        addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                         popup = ~paste(HBRC_Site_Name),
                         color = ~test_filtered_col()(present),
                         radius=7,
                         #           color = colour,
                         stroke = FALSE, fillOpacity = 1) %>% 
        addLegend("bottomright", pal = test_filtered_col(), values = ~present,
                  title = "Presence/absence",
                  opacity = 1)
    }
    if (input$map_background == "Grey-scale") {
      mp<-leaflet(test_filtered()) %>% 
        addProviderTiles(
          providers$Esri.WorldGrayCanvas,
          options = providerTileOptions(opacity = 0.90)
        ) %>%
        addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                         popup = ~paste(HBRC_Site_Name),
                         color = ~test_filtered_col()(present),
                         radius=7,
                         #           color = colour,
                         stroke = FALSE, fillOpacity = 1) %>% 
        addLegend("bottomright", pal = test_filtered_col(), values = ~present,
                  title = "Presence/absence",
                  opacity = 1)
    }
    return(mp)
  })
  
  output$table <- DT::renderDataTable({
    tab_out<-test_filtered() %>% rename(Status = present) %>% select(HBRC_Site_Name,Status,HBRC_Site_ID) %>% arrange(HBRC_Site_Name)
    tab_out %>% 
      DT::datatable()
  })
  
    }
