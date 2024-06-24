test_filtered_select<-test %>% filter(Family=='Nesameletidae')
range_test<-2021
if (range_test == 2021) {
  year2021_filtered<-test_filtered_select %>% filter(Year==2021 ) %>% select(HBRC_Site_Name) %>% distinct()
  year2021_all_sites<-test %>% filter(Year==2021 ) %>% select(HBRC_Site_Name) %>% distinct()
  df_test<-year2021_all_sites %>% mutate(present = case_when(HBRC_Site_Name %in% year2021_filtered$HBRC_Site_Name~"present_2021",TRUE ~ 'Absent'))
}

x<-left_join(df_test,meta_data_map,by='HBRC_Site_Name')

y<-colorFactor(
  palette = 'Dark2',
  domain = x$present)

leaflet(x) %>% 
  addTiles() %>%
  addCircleMarkers(~Longitude_hbrc, ~Latitude_hbrc,
                   popup = ~paste(HBRC_Site_Name),
 #                  color = ~ y(present),
                   stroke = FALSE, fillOpacity = 1) #%>% 
 # addLegend("bottomright", pal = y, values = ~present,
      #      title = "Presence/absence",
      #      opacity = 1)
