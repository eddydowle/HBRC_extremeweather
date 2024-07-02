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

df<-df %>% filter(!HBRC_Site_Name=='Wharerangi - Dan and Lindsay Bates') %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')

df %>% filter(!HBRC_Site_Name=='Mangaone Rv at Dartmoor Br')
