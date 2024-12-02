#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library('shiny')
#library('DT')
#library('data.table')
library('tidyverse')
library('devtools')
library('tidyverse')
#library('reshape2')
library('ggplot2') 
library('leaflet')
library('shinyfullscreen')
library('shinyWidgets')
#app to select species across years
fluidPage(

  ui <- fluidPage(
    tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
     sidebarLayout(
    sidebarPanel(
        selectInput("integer", "Select sampling year or Comparison",
                    choices =c("2020-2021","2021-2022","2022-2023","2023-2024","Pre-cyclone (2019 to March 2023) vs Post-cyclone (March 2023 to 2024)","2021-2022 Pre-cyclone vs 2023 (March-June) Post-cyclone","2023 (March-June) vs 2023-2024 (July 2023-April 2024)"), selected='2023-2024'),
        selectInput("level", "Select taxonomic level", choices = c("Family", "Genus", "Species"), selected='Family'),
    selectInput("filter_group", "Select taxa filter", choices = c("All", "Insecta", "Freshwater fish"), selected='All'),  
     uiOutput("secondSelection"),
 #  selectizeInput('secondSelection','Select Taxa', choices = NULL),
    selectInput("map_background",'Change background', choices = c('Satellite','Stylised','Grey-scale'),selected='Terrain')),
#   selectInput("taxa", "Select taxa", choices =sort((test %>% select(Family)  %>% unique()) [,1]))
#   selectInput("taxa", "Select taxa", choices = (test %>% select(Family)  %>% unique())))
   
   #,
   
#needs to be reactive to level
    mainPanel(leafletOutput("mymap"),DT::dataTableOutput('table'),
              fullscreen_those(items = list("mymap"))
              )
#mainPanel(tableOutput("table1"))
     )
    )
)
