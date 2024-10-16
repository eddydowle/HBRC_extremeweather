#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

#app to select species across years
fluidPage(

  ui <- fluidPage(
    tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
     sidebarLayout(
    sidebarPanel(
        selectInput("integer", "Select sampling year or Comparison",
                    choices =c("2020-2021","2021-2022","2022-2023","2023-2024","Pre-cyclone (2019 to March 2023) vs Post-cyclone (March 2023 to 2024)","2021-2022 Pre-cyclone vs 2022-2023 Post-cyclone","2022-2023 (Post-cyclone only) vs 2023-2024"), selected='2023-2024'),
        
        selectInput("level", "Select taxonomic level", choices = c("Family", "Genus", "Species"), selected='Family'),        uiOutput("secondSelection")),
#   selectInput("taxa", "Select taxa", choices =sort((test %>% select(Family)  %>% unique()) [,1]))
#   selectInput("taxa", "Select taxa", choices = (test %>% select(Family)  %>% unique())))
   
   #,
   
#needs to be reactive to level
    mainPanel(leafletOutput("mymap"))
#mainPanel(tableOutput("table1"))
     )
    )
)
