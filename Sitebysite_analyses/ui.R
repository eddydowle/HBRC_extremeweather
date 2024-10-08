#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyWidgets)
# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("eDNA by site"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          selectInput("analysis","Select Analysis",choices =c('All locations','By location'), selected = 'All locations'),
         selectInput("assay", "Select Assay",
                            choices =c("CI","BE","BU","MZ","RV","TP","UM","WV"), selected='CI'),
         conditionalPanel("input.analysis=='By location'",
                          pickerInput('site_choice',"Select site", unique(meta_data$HBRC_Site_Name),multiple=F)),
        conditionalPanel(
          "input.assay=='CI'",
          pickerInput("subset_data_CI", "Select otu's", c("All","Just freshwater"),multiple = F))
        ,
        conditionalPanel(
          "input.assay=='RV'",
          pickerInput("subset_data_RV", "Select otu's", c("All","Just freshwater"),multiple = F))
        ,
        conditionalPanel(
          "input.assay=='WV'",
          pickerInput("subset_data_WV", "Select otu's", c("All","Just freshwater"),multiple = F))
        ),
        
        mainPanel("main panel",
                  column(6,plotOutput(outputId="plot1", width="500px",height="400px"))
        )

        # Show a plot of the generated distribution
     #   mainPanel(
          #  plotOutput("distPlot")
      #    plotOutput("plot1")
          #tableOutput("table1")
    #    )
    )
)
