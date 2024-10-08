#better site analysis no all samples

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
library(shinythemes)
# Define UI for application that draws a histogram
fluidPage(theme = shinytheme("united"),
  
  # Application title
  titlePanel("eDNA by site"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("analysis","Select Analysis",choices =c('Site statistics','Site analysis'), selected = 'Site analysis'),
            selectInput("assay", "Select Assay",
                  choices =c("CI","BE","BU","MZ","RV","TP","UM","WV"), selected='CI'),
      selectInput('site_choice',"Select site", unique(meta_data$HBRC_Site_Name)),
      conditionalPanel(
        "input.analysis=='Site analysis'",
        pickerInput("diversity_measure", "Alpha diversity type", c("Observed", "Chao1",  "Shannon", "Simpson", "InvSimpson"),multiple = F)),
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
        pickerInput("subset_data_WV", "Select otu's", c("All","Just freshwater"),multiple = F)),
      selectInput("level", "Select classification level (barchart colours only)",
                  choices =c("Class","Genus"), selected='CI')
    ,width = 3),
   # mainPanel("main panel",
#              column(6,plotOutput(outputId="plot1", width="500px",height="400px"),width = 10)
              
    mainPanel("main panel",
              column(6,plotOutput(outputId="plot1", width="1000px",height="400px"))
    )
    
    # Show a plot of the generated distribution
    #   mainPanel(
    #  plotOutput("distPlot")
    #    plotOutput("plot1")
    #tableOutput("table1")
        )
  )

