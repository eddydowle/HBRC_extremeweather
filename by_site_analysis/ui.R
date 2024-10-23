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
      selectInput('site_choice',"Select site", sort(unique(meta_data$HBRC_Site_Name))),
      conditionalPanel(
        "input.analysis=='Site analysis'",
        pickerInput("diversity_measure", "Alpha diversity type", c("Observed", "Chao1",  "Shannon", "Simpson", "InvSimpson"),multiple = F)),
      conditionalPanel(
        "input.assay=='CI'",
        pickerInput("subset_data_CI", "Select otu's", c("All","Just freshwater"),multiple = F,selected='All'))
      ,
      conditionalPanel(
        "input.assay=='RV'",
        pickerInput("subset_data_RV", "Select otu's", c("All","Just freshwater"),multiple = F,selected='All'))
      ,
      conditionalPanel(
        "input.assay=='WV'",
        pickerInput("subset_data_WV", "Select otu's", c("All","Just freshwater"),multiple = F,selected='All')),
      selectInput("level", "Select classification level (barchart colours only)",
                  choices =c("Class","Genus"), selected='CI')
    ,width = 3),

    #mainPanel("main panel",
   #           column(6,plotOutput(outputId="plot1", width="1000px",height="400px"))
  #  )))
#mainPanel("main panel",
#          column(6,plotOutput(outputId="plot1", width="1000px",height="400px"),uiOutput("lm1"))
#)))
#make it conditional on input type
#need seperate plot labels to make conditional work (R doesnt like two thing with the same name)
mainPanel(
  conditionalPanel(
    "input.analysis=='Site analysis'",
    column(6,plotOutput(outputId="plot1", width="1000px",height="400px"),uiOutput("lm1"))
  ),
  conditionalPanel(
    "input.analysis=='Site statistics'",
    column(6,plotOutput(outputId="plot2", width="1000px",height="400px"))
  )
   )))

#that means I can probably modify how I display the plots ~ but maybe not as innext will have a dependency in it to skip if its not enough samples
