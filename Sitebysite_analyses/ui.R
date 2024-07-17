#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
      ,      
        selectInput("assay", "Select Assay",
                            choices =c("CI","BE","BU","MZ","RV","TP","UM","WV"), selected='CI')),
        

        # Show a plot of the generated distribution
        mainPanel(
          #  plotOutput("distPlot")
          plotOutput("plot1")
          #tableOutput("table1")
        )
    )
)
