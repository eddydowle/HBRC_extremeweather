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

  ui <- fluidPage(
    leafletOutput("mymap")
  
    # Application title
  #  titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
  #  sidebarLayout(
   #     sidebarPanel(
    #        sliderInput("bins",
     #                   "Number of bins:",
      #                  min = 1,
       #                 max = 50,
        #                value = 30)
        #),

        # Show a plot of the generated distribution
        #mainPanel(
         #   plotOutput("distPlot")
        #)

#        selectInput("Year_select", label = h3("Select year"), 
 #                   choices = c("All years", "2021" ,  "2022" ,"2023","2024"), 
  #                  selected = "All regions",multiple = TRUE),        

    )
)
