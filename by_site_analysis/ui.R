#better site analysis no all samples

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library('shiny')
library('shinyWidgets')
library('shinythemes')
library('rsconnect')
library('shiny')
library('microViz')
library('RColorBrewer')
library('tidyverse')
#library('readxl')
library('devtools')
#library('tidyverse')
#library('viridis')
#library('indicspecies') 
library('vegan')
#library('glue')
library('reshape2')
library('phyloseq')
library('iNEXT.3D') 
library('ggplot2') 
library('microbiome')
library('gridExtra')
library('ranacapa')
library('stargazer')

meta_data<-read.table('Wilderlab_meta_out_clean.txt',sep='\t',quote='',header =T)

#otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
#write.table(otu_table,'HBRC_records_Eddy_otu_table.txt',sep='\t',quote=F,row.names = F)

#taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
#write.table(taxa_table,'HBRC_records_Eddy_taxa_table.txt',sep='\t',quote=F,row.names = F)

#otu_table<-read.csv('HBRC_records_Eddy_otu_table.txt',header=T,row.names=NULL,sep='\t',quote='')
#otu_table<-read.csv('HBRC_records_Eddy_otu_table_cut.txt',header=T,row.names=NULL,sep='\t',quote='')
taxa_table<-read.csv('HBRC_records_Eddy_taxa_table.txt',header=T,row.names=NULL,sep='\t',quote='')
freshwater_taxa<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')

# Define UI for application that draws a histogram
fluidPage(theme = shinytheme("united"),
          
          # Application title
          titlePanel("eDNA by site"),
          
          # Sidebar with a slider input for number of bins
          sidebarLayout(
            sidebarPanel(
              selectInput("analysis","Select Analysis",choices =c('Site statistics','Site analysis'), selected = 'Site analysis'),
              selectInput("assay", "Select Assay",
                          choices =c("CI: Mostly insects","BE: General Eukaryote","BU: General Eukaryote","MZ: Vascular Plants","RV: Vertebrates","TP: Vascular Plants","UM: Microbe","WV: Vertebrates"), selected='RV: Vertebrates'),
              selectInput('site_choice',"Select site", sort(unique(meta_data$HBRC_Site_Name))),
              conditionalPanel(
                "input.analysis=='Site analysis'",
                pickerInput("diversity_measure", "Alpha diversity type", c("Observed", "Chao1",  "Shannon", "Simpson", "InvSimpson"),multiple = F)),
              conditionalPanel(
                "input.assay=='BE: General Eukayote'|input.assay=='BU: General Eukaryote'|input.assay=='MZ: Vascular Plants'|input.assay=='TP: Vascular Plants'|input.assay=='UM: Microbe'",
                pickerInput('subset_data_rest', "Select otu's",'All',multiple=F
                )),
              conditionalPanel(
#                "input.assay=='CI: Mostly insects'|input.assay=='RV: Vertebrates'|input.assay=='WV: Vertebrates'",
                  "input.assay=='RV: Vertebrates'|input.assay=='WV: Vertebrates'",
                pickerInput("subset_data_freshwater", "Select otu's", c("All","Just freshwater",'Fish (Order filter)','Fish (Genus filter)'),multiple = F))
              ,
              conditionalPanel(
                "input.assay=='CI: Mostly insects'",
                pickerInput("subset_data_insects", "Select otu's", c("All","Just freshwater",'EPT'),multiple = F))
              ,
               selectInput("level", "Select classification level (barchart colours only)",
                          choices =c("Class","Genus"), selected='CI') ,width = 3),
            
            #make it conditional on input type
            #need seperate plot labels to make conditional work (R doesnt like two thing with the same name)
            mainPanel(
              conditionalPanel(
                "input.analysis=='Site analysis'",
                column(6,plotOutput(outputId="plot1", width="1000px",height="600px"),uiOutput("lm1"))
                #    column(6,plotOutput(outputId="plot1"),uiOutput("lm1"))
              ),
              conditionalPanel(
                "input.analysis=='Site statistics'",
                column(6,plotOutput(outputId="plot2", width="1000px",height="600px"))
                #  column(6,plotOutput(outputId="plot2"))
              )
            )))

#that means I can probably modify how I display the plots ~ but maybe not as innext will have a dependency in it to skip if its not enough samples
