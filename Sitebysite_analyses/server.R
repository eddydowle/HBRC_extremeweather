#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library('shiny')
library('microViz')
library('RColorBrewer')
library('tidyverse')
library('readxl')
library('devtools')
library('tidyverse')
library('viridis')
library('indicspecies') 
library('vegan')
library('glue')
library('reshape2')
library('phyloseq')
library('iNEXT.3D') 
library('ggplot2') 
library('microbiome')

#input into phyloseq
setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
records<-read.csv('HBRC_records_Eddy.csv')
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()

#this is the table we will subset for choosing a primerset

# Define server logic required to draw a histogram
function(input, output, session) {

    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')

    })
  
    #filtered otu table
    otu_table_subset<- reactive({
      otu_table %>% filter(PrimerSet == input$assay) %>% select(-PrimerSet) 
    })
    #testing table
  #  output[["table1"]] <- renderTable({otu_table_subset()%>%  slice_head(n=10)})
   # output$table1 <- renderTable({ otu_table %>% filter(PrimerSet == input$assay) %>% select(-PrimerSet) %>% group_by(UID, TaxID) %>% slice_head(n=10)  })
    #meta data for phyloseq
    meta_data_clean<- reactive({
      meta_data<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))
      meta_data_cl <-meta_data %>% filter(UID %in% otu_table_subset()$UID)
      row.names(meta_data_cl) <- meta_data_cl$UID
      meta_data_cl[1] <- NULL
      return(meta_data_cl)
    })
    #taxa table for phyloseq
    taxa_table_clean_phyloseq <- reactive ({
      taxa_table_clean <- taxa_table %>% filter(TaxID %in% otu_table_subset()$TaxID)
      taxa_table_clean<-taxa_table_clean[rowSums(is.na(taxa_table_clean)) != ncol(taxa_table_clean), ]
      row.names(taxa_table_clean) <- taxa_table_clean$TaxID
      taxa_table_clean[1] <- NULL
      taxa_table_clean$root<-'Root'
      taxa_table_clean <- taxa_table_clean %>%
        select(root, everything())
      return(as.matrix(taxa_table_clean))
    })
    
    #otu_table out_table_MZ_cut_wide
    otu_table_wide_clean_phyloseq <- reactive ({
      out_table_MZ_wide<-dcast(otu_table_subset(), TaxID~UID,value.var = "Count", fun.aggregate = sum)
      out_table_MZ_wide<-out_table_MZ_wide[complete.cases(out_table_MZ_wide), ]
      row.names(out_table_MZ_wide) <- out_table_MZ_wide$TaxID
      out_table_MZ_wide[1] <- NULL
      return(out_table_MZ_wide)
      })
    
#    output[["table1"]] <- renderTable({meta_data_clean() %>%  slice_head(n=10)})
  #make a phyloseq object and create figure to see if that is the issue
    
    phyloseq_object<-reactive({
      otu <- otu_table(otu_table_wide_clean_phyloseq(), taxa_are_rows = TRUE) 
      taxa <- tax_table(taxa_table_clean_phyloseq())
      sample <- sample_data(meta_data_clean())
      test<-meta_data_clean()[order(as.Date(meta_data_clean()$CollectionDate, format="%d/%m/%Y")),]
      test_names<-row.names(test)
      FisheDNA<-phyloseq(otu, taxa, sample)
      FisheDNA<-ps_reorder(FisheDNA, test_names)
      #remove singletons
      FisheDNA.5 <- filter_taxa(FisheDNA, function(x){sum(x > 0) > 1}, prune=TRUE) 
      return(FisheDNA.5)
    })
    #then should be able to build plots from phyloseq_object
    
    output[["table1"]] <- renderTable({taxa_names(phyloseq_object())[1:10]})
    #otu_table(phyloseq_object())[1:5, 1:5]
    #tax_table(phyloseq_object())[1:5, 1:4]
    #taxa_names(phyloseq_object())[1:10]
    
    #plot across all 
    output$plot1<-renderPlot({
      plot_bar(phyloseq_object(),  fill="Class")+
      geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
        ggtitle("Presence/Absence across sites") +
         theme(legend.key.size = unit(0.05, 'cm'))+
      theme(legend.position="none")
 
})
    
}
