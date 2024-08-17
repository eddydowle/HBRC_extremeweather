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
library('gridExtra')

#input into phyloseq
setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
#records<-read.csv('HBRC_records_Eddy.csv')
#just cutting down to make it lighter for testing
#records<-read.csv('HBRC_records_Eddy_random50k.csv')
meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
#otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
#write.csv(otu_table,'HBRC_records_Eddy_random50k_otuapp.csv',row.names=F)
#write.csv(taxa_table,'HBRC_records_Eddy_random50k_taxaapp.csv',row.names=F)
#taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
otu_table<-read.csv('HBRC_records_Eddy_random50k_otuapp.csv')
taxa_table<-read.csv('HBRC_records_Eddy_random50k_taxaapp.csv')
freshwater_taxa<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')
#this is the table we will subset for choosing a primerset

# Define server logic required to draw a histogram
function(input, output, session) {

    #filtered otu table
    otu_table_subset<- reactive({
      otu_table_filtered<-otu_table %>% filter(PrimerSet == input$assay) %>% select(-PrimerSet) 
      if (input$subset_data_CI == "Just freshwater" | input$subset_data_RV == "Just freshwater" |input$subset_data_WV == "Just freshwater" ) {
        family<-freshwater_taxa %>% filter(Rank == 'Family')
        genus<-freshwater_taxa %>% filter(Rank == 'Genus')
        phylum<-freshwater_taxa %>% filter(Rank == 'Phylum')
        order<-freshwater_taxa %>% filter(Rank == 'Order')
        class<-freshwater_taxa %>% filter(Rank == 'Class')
        records_eukaryote_freshwater<-taxa_table %>% filter(Phylum %in% phylum$ID | Family %in% family$ID | Genus %in% genus$ID | Order %in% order$ID | Class %in% class$ID)
        otu_table_filtered_freshwater <- otu_table_filtered %>% filter(TaxID %in% records_eukaryote_freshwater$TaxID)
        return(otu_table_filtered_freshwater)
      }
      else{
        return(otu_table_filtered)
      }
      #maybe can do a if xx | xx | xx == 'Just Freshwater'
      #do filter
      #like this https://stackoverflow.com/questions/50926109/conditional-filtering-based-on-user-inputs-in-shiny ??
    })
    #testing table
  #  output[["table1"]] <- renderTable({otu_table_subset()%>%  slice_head(n=10)})
   # output$table1 <- renderTable({ otu_table %>% filter(PrimerSet == input$assay) %>% select(-PrimerSet) %>% group_by(UID, TaxID) %>% slice_head(n=10)  })
    #meta data for phyloseq
    meta_data_clean<- reactive({
      meta_data_2<-meta_data %>% mutate(location_date=paste(HBRC_Site_Name,CollectionDate,sep = '_')) %>% mutate(location_cyclone=paste(HBRC_Site_Name,Cyclone,sep = '_'))
      meta_data_cl <-meta_data_2 %>% filter(UID %in% otu_table_subset()$UID)
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
    
    count_samplesprecyclone<-reactive({
      nrow(meta_data_clean()[meta_data_clean()$Year == "2022"|meta_data_clean()$Year == "2021"|meta_data_clean()$Year == "2020",])
    })
    
    #plot across all 
   # output$plot1<-renderPlot({
      pt1<-reactive({
      if (input$analysis == "All locations" ) {
      p<-plot_bar(phyloseq_object(),  fill="Class")+
      geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
        ggtitle("Presence/Absence across sites") +
         theme(legend.key.size = unit(0.05, 'cm'))+
        geom_vline(aes(xintercept = count_samplesprecyclone()), linetype = "dashed", colour = "red",size = 1)+
      theme(legend.position="none")
      }
      if (input$analysis == "By location" ) {
        to_prune<-meta_data_clean() %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==input$site_choice) %>% pull(sample_id)
        phyloseq_object_prune<-prune_samples(to_prune,phyloseq_object())
#        phyloseq_object_prune<-prune_samples(phyloseq_object(), HBRC_Site_Name == input$site_choice)
        phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object_prune, 'pa') #pa dataset
        FisheDNA.ord <- ordinate(phyloseq_object_prune.pa, "NMDS", "jaccard") 
        p<-plot_ordination(phyloseq_object_prune, FisheDNA.ord, type="samples", color="CollectionDate")
        p + stat_ellipse(aes(x=NMDS1, y=NMDS2, colour= CollectionDate))+
          geom_point(size=2) +  coord_equal() +
          ggtitle("NMDS") +
          theme_classic()   
      }
      p
})
      pt2 <- reactive({
        if (input$analysis == "All locations") return(NULL)
        to_prune<-meta_data_clean() %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==input$site_choice)  %>% .[order(as.Date(.$CollectionDate, format="%d/%m/%Y")),] %>% pull(sample_id)
        phyloseq_object_prune<-prune_samples(to_prune,phyloseq_object())
        phyloseq_object_prune<-ps_reorder(phyloseq_object_prune, to_prune)
     #   https://github.com/joey711/phyloseq/issues/1075
   #reordering not working must be something to do with this but not sure
           pbar = plot_bar(phyloseq_object_prune, "CollectionDate", fill="Class")
        #do presence/absence plot
        #p1 = plot_bar(Fish.eDNA.pa, "location_date", fill="Class")
          #relevel to order chart
   pbar + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
          ggtitle("Presence/absence") +
          #  ggtitle("Read abundance across HBRC_Sites")+
          theme(legend.key.size = unit(0.03, 'cm'))+theme(axis.text=element_text(size=4),axis.title=element_text(size=10,face="bold"))
   
      })
      output$plot1 = renderPlot({
        ptlist <- list(pt1(),pt2())
        # remove the null plots from ptlist and wtlist
        to_delete <- !sapply(ptlist,is.null)
        ptlist <- ptlist[to_delete] 
        if (length(ptlist)==0) return(NULL)
        grid.arrange(grobs=ptlist,ncol=length(ptlist))
      })
    
}
