#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#better site analysis no all samples

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
library('ranacapa')
library('stargazer')
#think use stargazer to render the beta diversity table
#https://groups.google.com/g/shiny-discuss/c/3VYuD6gjdk4?pli=1

#input into phyloseq
setwd("~/Dropbox (Personal)/Otago_2023 (1)/BiomeGene/extremeweatherevents/ASV_files/OneDrive_1_7-05-2024/")
samples<-read.csv('HBRC_samples_Eddy.csv')
records<-read.csv('HBRC_records_Eddy.csv')
#just cutting down to make it lighter for testing
#records<-read.csv('HBRC_records_Eddy_random50k.csv')
#meta_data<-read.table('Wilderlab_meta_out.txt',sep='\t',quote='',header =T)
meta_data<-read.table('Wilderlab_meta_out_clean.txt',sep='\t',quote='',header =T)
otu_table <- records %>% select(UID, PrimerSet, Count, TaxID)
#write.csv(otu_table,'HBRC_records_Eddy_random50k_otuapp.csv',row.names=F)
#write.csv(taxa_table,'HBRC_records_Eddy_random50k_taxaapp.csv',row.names=F)
taxa_table<-records %>% select(TaxID,Phylum,Class,Order,Family,Genus,Species) %>% distinct()
#otu_table<-read.csv('HBRC_records_Eddy_random50k_otuapp.csv')
#taxa_table<-read.csv('HBRC_records_Eddy_random50k_taxaapp.csv')
freshwater_taxa<-read.table('freshwater_invertsfish.txt',header=T,row.names=NULL,sep='\t')
#this is the table we will subset for choosing a primerset

# Define server logic required to draw a histogram
function(input, output, session) {
  
  #filtered otu table
  otu_table_subset<- reactive({
    otu_table_filtered<-otu_table %>% filter(PrimerSet == input$assay) %>% select(-PrimerSet) 
    if (input$subset_data_CI == "Just freshwater" | input$subset_data_RV == "Just freshwater" |input$subset_data_WV == "Just freshwater" ) {
      print('just freshwater')
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
    #print(meta_data_cl)
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
    #FisheDNA<-ps_reorder(FisheDNA, test_names)
    FisheDNA.5<-ps_reorder(FisheDNA, test_names)
    #remove singletons (moving this down to after other filters)
#    FisheDNA.5 <- filter_taxa(FisheDNA, function(x){sum(x > 0) > 1}, prune=TRUE)          
    to_prune<-meta_data_clean() %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==input$site_choice) %>% pull(sample_id)
    FisheDNA.5<-prune_samples(to_prune,FisheDNA.5)
    #filter those with less than 5 replicates
    #this is important when the taxa filters come on as it causes some samples to drop out if there is no target taxa ~ which means some samples will drop between taxa filters but this seems the best way around it. 
    filter<-as.data.frame(sample_data(FisheDNA.5)) %>% group_by(CollectionDate) %>% summarise(count=n()) %>% filter(count >=5) 
    remove_idx_test = as.character(get_variable(FisheDNA.5, "CollectionDate")) %in% filter$CollectionDate   
 #   print(remove_idx_test)
    FisheDNA.5 <- prune_samples(remove_idx_test, FisheDNA.5)
    #remove singletons
    FisheDNA.5 <- filter_taxa(FisheDNA.5, function(x){sum(x > 0) > 1}, prune=TRUE)          
    #for whatever reason subset_samples doesnt work in shinny so have to use the prune samples work around
        #FisheDNA.test<-subset_samples(FisheDNA.5, CollectionDate %in% filter$CollectionDate)
    #print(sample_names(FisheDNA.5))
   # print(otu_table(FisheDNA.5))
    return(FisheDNA.5)
  })
  #then should be able to build plots from phyloseq_object

#im probably going to have to put a minimum taxa count on that switches innext on or off as it seems to bug out when there is <5ish taxa?
  
#generate a innext object
#if we skip when there is <5 samples passing filter...then would need to fix the colours later on
  inext_object<-reactive({
    if (nrow(otu_table(phyloseq_object())) < 10) {
      return(NULL)
    }
    else {
    phyloseq.pa <- microbiome::transform(phyloseq_object(), 'pa')
   # print(otu_table(phyloseq_object()))
    print(nrow(otu_table(phyloseq_object())))
   categories <- unique(sample_data(phyloseq.pa)$CollectionDate)
  split_physeq_list <- list()
  print(categories)
  #subselected_categories<-NULL
for (category in categories) {
  #  print(category)
    remove_idx = as.character(get_variable(phyloseq.pa, "CollectionDate")) == category
    sub_physeq.100 <- prune_samples(remove_idx, phyloseq.pa)
 #   print(otu_table(sub_physeq.100))
    split_physeq_list[[category]] <- otu_table(sub_physeq.100)
}
  matrix_list <- lapply(split_physeq_list, function(x) {
    otu_table <- as(x, "matrix")
    return(otu_table)
  })
  
  matrix_list<-NULL
  matrix_list <- list(data = list())
  for (category in categories) {
    otu_table <- as(otu_table(split_physeq_list[[category]]), "matrix")
    matrix_list[["data"]][[category]] <- otu_table
  }
 # print(matrix_list$data)
  out.raw <- iNEXT3D(data = matrix_list$data, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_raw', nboot = 50)
  print('endinnext')
  return(out.raw)}
  #type 1 for site stats
 #ggiNEXT3D(out.raw, type = 1, facet.var = 'Assemblage') + facet_wrap(~Assemblage, nrow = 3)
#  ggiNEXT3D(out.raw, type = 1, facet.var = "Order.q")
 #type 2 for sample completeness curve
  # ggiNEXT3D(out.raw, type = 2, facet.var = "Order.q", color.var = "Assemblage")
  })
  
  #plot across all 
  # output$plot1<-renderPlot({
  pt1<-reactive({
    if (input$analysis == "Site statistics" ) {
      p<-ggrare(phyloseq_object(), color = "CollectionDate", label = "Sample",step = 100, se = FALSE) + theme_bw()
      p$data$CollectionDate <- factor(p$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(p$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
      p<-p+ theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_text(size = 8),legend.text = element_text(size = 7)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Rarefraction Curve")+ labs(colour = "Sampling Date")#guides(fill=guide_legend(title='Sampling Date'))

    }
    if (input$analysis == "Site analysis" ) {
      phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa') #pa dataset
      FisheDNA.ord <- ordinate(phyloseq_object_prune.pa, "NMDS", "jaccard") 
      p<-plot_ordination(phyloseq_object(), FisheDNA.ord, type="samples", color="CollectionDate") + theme_bw()+ stat_ellipse(type = "norm", linetype = 2) +ggtitle("NMDS")+ labs(colour = "Sampling Date")+
      theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 8),axis.title=element_text(size=11)) 
      p$data$CollectionDate <- factor(p$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(p$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
    }
    return(p)
  })
  pt2 <- reactive({
    if (input$analysis == "Site statistics") 
    {
      if (is.null(inext_object())){
        print('isnull')
        pbar<-NULL
      }
      else{
     pbar<- ggiNEXT3D(inext_object(), type = 2, facet.var = "Order.q", color.var = "Assemblage")
     #this does not solve all the issues for some reason I cant change the factor of the shape but can the colour so shifting the colour over to ordered and just turning of the shape legend. 
     pbar$data$col <- factor(pbar$data$col, levels = str_remove(format(sort(as.Date(unique(pbar$data$shape), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
     #have tried the same on Assemblage and shape columns in data table but no dice, not sure why its not working. Also tried changing factor of inext.3d object out.raw and didnt work either
    # pbar<-pbar+ theme(plot.title = element_text(size=8),text = element_text(size = 8))+theme(legend.position = "right") +theme(legend.text=element_text(size=6))+guides(shape = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +theme(legend.key.size = unit(0.5,"line"))
     
     pbar<-pbar+ theme(legend.text = element_text(size = 8),title =element_text(size=12)) + ggtitle("Sampling Curve")+theme(legend.position = "right")+ theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 8),axis.title=element_text(size=11)) +theme(legend.title=element_blank())+theme(legend.key.size = unit(0.1,"line"))+theme(legend.key.height=unit(0.1, "cm"))+guides(shape = "none")
     # pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
     #not sure how to reorder by date for innext plot object
     #trying to reset the colours to ggplots default colours
     gg_color_hue <- function(n) {
       hues = seq(15, 375, length = n + 1)
       hcl(h = hues, l = 65, c = 100)[1:n]
     }
     n = length(unique(pbar$data$shape))
     cols = gg_color_hue(n)
     pbar<-pbar+scale_color_manual(values = cols)+scale_fill_manual(values = cols)
      }
    }

    if (input$analysis == "Site analysis") 
    {    
      phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa') #pa dataset
      rich<-estimate_richness(phyloseq_object_prune.pa, measures=c(input$diversity_measure))
      Fish_Sample = as(sample_data(phyloseq_object_prune.pa), "matrix")
      rich<-cbind(rich, Fish_Sample)
      pbar<-ggplot(rich, aes(x=CollectionDate, y=!! sym(input$diversity_measure), colour = CollectionDate)) +
        geom_boxplot() + ggtitle("Species Richness") +
        ylab(paste(input$diversity_measure,"richness estimate"))+
        theme_bw() +
        theme(axis.text.x = element_text(size = 8,angle = 90, hjust=1))+ labs(colour = "Sampling Date",x="Sampling Date")+theme(        axis.text.y = element_text(size = 8),axis.title=element_text(size=11))
      #+theme(axis.text=element_text(size=4),axis.title=element_text(size=12))
     pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
    }
    return(pbar)
      })
  pt3 <- reactive({
    if (input$analysis == "Site statistics") {    
      to_prune<-meta_data_clean() %>% tibble::rownames_to_column(., "sample_id") %>% filter(HBRC_Site_Name==input$site_choice)  %>% .[order(as.Date(.$CollectionDate, format="%d/%m/%Y")),] %>% pull(sample_id)
    phyloseq_object_prune<-ps_reorder(phyloseq_object(), to_prune)
    ptest = plot_bar(phyloseq_object_prune, "CollectionDate", fill=input$level) + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") + ggtitle("Read Abundance") + theme_bw()+  theme(legend.position="none")
    test_data<-ptest$data
    #labeling only the top ten classes or genera
    top_10class<-test_data %>% group_by(!! sym(input$level)) %>% drop_na(!! sym(input$level)) %>%  summarise(Abundance = sum(Abundance,na.rm=TRUE)) %>% arrange(desc(Abundance)) %>% slice_max(Abundance, n = 10,with_ties=FALSE) %>% mutate(top10group= !! sym(input$level))
#    print(top_10class)
    #maybe change the column name then dont have the issue of selecting a column using input
    data_fish_top10 <- test_data %>% 
      mutate(top10 = ifelse (! (!! sym(input$level)) %in% top_10class$top10group, "Other", !! sym(input$level)))
 #   print(data_fish_top10)
    #make a colour pallet of 11 colours (grey at the end for other)
    col_brew<-c(brewer.pal(n = 10, name = "Paired"),"#808080")
    #relevel other to end
    data_fish_top10$top10 <- forcats::fct_relevel(data_fish_top10$top10, "Other", after = Inf)
    pbar<-ggplot(data_fish_top10, aes(fill=top10, y=Abundance, x=CollectionDate)) + 
      geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=col_brew)+ guides(fill=guide_legend(title=paste("Top 10",input$level)))+theme_bw()+ theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_text(size = 8),legend.text = element_text(size = 7)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Read Abundance")+ labs(x = "Sample")
    
   # +theme(legend.title = element_text(size = 8),legend.text = element_text(size = 6)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Read Abundance")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
    }
    
    if (input$analysis == "Site analysis") {
      if (is.null(inext_object())){
        print('isnull')
        pbar<-NULL
      }
      else{
    #  pbar<- ggiNEXT3D(inext_object(), type = 1, facet.var = 'Assemblage') + facet_wrap(~Assemblage, nrow = 3)
      pbar<- ggiNEXT3D(inext_object(), type = 1, facet.var = "Order.q")
      #this does not solve all the issues for some reason I cant change the factor of the shape but can the colour so shifting the colour over to ordered and just turning of the shape legend. 
      pbar$data$col <- factor(pbar$data$col, levels = str_remove(format(sort(as.Date(unique(pbar$data$shape), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
      #have tried the same on Assemblage and shape columns in data table but no dice, not sure why its not working. Also tried changing factor of inext.3d object out.raw and didnt work either
      pbar<-pbar + theme(legend.text = element_text(size = 8),title =element_text(size=12)) + ggtitle("Hill Numbers")+ theme(strip.text.x = element_text(size = 6))+theme(legend.position = "right")+ theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 8),axis.title=element_text(size=11)) +theme(legend.title=element_blank())+theme(legend.key.size = unit(0.1,"line"))+theme(legend.key.height=unit(0.1, "cm"))+guides(shape = "none")
      #legend.title = element_text(size = 8),
    #  pbar<-pbar+ theme(plot.title = element_text(size=8),text = element_text(size = 8))+theme(legend.position = "right") +theme(legend.text=element_text(size=6))+guides(shape = "none") + theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1)) +theme(legend.key.size = unit(0.5,"line")) 
     # + theme(legend.title = element_text(size = 8),legend.text = element_text(size = 6)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Relative Read Abundance")+ theme(strip.text.x = element_text(size = 6))
      
      #cant seem to change the legend size without stuffing up the rest of the legend
      # pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
      #not sure how to reorder by date for innext plot object
      #trying to reset the colours to ggplots default colours
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      n = length(unique(pbar$data$shape))
      cols = gg_color_hue(n)
      pbar<-pbar+scale_color_manual(values = cols)+scale_fill_manual(values = cols)
      }
   #   return(pbar)
    }
    return(pbar)
     })
  
  pt4 <- reactive({
    if (input$analysis == "Site statistics") { return(NULL)}
    if (input$analysis == "Site analysis") {
      #bar chart of species proportions per site
      #abundances from P/A dataset
 #     phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa') #pa dataset
 #     FisheDNA.6 = transform_sample_counts(phyloseq_object_prune.pa, function(x) x / sum(x) )
     #raw abundances not P/A
      FisheDNA.6 = transform_sample_counts(phyloseq_object(), function(x) x / sum(x) )
      data_fish <- psmelt(FisheDNA.6)
      #find top 10 classes in dataset
      top_10class<-data_fish %>% group_by(!! sym(input$level)) %>% drop_na(!! sym(input$level)) %>%  summarise(Abundance = sum(Abundance,na.rm=TRUE)) %>% arrange(desc(Abundance)) %>% slice_max(Abundance, n = 10,with_ties=FALSE) %>% mutate(top10group= !! sym(input$level))
      #    print(top_10class)
      #maybe change the column name then dont have the issue of selecting a column using input
      data_fish_top10 <- data_fish %>% 
        mutate(top10 = ifelse (! (!! sym(input$level)) %in% top_10class$top10group, "Other", !! sym(input$level)))
      
 #     top_10class<-data_fish %>% group_by(Class) %>% drop_na(Class) %>%  summarise(Abundance = sum(Abundance,na.rm=TRUE)) %>% arrange(desc(Abundance)) %>% slice_max(Abundance, n = 10,with_ties=FALSE)
      #then everything that is not in top_10class gets mutated into 'other' category
      #https://stackoverflow.com/questions/71595962/conditional-mutate-by-matching-strings-or-characters
    #  data_fish_top10 <- data_fish %>% 
     #   mutate(Class_top10 = ifelse (!Class %in% top_10class$Class, "Other", Class))
      #make a colour pallet of 11 colours (grey at the end for other)
      col_brew<-c(brewer.pal(n = 10, name = "Paired"),"#808080")
      #relevel other to end
      data_fish_top10$top10 <- forcats::fct_relevel(data_fish_top10$top10, "Other", after = Inf)
     pbar <- ggplot(data=data_fish_top10, aes(x=Sample, y=Abundance, fill=top10)) + facet_grid(~CollectionDate, scales = "free")
      pbar<-pbar + geom_bar(aes(), stat="identity", position="fill") +theme_bw()+ theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=col_brew)+ guides(fill=guide_legend(title=paste("Top 10",input$level)))+ theme(legend.title = element_text(size = 8),legend.text = element_text(size = 7)) +theme(legend.key.height=unit(0.3, "cm"))+ ggtitle("Relative Read Abundance")+ theme(strip.text.x = element_text(size = 6,angle = 90))
      
    #  pbar <- ggplot(data=data_fish, aes(x=Sample, y=Abundance, fill=Class)) + facet_grid(~CollectionDate, scales = "free") + geom_bar(aes(), stat="identity", position="fill") + ggtitle("Relative abundance P/A") + theme(legend.key.size = unit(0.03, 'cm')) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=6))
        #theme(axis.text=element_text(size=4),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=10,face="bold"))
 pbar$data$CollectionDate <- factor(pbar$data$CollectionDate, levels = str_remove(format(sort(as.Date(unique(pbar$data$CollectionDate), format="%d/%m/%Y")),"%d/%m/%y"),"^0+") )
    }
    return(pbar)
    })
  
mod <- reactive({
 #   print('test')
    if (input$analysis == "Site statistics") { 
      return(NULL)}
    if (input$analysis == "Site analysis") {
      phyloseq_object_prune.pa <- microbiome::transform(phyloseq_object(), 'pa')
      phyloseq_object_prune.pa.OTU <- as(otu_table(phyloseq_object_prune.pa), "matrix")
      phyloseq_object_prune.pa_OTU_df = as.data.frame(phyloseq_object_prune.pa.OTU)
      phyloseq_object_prune.pa_OTU_df= t(phyloseq_object_prune.pa_OTU_df)
      phyloseq_object_prune.pa_Sample = as(sample_data(phyloseq_object_prune.pa), "matrix")
      phyloseq_object_prune.pa_Sample_df =  as.data.frame(phyloseq_object_prune.pa_Sample)
      OTU_Dist <- vegdist(phyloseq_object_prune.pa_OTU_df, method="jaccard")
      mod_beta<-adonis2(OTU_Dist ~ CollectionDate, method = "jaccard", data = phyloseq_object_prune.pa_Sample_df, permutations = 9999)
#      print(mod_beta)
 #     print(OTU_Dist)
    }
    return(mod_beta)
  })

# print(pt5)
#?stargazer
output$lm1 <- renderText(HTML(stargazer(mod(), type="html",title="PERMOVA OTU distance matrix ~ Samples")))  
#output$modelSummary <- renderPrint({
#  mod()
#})

#cause its now conditional on the table I need to have seperate plot labels for the UI side
output$plot1 = renderPlot({
    ptlist <- list(pt1(),pt2(),pt3(),pt4())
    # remove the null plots from ptlist and wtlist
    to_delete <- !sapply(ptlist,is.null)
    ptlist <- ptlist[to_delete] 
    if (length(ptlist)==0) return(NULL)
    grid.arrange(grobs=ptlist,ncol=length(ptlist))
  })
  output$plot2 = renderPlot({
    ptlist <- list(pt1(),pt2(),pt3(),pt4())
    # remove the null plots from ptlist and wtlist
    to_delete <- !sapply(ptlist,is.null)
    ptlist <- ptlist[to_delete] 
    if (length(ptlist)==0) return(NULL)
    grid.arrange(grobs=ptlist,ncol=length(ptlist))
  })
}

