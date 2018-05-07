# This code formats the data from the plot_core function. 
# As otu_tables have OTU ids as rownames, we need to add taxonomic infomration.
# This function does exactly that. 
# plotdata is that dataframe in the output of plot_core. 
# tax.table is the taxonomy table in phyloseq object. tax_table(phyloseqobject)
# levels is from which taxonomic level you want information. options if "ALL" then all levels are given, 
# else only till specified level for instance "Genus".
# p.core <- plot_core(ps1.stool.rel2, plot.type = "heatmap", 
#                     colours = rev(brewer.pal(5, "Spectral")),
#                     prevalences = prevalences, 
#                     detections = detections, 
#                     min.prevalence = .5) + 
#                     xlab("Detection Threshold (Relative Abundance (%))")
#
# plot_data <- p.core$data
#
# df_plot <- plot_core_assist(plot_data, tax_table(ps1.stool.rel2), levels = "Family")
# p.core$data <- df_plot

# plot(p.core + theme(axis.text.y = element_text(face="italic")))


plot_core_assist <- function(plotdata, tax.table, levels){
  
  df <- plotdata
  # get the list of OTUs
  list <- df$Taxa 
  
  tax <- as.data.frame(tax.table)
  # get the taxonomy data
  
  # add the OTus to last column
  tax$OTU <- rownames(tax)
  
  # select taxonomy of only 
  # those OTUs that are used in the plot
  tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
  
  if (levels == "ALL"){
    # We will merege all the column into one except the Doamin as all is bacteria in this case
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else if (levels == "Phylum") {
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else if (levels == "Class") {
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Class", "Order", "Family", "Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else if (levels == "Order") {
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Order", "Family", "Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else if (levels == "Family") {
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Family", "Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else if (levels == "Genus") {
    tax.unit <- tidyr::unite(tax2, Taxa_level,c("Genus", "Species", "OTU"), sep = "_;", remove = TRUE)
    
  } else {
    stop("Please specifiy a level: ALL, Phylum, Class, Order, Family or Genus")
  }

  tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
  
  # add this new information into the plot data df
  
  df$Taxa <- tax.unit$Taxa_level

  # replace the data in the plot object
  return(df)
 
}



