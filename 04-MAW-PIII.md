---
title: "OPEN & REPRODUCIBLE MICROBIOME DATA ANALYSIS SPRING SCHOOL 2018"
author: "Sudarshan"
date: "2018-05-28"
output: bookdown::gitbook
site: bookdown::bookdown_site
---

# Composition plots  

Barplots are a simple way of visualising the composition of your samples.    

We will use the filtered phyloseq object from **Set-up and Pre-processing** section.  

**Load packages** 


```r
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  
```

  

```r
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")

# use print option to see the data saved as phyloseq object.

print(ps1)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4710 taxa and 474 samples ]
## sample_data() Sample Data:       [ 474 samples by 30 sample variables ]
## tax_table()   Taxonomy Table:    [ 4710 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 4710 tips and 4709 internal nodes ]
```

## Barplot counts 


```r
ps1.com <- ps1
taxa_names(ps1.com) <- paste0("Seq_", rownames(tax_table(ps1.com)))

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table)  # this will help in setting large color options

#colourCount = length(unique(taxic$Family))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic)  # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic)  # You can see that we now have extra taxonomy levels.
```

```
## [1] "Domain" "Phylum" "Class"  "Order"  "Family" "Genus"  "OTU"
```

```r
taxmat <- as.matrix(taxic)  # convert it into a matrix.
new.tax <- tax_table(taxmat)  # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax  # incroporate into phyloseq Object



# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"


# We will also remove the 'f__' patterns for cleaner labels
# tax_table(ps1.com)[, colnames(tax_table(ps1.com))] <- gsub(tax_table(ps1.com)[, 
#    colnames(tax_table(ps1.com))], pattern = "[a-z]__", replacement = "")

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
    face = "italic", colour = "Black", angle = 0)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL

# Second merge at family level 

ps1.com.fam <- microbiome::aggregate_taxa(ps1.com, "Family", top = 10)

plot.composition.COuntAbun <- plot_composition(ps1.com.fam) + theme(legend.position = "bottom") + 
    scale_fill_brewer("Family", palette = "Paired") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size=18))
  
plot.composition.COuntAbun
```

<img src="04-MAW-PIII_files/figure-html/unnamed-chunk-3-1.png" width="1920" />

```r
#ggsave("./Test_Outputfiles/Family_barplot_CountAbundance.pdf", height = 6, width = 8)
```

This plot is based on the reads per sample. In the next step, we plot the relative abundance.

## Barplot relative abundance 

Make it relative abundance


```r
# the previous pseq object ps1.com.fam is only counts.

# Use traqnsform function of microbiome to convert it to rel abun.
ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")

plot.composition.relAbun <- plot_composition(ps1.com.fam.rel, 
                                             sample.sort = "scientific_name", 
                                             x.label = "env_material") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size=18))
  
plot.composition.relAbun
```

<img src="04-MAW-PIII_files/figure-html/unnamed-chunk-4-1.png" width="1920" />

```r
#ggsave("./figures/Family_barplot_RelAbundance.pdf", height = 6, width = 8)
```

### Barplot customize 


```r
data.com <- plot.composition.relAbun$data
colnames(data.com)
```

```
## [1] "OTU"       "Sample"    "Abundance" "xlabel"
```


```r
p.com <- ggplot(data.com, aes(x = Sample, y = Abundance, fill = OTU))
p.com <- p.com + geom_bar(position = "stack", stat = "identity")
p.com <- p.com + scale_x_discrete(labels = data.com$xlabel, breaks = data.com$Sample)
p.com <- p.com + facet_grid(~xlabel, scales = "free") + theme_bw()
p.com <- p.com + scale_fill_brewer("Family", palette = "Paired") 
p.com <- p.com + rremove("x.text") 

ggsave("./figures/Composition plots.pdf", height = 4, width = 6)
```

For more information [Microbiome tutorial](http://microbiome.github.io/microbiome/Composition.html)   
## Heatmaps  
These are a good alternative to barplots (if done right).  


```r
# base plot
p.heat <- ggplot(data.com, aes(x = Sample, y = OTU)) + geom_tile(aes(fill = Abundance)) 

# Change color
p.heat <- p.heat + scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw() 

# Make bacterial names italics
p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black', 
                                                    size = 10, 
                                                    face = 'italic')) 
# Make seperate samples based on main varaible
p.heat <- p.heat + facet_grid(~xlabel, 
                              scales = "free") + rremove("x.text") 

p.heat <- p.heat + ylab("Family")

#Clean the x-axis
p.heat <- p.heat + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) 

# Clean the facet label box
p.heat <- p.heat + theme(legend.key = element_blank(), 
                     strip.background = element_rect(colour="black", fill="white"))

print(p.heat)
```

<img src="04-MAW-PIII_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
ggsave("./figures/Heatmap.pdf", height = 4, width = 6)


# + geom_text(aes(label = round(Abundance)), size = 0.4)
```


**Extra**  

Following is an example of customizing the plot using ggpubr.  


```r
# we use count data at family level from the barplot for counts
ps_df <- microbiomeutilities::phy_to_ldf(ps1.com.fam, transform.counts = "compositional")
```

```
## An additonal column Sam_rep with sample names is created for reference purpose
```

```r
colnames(ps_df)
```

```
##  [1] "OTUID"                       "Family"                     
##  [3] "OTU"                         "Sam_rep"                    
##  [5] "Abundance"                   "BarcodeSequence"            
##  [7] "LinkerPrimerSequence"        "run_prefix"                 
##  [9] "body_habitat"                "body_product"               
## [11] "body_site"                   "bodysite"                   
## [13] "dna_extracted"               "elevation"                  
## [15] "env"                         "env_biome"                  
## [17] "env_feature"                 "env_material"               
## [19] "env_package"                 "geo_loc_name"               
## [21] "host_common_name"            "host_scientific_name"       
## [23] "host_subject_id"             "host_taxid"                 
## [25] "latitude"                    "longitude"                  
## [27] "physical_specimen_location"  "physical_specimen_remaining"
## [29] "psn"                         "public"                     
## [31] "sample_type"                 "scientific_name"            
## [33] "sequencecenter"              "title"                      
## [35] "Description"
```

```r
# this data.frame can be used to customize several plots.  

# example boxplot at family level

p.box <- ggstripchart(ps_df, "scientific_name", "Abundance", 
                      facet.by = "Family", color = "body_product",
                      palette = "jco"
                      )

p.box + rremove("x.text")
```

<img src="04-MAW-PIII_files/figure-html/unnamed-chunk-8-1.png" width="672" />



```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 16299)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] bindrcpp_0.2.2       dplyr_0.7.5          ggpubr_0.1.6        
## [4] magrittr_1.5         RColorBrewer_1.1-2   microbiome_1.1.10013
## [7] ggplot2_2.2.1.9000   phyloseq_1.22.3     
## 
## loaded via a namespace (and not attached):
##  [1] ggrepel_0.8.0              Rcpp_0.12.17              
##  [3] ape_5.1                    lattice_0.20-35           
##  [5] tidyr_0.8.1                Biostrings_2.46.0         
##  [7] assertthat_0.2.0           rprojroot_1.3-2           
##  [9] digest_0.6.15              foreach_1.4.4             
## [11] R6_2.2.2                   plyr_1.8.4                
## [13] backports_1.1.2            stats4_3.4.4              
## [15] evaluate_0.10.1            pillar_1.2.2              
## [17] zlibbioc_1.24.0            rlang_0.2.0               
## [19] lazyeval_0.2.1             data.table_1.11.2         
## [21] vegan_2.5-2                S4Vectors_0.16.0          
## [23] Matrix_1.2-12              rmarkdown_1.9             
## [25] labeling_0.3               splines_3.4.4             
## [27] Rtsne_0.13                 stringr_1.3.1             
## [29] pheatmap_1.0.10            igraph_1.2.1              
## [31] munsell_0.4.3              compiler_3.4.4            
## [33] xfun_0.1                   pkgconfig_2.0.1           
## [35] BiocGenerics_0.24.0        multtest_2.34.0           
## [37] mgcv_1.8-23                htmltools_0.3.6           
## [39] biomformat_1.6.0           tidyselect_0.2.4          
## [41] gridExtra_2.3              tibble_1.4.2              
## [43] bookdown_0.7               IRanges_2.12.0            
## [45] codetools_0.2-15           viridisLite_0.3.0         
## [47] permute_0.9-4              withr_2.1.2               
## [49] MASS_7.3-49                grid_3.4.4                
## [51] nlme_3.1-131.1             jsonlite_1.5              
## [53] gtable_0.2.0               scales_0.5.0              
## [55] stringi_1.2.2              XVector_0.18.0            
## [57] reshape2_1.4.3             viridis_0.5.1             
## [59] ggsci_2.9                  iterators_1.0.9           
## [61] tools_3.4.4                microbiomeutilities_0.99.0
## [63] ade4_1.7-11                Biobase_2.38.0            
## [65] glue_1.2.0                 purrr_0.2.4               
## [67] parallel_3.4.4             survival_2.41-3           
## [69] yaml_2.1.19                colorspace_1.3-2          
## [71] rhdf5_2.22.0               cluster_2.0.6             
## [73] knitr_1.20                 bindr_0.1.1
```


