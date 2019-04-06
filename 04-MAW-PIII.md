---
title: "OPEN & REPRODUCIBLE MICROBIOME DATA ANALYSIS SPRING SCHOOL 2018"
author: "Sudarshan"
date: "2019-04-06"
output: bookdown::gitbook
site: bookdown::bookdown_site
---

# Composition plots  

Barplots are a one way of visualising the composition of your samples.    

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

# if you have dada2/deblur output and sequences as taxa names, then you can change them as follows
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

# colourCount = length(unique(taxic$Family))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.
```

```
## [1] "Domain" "Phylum" "Class"  "Order"  "Family" "Genus"  "OTU"
```

```r
taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object


# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL

# Second merge at family level

ps1.com.fam <- microbiome::aggregate_top_taxa(ps1.com, "Family", top = 10)

plot.composition.COuntAbun <- plot_composition(ps1.com.fam) + theme(legend.position = "bottom") +
  scale_fill_brewer("Family", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

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
                                             x.label = "env_material") 
plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Family", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)
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
## [1] "Tax"       "Sample"    "Abundance" "xlabel"
```


```r
p.com <- ggplot(data.com, aes(x = Sample, y = Abundance, fill = Tax))
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
p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance)) 

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
ps_df <- microbiomeutilities::phy_to_ldf(ps1.com.fam, 
                                         transform.counts = "compositional")
```

```
## Warning: replacing previous import 'ggplot2::alpha' by 'microbiome::alpha'
## when loading 'microbiomeutilities'
```

```
## An additonal column Sam_rep with sample names is created for reference purpose
```

```r
colnames(ps_df)
```

```
##  [1] "OTUID"                       "Family"                     
##  [3] "unique"                      "Sam_rep"                    
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
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_Netherlands.1252   
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] bindrcpp_0.2.2     dplyr_0.7.7        ggpubr_0.1.8      
## [4] magrittr_1.5       RColorBrewer_1.1-2 microbiome_1.5.28 
## [7] ggplot2_3.1.0      phyloseq_1.24.2   
## 
## loaded via a namespace (and not attached):
##  [1] Biobase_2.40.0             viridis_0.5.1             
##  [3] tidyr_0.8.2                jsonlite_1.5              
##  [5] viridisLite_0.3.0          splines_3.5.1             
##  [7] foreach_1.4.4              assertthat_0.2.0          
##  [9] stats4_3.5.1               yaml_2.2.0                
## [11] ggrepel_0.8.0              pillar_1.3.0              
## [13] backports_1.1.2            lattice_0.20-35           
## [15] glue_1.3.0                 digest_0.6.18             
## [17] XVector_0.20.0             colorspace_1.3-2          
## [19] htmltools_0.3.6            Matrix_1.2-15             
## [21] plyr_1.8.4                 microbiomeutilities_0.99.0
## [23] pkgconfig_2.0.2            pheatmap_1.0.10           
## [25] bookdown_0.7               zlibbioc_1.26.0           
## [27] purrr_0.2.5                scales_1.0.0              
## [29] Rtsne_0.15                 tibble_1.4.2              
## [31] mgcv_1.8-25                IRanges_2.14.12           
## [33] withr_2.1.2                BiocGenerics_0.26.0       
## [35] lazyeval_0.2.1             survival_2.43-1           
## [37] crayon_1.3.4               evaluate_0.12             
## [39] nlme_3.1-137               MASS_7.3-51.1             
## [41] vegan_2.5-3                tools_3.5.1               
## [43] data.table_1.11.8          stringr_1.3.1             
## [45] Rhdf5lib_1.2.1             S4Vectors_0.18.3          
## [47] munsell_0.5.0              ggsci_2.9                 
## [49] cluster_2.0.7-1            Biostrings_2.48.0         
## [51] ade4_1.7-13                compiler_3.5.1            
## [53] rlang_0.3.0.1              rhdf5_2.24.0              
## [55] grid_3.5.1                 iterators_1.0.10          
## [57] biomformat_1.8.0           igraph_1.2.2              
## [59] labeling_0.3               rmarkdown_1.10            
## [61] gtable_0.2.0               codetools_0.2-15          
## [63] multtest_2.36.0            reshape2_1.4.3            
## [65] R6_2.3.0                   gridExtra_2.3             
## [67] knitr_1.20                 bindr_0.1.1               
## [69] rprojroot_1.3-2            permute_0.9-4             
## [71] ape_5.2                    stringi_1.2.4             
## [73] parallel_3.5.1             Rcpp_0.12.19              
## [75] tidyselect_0.2.5           xfun_0.4
```


