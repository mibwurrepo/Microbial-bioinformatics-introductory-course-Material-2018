---
title: "OPEN & REPRODUCIBLE MICROBIOME DATA ANALYSIS SPRING SCHOOL 2018"
author: "Sudarshan"
date: "2020-04-11"
output: bookdown::gitbook
site: bookdown::bookdown_site
---

# Beta diversity metrics  

Beta-diversity: Measures for differences between samples from different groups to identify if there are differences in the overall community composition and structure.   

**Load packages and data**  


```r
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  
```


For more information:

[Ordination](http://ordination.okstate.edu/overview.htm).  

[Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531).  

[Normalisation and data transformation](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y).  

[What is Constrained and Unconstrained Ordination](http://www.davidzeleny.net/anadat-r/doku.php/en:ordination).  

[Microbiome Datasets Are Compositional: And This Is Not Optional](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)  

[Compositional analysis: a valid approach to analyze microbiome high-throughput sequencing data](http://www.nrcresearchpress.com/doi/full/10.1139/cjm-2015-0821#.WvIJZIiFM2w)  



```r
# read non rarefied data
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")

# read rarefied data
ps0.rar.rds <- readRDS("./phyobjects/ps0.rar.rds")
# use print option to see the data saved as phyloseq object.
```

## Phylogenetic beta-diversity metrics  

### Unweighted Unifrac  

Unweighted Unifrac is based on presence/absence of different taxa and abundance is not important. However, it is sensitive to the sequencing depth. If a sample is sequenced more than the others then it may have many OTUs (most of them unique) consequently affecting the unifrac dissimilarity estimation.  

Usually, using subOTU/ASV approaches many singletons/OTUs with very low reads are discarded. If you have you won data you try the following code. For data from NG-tax we will skip this step as we have no singletons.    


```r
# if we remove OTUs that are detected atleast 10 times in 5% of the samples
ps0.rar.filtered <- core(ps0.rar.rds, detection = 10, prevalence = 0.05)

summarize_phyloseq(ps0.rar.filtered)

# we reduce the sparsity considerably. 
```

*For data from OTU picking*  
Since the data used here consists of different body sites with distinct biological properties, the results of ordination do not change a lot by filtering "rare" OTUs. Once again, knowing the biology of your samples and making choices rationally and documenting them is crucial.  
Feel free to use the OTU-picking strategy phyloseq object to investigate yourself.  


```r
ordu.unwt.uni <- ordinate(ps0.rar.rds, "PCoA", "unifrac", weighted=F)

# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

unwt.unifrac <- plot_ordination(ps0.rar.rds, 
                                     ordu.unwt.uni, color="scientific_name") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(unwt.unifrac)
```

<img src="05-MAW-PIV_files/figure-html/unnamed-chunk-4-1.png" width="672" />

Try repeating the above ordination using non-filtered phyloseq object.  

### Weighted Unifrac  

Weighted Unifrac will consider the abundances of different taxa.  


```r
ps1.rel <- microbiome::transform(ps1, "compositional")

ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)

# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

wt.unifrac <- plot_ordination(ps1.rel, 
                                     ordu.wt.uni, color="scientific_name") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(wt.unifrac)
```

<img src="05-MAW-PIV_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
print(wt.unifrac + stat_ellipse())
```

<img src="05-MAW-PIV_files/figure-html/unnamed-chunk-5-2.png" width="672" />

The figure brings forward an important characteristics of microbiome data called the 'Horse-shoe effect'. An investigation and explaination for this can be found in the article by Morton JT., et al. 2017 [Uncovering the Horseshoe Effect in Microbial Analyses](http://msystems.asm.org/content/2/1/e00166-16).   

You can repeat this analysis with phyloseq object from OTU-picking approach. 

Another important aspect regarding weighted unifrac is its property of having heavier weights for abundant taxa. To detect changes in moderately abundant lineages, an extenstion called generalized (UniFrac distance)(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413390/) has been developed. In this test data, we expect sufficient biological variation in composition between sites and hence, we do not apply GUniFrac.  

**To reiterate:** It is crucial to understand the biological features of the samples. Although these are exploratory approaches, it is important to differentiate between biological signal and technical artifacts.  

## Population-level Density landscapes    


```r
p <- plot_landscape(ps1.rel, 
                    "NMDS", 
                    "bray", 
                    col = "scientific_name") +
       labs(title = paste("NMDS / Bray-Curtis"))   

p <- p + scale_color_brewer(palette = "Dark2")+ scale_fill_gradient(low = "#e0ecf4", high = "#6e016b") 
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p 
```

```
## Warning: Position guide is perpendicular to the intended axis. Did you mean to
## specify a different guide `position`?
```

<img src="05-MAW-PIV_files/figure-html/unnamed-chunk-6-1.png" width="672" />

Bray-Curtis dissimilarity does not consider phylogenetic relationships between ASVs. There are several distance methods and a list can be obtained by typing `?distanceMethodList` in the console pane.   

Section on multivariate analysis will be discussed on Day3.  

## PERMANOVA  

Permutational multivariate analysis of variance [further reading](https://onlinelibrary.wiley.com/doi/10.1002/9781118445112.stat07841)  


```r
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.5-6
```

```
## 
## Attaching package: 'vegan'
```

```
## The following object is masked from 'package:microbiome':
## 
##     diversity
```

```r
metadf <- data.frame(sample_data(ps1.rel))

unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ scientific_name, data = metadf)

permanova
```

```
## 
## Call:
## adonis(formula = unifrac.dist ~ scientific_name, data = metadf) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
## scientific_name   3    71.912 23.9705  376.03 0.7059  0.001 ***
## Residuals       470    29.960  0.0637         0.2941           
## Total           473   101.872                 1.0000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Checking the homogeneity condition  

Type `?betadisper` in R console for more information.  


```r
ps.disper <- betadisper(unifrac.dist, metadf$scientific_name)
permutest(ps.disper, pairwise = TRUE)
```

```
## 
## Permutation test for homogeneity of multivariate dispersions
## Permutation: free
## Number of permutations: 999
## 
## Response: Distances
##            Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
## Groups      3 1.0014 0.33381 31.889    999  0.001 ***
## Residuals 470 4.9200 0.01047                         
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Pairwise comparisons:
## (Observed p-value below diagonal, permuted p-value above diagonal)
##                          human gut metagenome human oral metagenome
## human gut metagenome                                     1.0000e-03
## human oral metagenome              8.3027e-05                      
## human skin metagenome              2.8177e-06            2.4596e-02
## human vaginal metagenome           3.1013e-18            1.4078e-10
##                          human skin metagenome human vaginal metagenome
## human gut metagenome                1.0000e-03                    0.001
## human oral metagenome               2.5000e-02                    0.001
## human skin metagenome                                             0.006
## human vaginal metagenome            7.0275e-03
```




```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_Netherlands.1252  LC_CTYPE=English_Netherlands.1252   
## [3] LC_MONETARY=English_Netherlands.1252 LC_NUMERIC=C                        
## [5] LC_TIME=English_Netherlands.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] vegan_2.5-6        lattice_0.20-40    permute_0.9-5      dplyr_0.8.5       
##  [5] ggpubr_0.2.5       magrittr_1.5       RColorBrewer_1.1-2 microbiome_1.8.0  
##  [9] ggplot2_3.3.0      phyloseq_1.30.0   
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.0.0    xfun_0.12           reshape2_1.4.3     
##  [4] purrr_0.3.3         splines_3.6.3       rhdf5_2.30.1       
##  [7] vctrs_0.2.4         colorspace_1.4-1    htmltools_0.4.0    
## [10] stats4_3.6.3        mgcv_1.8-31         yaml_2.2.1         
## [13] survival_3.1-11     rlang_0.4.5         pillar_1.4.3       
## [16] withr_2.1.2         glue_1.3.2          BiocGenerics_0.32.0
## [19] foreach_1.4.8       lifecycle_0.2.0     plyr_1.8.6         
## [22] stringr_1.4.0       zlibbioc_1.32.0     Biostrings_2.54.0  
## [25] ggsignif_0.6.0      munsell_0.5.0       gtable_0.3.0       
## [28] codetools_0.2-16    evaluate_0.14       labeling_0.3       
## [31] Biobase_2.46.0      knitr_1.28          IRanges_2.20.2     
## [34] biomformat_1.14.0   parallel_3.6.3      Rcpp_1.0.4         
## [37] scales_1.1.0        S4Vectors_0.24.3    jsonlite_1.6.1     
## [40] XVector_0.26.0      farver_2.0.3        digest_0.6.25      
## [43] Rtsne_0.15          stringi_1.4.6       bookdown_0.18      
## [46] grid_3.6.3          ade4_1.7-15         tools_3.6.3        
## [49] tibble_2.1.3        cluster_2.1.0       tidyr_1.0.2        
## [52] crayon_1.3.4        ape_5.3             pkgconfig_2.0.3    
## [55] MASS_7.3-51.5       Matrix_1.2-18       data.table_1.12.8  
## [58] assertthat_0.2.1    rmarkdown_2.1       iterators_1.0.12   
## [61] Rhdf5lib_1.8.0      R6_2.4.1            multtest_2.42.0    
## [64] igraph_1.2.4.2      nlme_3.1-145        compiler_3.6.3
```





