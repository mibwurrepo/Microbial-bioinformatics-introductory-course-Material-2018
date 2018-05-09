---
title: "OPEN & REPRODUCIBLE MICROBIOME DATA ANALYSIS SPRING SCHOOL 2018"
author: "Sudarshan"
date: "2018-05-09"
output: bookdown::gitbook
site: bookdown::bookdown_site
---

# Core microbiota  

For more information:

[The adult intestinal core microbiota is determined by analysis depth and health status](https://www.sciencedirect.com/science/article/pii/S1198743X14609629?via%3Dihub).  

[Intestinal microbiome landscaping: insight in community assemblage and implications for microbial modulation strategies](https://academic.oup.com/femsre/article/41/2/182/2979411).  

[Intestinal Microbiota in Healthy Adults: Temporal Analysis Reveals Individual and Common Core and Relation to Intestinal Symptoms](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0023035).  


```r
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  
```

## Core microbiota anlaysis   

We will use the filtered phyloseq object from previous tutorial. We will use the filtered phyloseq object from the first section for pre-processioning.


```r
# read non rarefied data
ps1 <- readRDS("./phyobjects/ps1.rds")

# use print option to see the data saved as phyloseq object.
```

Subset the data to keep only stool samples.  


```r
ps1.stool <- subset_samples(ps1, bodysite == "Stool")

# convert to relative abundance  
ps1.stool.rel <- microbiome::transform(ps1.stool, "compositional")
print(ps1.stool.rel)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 3690 taxa and 169 samples ]
## sample_data() Sample Data:       [ 169 samples by 31 sample variables ]
## tax_table()   Taxonomy Table:    [ 3690 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 3690 tips and 3689 internal nodes ]
```

```r
ps1.stool.rel2 <- prune_taxa(taxa_sums(ps1.stool.rel) > 0, ps1.stool.rel)

print(ps1.stool.rel2)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2717 taxa and 169 samples ]
## sample_data() Sample Data:       [ 169 samples by 31 sample variables ]
## tax_table()   Taxonomy Table:    [ 2717 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 2717 tips and 2716 internal nodes ]
```

Check of the core OTUs  


```r
core.taxa.standard <- core_members(ps1.stool.rel2, detection = 0.001, prevalence = 50/100)

print(core.taxa.standard)
```

```
##  [1] "191483"  "194909"  "197214"  "177343"  "184753"  "560336"  "192684" 
##  [8] "535375"  "193233"  "772282"  "4325533" "850586"  "196219"  "180082" 
## [15] "330458"  "173996"
```

There are 16 OTUs that are core based on the cut-offs for prevalence and detection we choose. However, we only see IDs, not very informative. We can get the classification of these as below.    


```r
# Extract the taxonomy table

taxonomy <- as.data.frame(tax_table(ps1.stool.rel2))

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)

DT::datatable(core_taxa_id)
```

<!--html_preserve--><div id="htmlwidget-742f8a7241d4fec4e9c1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-742f8a7241d4fec4e9c1">{"x":{"filter":"none","data":[["191483","194909","197214","177343","184753","560336","192684","535375","193233","772282","4325533","850586","196219","180082","330458","173996"],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Bacteroidetes","Firmicutes","Firmicutes"],["Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Bacteroidia","Clostridia","Clostridia"],["Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Bacteroidales","Clostridiales","Clostridiales"],["Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Bacteroidaceae","Rikenellaceae","Rikenellaceae","[Odoribacteraceae]","Porphyromonadaceae","Porphyromonadaceae","Lachnospiraceae","Lachnospiraceae"],["Bacteroides","Bacteroides","Bacteroides","Bacteroides","Bacteroides","Bacteroides","Bacteroides","Bacteroides","Bacteroides",null,null,"Odoribacter","Parabacteroides","Parabacteroides",null,null],[null,null,null,null,null,null,null,"ovatus",null,null,null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>Family<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


## Core abundance and diversity  
Total core abundance in each sample (sum of abundances of the core members):


```r
core.abundance <- sample_sums(core(ps1.stool.rel2, detection = 0.001, prevalence = 50/100))

DT::datatable(as.data.frame(core.abundance))
```

<!--html_preserve--><div id="htmlwidget-4afee484723c7cf3ad7c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-4afee484723c7cf3ad7c">{"x":{"filter":"none","data":[["1927.SRS014659.SRX020511.SRR045649","1927.SRS015854.SRX020554.SRR043727","1927.SRS015599.SRX020539.SRR044148","1927.SRS015452.SRX020545.SRR044740","1927.SRS014442.SRX020511.SRR045656","1927.SRS014287.SRX020511.SRR045664","1927.SRS049428.SRX020670.SRR047663","1927.SRS019910.SRX020516.SRR044347","1927.SRS050422.SRX020511.SRR045628","1927.SRS064276.SRX020670.SRR046347","1927.SRS052697.SRX020550.SRR046441","1927.SRS063827.SRX020532.SRR046984","1927.SRS062847.SRX020532.SRR047036","1927.SRS014659.SRX020543.SRR045556","1927.SRS055482.SRX020550.SRR047716","1927.SRS021910.SRX020675.SRR047783","1927.SRS011271.SRX020661.SRR044827","1927.SRS011529.SRX020659.SRR044994","1927.SRS011271.SRX020659.SRR045051","1927.SRS018733.SRX020516.SRR047425","1927.SRS015281.SRX020554.SRR043805","1927.SRS011415.SRX020661.SRR044786","1927.SRS014287.SRX020543.SRR045558","1927.SRS046382.SRX020540.SRR047262","1927.SRS064321.SRX020519.SRR046501","1927.SRS056259.SRX020578.SRR047235","1927.SRS020470.SRX022097.SRR057663","1927.SRS015247.SRX020577.SRR044079","1927.SRS015724.SRX020539.SRR044132","1927.SRS015782.SRX020539.SRR044139","1927.SRS018872.SRX020572.SRR051576","1927.SRS016203.SRX020516.SRR044343","1927.SRS019089.SRX020564.SRR044731","1927.SRS011621.SRX020661.SRR044776","1927.SRS011413.SRX020661.SRR044792","1927.SRS011653.SRX020659.SRR044997","1927.SRS011413.SRX020659.SRR045020","1927.SRS019582.SRX020543.SRR045562","1927.SRS050422.SRX020543.SRR045575","1927.SRS048819.SRX020511.SRR045612","1927.SRS054461.SRX020550.SRR046442","1927.SRS042703.SRX020540.SRR047267","1927.SRS046639.SRX020670.SRR047659","1927.SRS023047.SRX020675.SRR047826","1927.SRS016056.SRX020536.SRR049357","1927.SRS065725.SRX020519.SRR049530","1927.SRS011410.SRX020661.SRR044800","1927.SRS014442.SRX020543.SRR045553","1927.SRS048819.SRX020543.SRR045568","1927.SRS045877.SRX020523.SRR046406","1927.SRS057447.SRX020523.SRR046411","1927.SRS046313.SRX020550.SRR046437","1927.SRS014613.SRX020546.SRR047386","1927.SRS021601.SRX020675.SRR047748","1927.SRS022987.SRX020675.SRR047766","1927.SRS022924.SRX020675.SRR047811","1927.SRS020176.SRX020676.SRR047859","1927.SRS024511.SRX020676.SRR047887","1927.SRS063807.SRX022232.SRR058094","1927.SRS063968.SRX022156.SRR058115","1927.SRS014823.SRX020546.SRR043699","1927.SRS018712.SRX020572.SRR046003","1927.SRS014885.SRX020554.SRR043822","1927.SRS014369.SRX020554.SRR043823","1927.SRS065466.SRX020532.SRR047056","1927.SRS050733.SRX020550.SRR047710","1927.SRS048870.SRX020577.SRR044064","1927.SRS015578.SRX020577.SRR044021","1927.SRS016018.SRX020539.SRR044156","1927.SRS015911.SRX020520.SRR044615","1927.SRS019089.SRX020545.SRR044757","1927.SRS011452.SRX020661.SRR044798","1927.SRS011586.SRX020661.SRR044826","1927.SRS053335.SRX020531.SRR046056","1927.SRS043411.SRX020670.SRR046345","1927.SRS045526.SRX020670.SRR046355","1927.SRS043299.SRX020523.SRR046410","1927.SRS049823.SRX020523.SRR046414","1927.SRS045607.SRX020523.SRR046421","1927.SRS063275.SRX020532.SRR046899","1927.SRS049164.SRX020540.SRR047273","1927.SRS042628.SRX020578.SRR047303","1927.SRS045414.SRX020563.SRR047697","1927.SRS015389.SRX020546.SRR043661","1927.SRS015190.SRX020554.SRR043731","1927.SRS019787.SRX020577.SRR043974","1927.SRS015332.SRX020577.SRR044071","1927.SRS015911.SRX020510.SRR044585","1927.SRS019013.SRX020510.SRR044597","1927.SRS019013.SRX020520.SRR044627","1927.SRS011410.SRX020659.SRR045027","1927.SRS048083.SRX020677.SRR047956","1927.SRS023488.SRX020678.SRR048028","1927.SRS050599.SRX020678.SRR048036","1927.SRS065665.SRX020519.SRR049404","1927.SRS016267.SRX020577.SRR044054","1927.SRS011452.SRX020659.SRR047485","1927.SRS012191.SRX020675.SRR047747","1927.SRS016095.SRX020516.SRR044317","1927.SRS014999.SRX020577.SRR044050","1927.SRS011653.SRX020661.SRR044770","1927.SRS019582.SRX020511.SRR045623","1927.SRS018559.SRX020516.SRR044352","1927.SRS019267.SRX020516.SRR044360","1927.SRS015960.SRX020577.SRR044087","1927.SRS063324.SRX020519.SRR046503","1927.SRS054590.SRX020550.SRR046427","1927.SRS011157.SRX020659.SRR045036","1927.SRS015663.SRX020539.SRR044161","1927.SRS020119.SRX020676.SRR047844","1927.SRS015815.SRX020562.SRR044705","1927.SRS046349.SRX020531.SRR046051","1927.SRS015452.SRX020564.SRR047456","1927.SRS011529.SRX020661.SRR044768","1927.SRS011157.SRX020661.SRR044808","1927.SRS011415.SRX020659.SRR045013","1927.SRS011586.SRX020659.SRR045050","1927.SRS018920.SRX020572.SRR045810","1927.SRS043804.SRX020531.SRR046050","1927.SRS052326.SRX020670.SRR046232","1927.SRS054928.SRX020670.SRR046262","1927.SRS065263.SRX020670.SRR046351","1927.SRS055563.SRX020523.SRR046412","1927.SRS042387.SRX020523.SRR046416","1927.SRS013543.SRX020677.SRR047980","1927.SRS064557.SRX020569.SRR046542","1927.SRS063961.SRX020532.SRR046907","1927.SRS062464.SRX020532.SRR046982","1927.SRS063214.SRX020532.SRR046986","1927.SRS051031.SRX020578.SRR047244","1927.SRS042669.SRX020540.SRR047279","1927.SRS058416.SRX020523.SRR047701","1927.SRS055665.SRX020523.SRR047702","1927.SRS049748.SRX020523.SRR047703","1927.SRS021541.SRX020675.SRR047797","1927.SRS013762.SRX020677.SRR048007","1927.SRS019381.SRX020572.SRR051585","1927.SRS021109.SRX022097.SRR057663","1927.SRS020584.SRX022097.SRR057663","1927.SRS020641.SRX022097.SRR057663","1927.SRS042803.SRX022156.SRR058115","1927.SRS063918.SRX022232.SRR058094","1927.SRS064647.SRX022232.SRR058094","1927.SRS049766.SRX022241.SRR058097","1927.SRS056832.SRX022241.SRR058097","1927.SRS063608.SRX022156.SRR058115","1927.SRS011621.SRX020659.SRR045004","1927.SRS011405.SRX020659.SRR046237","1927.SRS015518.SRX020539.SRR044146","1927.SRS011405.SRX020661.SRR044803","1927.SRS018968.SRX020572.SRR045817","1927.SRS063921.SRX020569.SRR046605","1927.SRS023422.SRX020683.SRR048341","1927.SRS019534.SRX020572.SRR051591","1927.SRS021853.SRX022232.SRR058094","1927.SRS065421.SRX022232.SRR058094","1927.SRS063307.SRX022232.SRR058094","1927.SRS042572.SRX022241.SRR058097","1927.SRS049949.SRX022241.SRR058097","1927.SRS057676.SRX022156.SRR058115","1927.SRS014923.SRX020546.SRR043675","1927.SRS051892.SRX022156.SRR058115","1927.SRS048838.SRX020523.SRR046418","1927.SRS057122.SRX020677.SRR047969","1927.SRS014345.SRX020546.SRR047368","1927.SRS023851.SRX020676.SRR047925","1927.SRS046341.SRX020554.SRR043828","1927.SRS023422.SRX020676.SRR047874","1927.SRS020811.SRX022097.SRR057663"],[0.0465431540894713,0.401378326996198,0.282029950083195,0.392269280562234,0.567991933160472,0.4,0.40018176310209,0.249308118081181,0.464148033924441,0.36194221234399,0.260434569113675,0.233664772727273,0.421884317233154,0.0711857203569911,0.172183047920564,0.255957634598411,0.146615430091635,0.102931870913309,0.135868120180803,0.00842696629213483,0.24501708428246,0.20120572720422,0.3818634313268,0.473793512658228,0.411666357049972,0.361031518624642,0.367285499247366,0.234413407821229,0.456430666108085,0.278208278208278,0.280252764612954,0.525851938895417,0.341131534972547,0.38109756097561,0.323096129837703,0.0810419681620839,0.316417212347989,0.230706742485784,0.429861529199278,0.469484629294756,0.218311213694941,0.354477180820335,0.403169014084507,0.224445812807882,0.405931242916509,0.600860178782257,0.273193234238852,0.517786069651741,0.39090523937297,0.362703467799009,0.315261789825473,0.714163090128755,0.323939881910896,0.247454751131222,0.389038031319911,0.367110887603721,0.265678449258837,0.70744935396395,0.491539970573811,0.22855773838045,0.477044694435999,0.487407407407407,0.323519682434668,0.210704862419902,0.403456048084147,0.11132712575463,0.47650827549386,0.52166776099803,0.374000761324705,0.219600725952813,0.326815642458101,0.187238723872387,0.441889763779528,0.445383912248629,0.588142178730753,0.415516932449382,0.42611060743427,0.28087044534413,0.576525821596244,0.157524322086525,0.234841035725991,0.26274328081557,0.158908822806586,0.388976377952756,0.671440606571188,0.296296296296296,0.302752293577982,0.225370418507928,0.258557660352277,0.252818585582508,0.270007124198528,0.270133936336754,0.396285979572888,0.527178404186943,0.0671641791044776,0.139395338102931,0.217722534081796,0.157066189624329,0.179702048417132,0.187580437580438,0.0859349855165755,0.24083269671505,0.403130613292276,0.799706529713866,0.332867458653622,0.421744767549546,0.00679953106682298,0.158000574547544,0.151909150757077,0.374622356495468,0.182673942701228,0.173104850421144,0.381830297645491,0.102909713307659,0.166075650118203,0.184064237183447,0.446985446985447,0.118605315804516,0.45103598691385,0.787892113749634,0.460692951015532,0.361114745518775,0.463103122043519,0.104176668266923,0.383466304532337,0.584875040571243,0.337055606198724,0.190838741887169,0.113517745302714,0.285714285714286,0.13189750039302,0.458975940777298,0.304897314375987,0.153448275862069,0.252950643776824,0.779561539789401,0.625614339404452,0.54324942791762,0.376998050682261,0.633716475095785,0.281361097003555,0.511503602138043,0.216030534351145,0.608873505349276,0.233232366940718,0.159851301115242,0.386591962905719,0.466666666666667,0.606121537086685,0.445582281495594,0.182204590535562,0.94693787309308,0.661422708618331,0.632988435788192,0.432082794307891,0.432393693263258,0.0233545647558386,0.24816280384398,0.245683930942895,0.439284485724114,0.0395156150414277,0.512964563526361,0.429738101943115,0.109796784956021,0.536155666579017,0.37303970637304,0.173579801623084,0.759483770043019,0.120355411954766]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>core.abundance<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":1},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


## Core visualization  

### Core heatmaps  

This visualization method has been used for instance in [Intestinal microbiome landscaping: insight in community assemblage and implications for microbial modulation strategies](https://academic.oup.com/femsre/article/41/2/182/2979411).  

Note that you can order the taxa on the heatmap with the order.taxa argument.



```r
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(ps1.stool.rel2, plot.type = "heatmap", colours = gray,
    prevalences = prevalences, detections = detections, min.prevalence = .5) +
    xlab("Detection Threshold (Relative Abundance (%))")
print(p.core)    
```

<img src="06-MAW-PV_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
# Same with the viridis color palette
# color-blind friendly and uniform
# options: viridis, magma, plasma, inferno
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
# Also discrete=TRUE versions available
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
print(p.core + scale_fill_viridis())
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill',
## which will replace the existing scale.
```

<img src="06-MAW-PV_files/figure-html/unnamed-chunk-7-2.png" width="672" />

Color change 


```r
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(ps1.stool.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

print(p.core) 
```

<img src="06-MAW-PV_files/figure-html/unnamed-chunk-8-1.png" width="672" />

We have a custom script to format this figure which can be found in `scripts` folder in the `RProject`.  


```r
source("./scripts/plot_core_assist.R")
plot_data <- p.core$data

df_plot <- plot_core_assist(plot_data, tax_table(ps1.stool.rel2), levels = "Family")

p.core$data <- df_plot

plot(p.core + theme(axis.text.y = element_text(face="italic")))
```

<img src="06-MAW-PV_files/figure-html/unnamed-chunk-9-1.png" width="672" />



```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Dutch_Netherlands.1252  LC_CTYPE=Dutch_Netherlands.1252   
## [3] LC_MONETARY=Dutch_Netherlands.1252 LC_NUMERIC=C                      
## [5] LC_TIME=Dutch_Netherlands.1252    
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2     viridis_0.5.1      viridisLite_0.3.0 
##  [4] dplyr_0.7.4        ggpubr_0.1.6       magrittr_1.5      
##  [7] RColorBrewer_1.1-2 microbiome_1.0.2   ggplot2_2.2.1     
## [10] phyloseq_1.23.1   
## 
## loaded via a namespace (and not attached):
##  [1] Biobase_2.38.0      tidyr_0.8.0         jsonlite_1.5       
##  [4] splines_3.4.4       foreach_1.4.4       shiny_1.0.5        
##  [7] assertthat_0.2.0    stats4_3.4.4        yaml_2.1.18        
## [10] pillar_1.2.2        backports_1.1.2     lattice_0.20-35    
## [13] glue_1.2.0          digest_0.6.15       promises_1.0.1     
## [16] XVector_0.18.0      colorspace_1.3-2    htmltools_0.3.6    
## [19] httpuv_1.4.1        Matrix_1.2-14       plyr_1.8.4         
## [22] pkgconfig_2.0.1     bookdown_0.7        zlibbioc_1.24.0    
## [25] purrr_0.2.4         xtable_1.8-2        scales_0.5.0       
## [28] later_0.7.1         tibble_1.4.2        mgcv_1.8-23        
## [31] IRanges_2.12.0      DT_0.4              BiocGenerics_0.24.0
## [34] lazyeval_0.2.1      survival_2.42-3     mime_0.5           
## [37] evaluate_0.10.1     nlme_3.1-137        MASS_7.3-49        
## [40] vegan_2.5-1         tools_3.4.4         data.table_1.10.4-3
## [43] stringr_1.3.0       S4Vectors_0.16.0    munsell_0.4.3      
## [46] cluster_2.0.7-1     Biostrings_2.46.0   ade4_1.7-11        
## [49] compiler_3.4.4      rlang_0.2.0         rhdf5_2.22.0       
## [52] grid_3.4.4          iterators_1.0.9     biomformat_1.7.0   
## [55] htmlwidgets_1.2     crosstalk_1.0.0     igraph_1.2.1       
## [58] labeling_0.3        rmarkdown_1.9       gtable_0.2.0       
## [61] codetools_0.2-15    multtest_2.34.0     reshape2_1.4.3     
## [64] R6_2.2.2            gridExtra_2.3       knitr_1.20         
## [67] bindr_0.1.1         rprojroot_1.3-2     permute_0.9-4      
## [70] ape_5.1             stringi_1.1.7       parallel_3.4.4     
## [73] Rcpp_0.12.16        tidyselect_0.2.4    xfun_0.1
```



