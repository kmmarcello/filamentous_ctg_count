---
title: "Filamentous CTG Count Heatmap"
author: "Kaylah Marcello"
date: '2023-01-19'
output: 
  html_document: 
    keep_md: yes
---




```r
#install.packages("viridisLite")
#install.packages("svglite")
#install.packages("factoextra")
#install.packages("cluster")
#install.packages("ggdendro")
#install.packages("grid")
```

## Load the Libraries


```r
library(tidyverse)
library(RColorBrewer)
library(paletteer)
library(janitor)
library(here)
library(skimr)
#library(viridis)
#library(shiny)
#library(shinydashboard)
library(devtools)
library(svglite)
#library(cluster)
#library(factoextra)
library(ggdendro)
library(grid)
library(gplots)
```


```r
filament_cyanos <- readr::read_csv("data/filamentous_cyano_meta.csv")
```
## Get an Idea of the Structure

```r
# summary(filament_cyanos)
```


```r
head(filament_cyanos)
```

```
## # A tibble: 6 × 40
##   GenBank…¹ RefSe…² Acces…³ Organ…⁴ Genus Species Strain fila_…⁵ envir…⁶ Stora…⁷
##   <chr>     <chr>   <chr>   <chr>   <chr> <chr>   <chr>  <chr>   <chr>   <chr>  
## 1 GCA_0000… GCF_00… ASM200… Nostoc… Nost… puncti… PCC 7… filame… <NA>    ATCC 2…
## 2 GCA_0003… GCF_00… ASM317… Oscill… Osci… acumin… PCC 6… filame… <NA>    JGI    
## 3 GCA_0003… GCF_00… ASM317… Oscill… Osci… nigro-… <NA>   filame… <NA>    <NA>   
## 4 GCA_0125… GCF_01… ASM125… Limnos… Limn… fusifo… SAG 8… filame… x_hot   <NA>   
## 5 GCA_0028… GCF_00… ASM281… Nostoc… Nost… flagel… CCNUN1 filame… <NA>    <NA>   
## 6 GCA_0000… GCF_00… ASM142… Tricho… Tric… erythr… IMS101 filame… <NA>    <NA>   
## # … with 30 more variables: `Regional Loaction` <chr>,
## #   `geographic feature` <chr>, `Environment Detail` <chr>,
## #   `Tempurature (avg)` <dbl>, `Lat/Long` <chr>, `Collection date` <chr>,
## #   `Storage/Collection...17` <chr>, genome_or_bin <chr>, gene_gyrA <dbl>,
## #   gene_nusA <dbl>, gene_infC <dbl>, gene_infA <dbl>, gene_otsA <dbl>,
## #   gene_dnaK <dbl>, gene_recA <dbl>, gene_dnaJ <dbl>, gene_aceF <dbl>,
## #   gene_deaD <dbl>, gene_infB <dbl>, gene_tig <dbl>, gene_rnr <dbl>, …
```


```r
filament_cyanos <- clean_names(filament_cyanos)
```


```r
names(filament_cyanos)
```

```
##  [1] "gen_bank_assembly_id_accession_version"
##  [2] "ref_seq_assembly_id_accession_version" 
##  [3] "accession_number"                      
##  [4] "organism"                              
##  [5] "genus"                                 
##  [6] "species"                               
##  [7] "strain"                                
##  [8] "fila_single"                           
##  [9] "environment"                           
## [10] "storage_collection_10"                 
## [11] "regional_loaction"                     
## [12] "geographic_feature"                    
## [13] "environment_detail"                    
## [14] "tempurature_avg"                       
## [15] "lat_long"                              
## [16] "collection_date"                       
## [17] "storage_collection_17"                 
## [18] "genome_or_bin"                         
## [19] "gene_gyr_a"                            
## [20] "gene_nus_a"                            
## [21] "gene_inf_c"                            
## [22] "gene_inf_a"                            
## [23] "gene_ots_a"                            
## [24] "gene_dna_k"                            
## [25] "gene_rec_a"                            
## [26] "gene_dna_j"                            
## [27] "gene_ace_f"                            
## [28] "gene_dea_d"                            
## [29] "gene_inf_b"                            
## [30] "gene_tig"                              
## [31] "gene_rnr"                              
## [32] "gene_dna_a"                            
## [33] "gene_hup_b"                            
## [34] "gene_rbf_a"                            
## [35] "gene_yfl_a"                            
## [36] "gene_pnp"                              
## [37] "gene_csp"                              
## [38] "gene_ace_e"                            
## [39] "gene_des_a"                            
## [40] "temp_source"
```
## Select needed columns (values) from metadata


```r
gene_data_organism <- filament_cyanos %>% 
  select(organism, contains("gene_")) %>% 
  pivot_longer(-organism,
               names_to = "gene",
               values_to = "gene_count") %>% 
  filter(!is.na(gene_count))
gene_data_organism
```

```
## # A tibble: 861 × 3
##    organism                     gene       gene_count
##    <chr>                        <chr>           <dbl>
##  1 Nostoc punctiforme PCC 73102 gene_gyr_a          2
##  2 Nostoc punctiforme PCC 73102 gene_nus_a          1
##  3 Nostoc punctiforme PCC 73102 gene_inf_c          5
##  4 Nostoc punctiforme PCC 73102 gene_inf_a          1
##  5 Nostoc punctiforme PCC 73102 gene_ots_a          0
##  6 Nostoc punctiforme PCC 73102 gene_dna_k         10
##  7 Nostoc punctiforme PCC 73102 gene_rec_a          1
##  8 Nostoc punctiforme PCC 73102 gene_dna_j         12
##  9 Nostoc punctiforme PCC 73102 gene_ace_f          3
## 10 Nostoc punctiforme PCC 73102 gene_dea_d         19
## # … with 851 more rows
```

## Heatmap


```r
gene_data_organism$gene_count <- as.numeric(gene_data_organism$gene_count) # needs to be numeric for scale()
lapply(gene_data_organism, class)
```

```
## $organism
## [1] "character"
## 
## $gene
## [1] "character"
## 
## $gene_count
## [1] "numeric"
```


```r
clust_map_org <- gene_data_organism %>% # make pivot_wider for scale()
  pivot_wider(names_from = "gene",
              values_from = "gene_count")
  
names(clust_map_org) = gsub(pattern = "gene_", replacement = "", x = names(clust_map_org))
names(clust_map_org) = gsub(pattern = "_", replacement = " ", x = names(clust_map_org))

clust_map_org
```

```
## # A tibble: 41 × 22
##    organism      `gyr a` `nus a` `inf c` `inf a` `ots a` `dna k` `rec a` `dna j`
##    <chr>           <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
##  1 Nostoc punct…       2       1       5       1       0      10       1      12
##  2 Oscillatoria…       2       2       1       1       0       9       3       7
##  3 Oscillatoria…       2       1       1       1       0       7       2       7
##  4 Nostoc flage…       2       1       2       0       1       7       1      10
##  5 Trichodesmiu…       2       1       1       1       0      11       1       7
##  6 Nostoc sphae…       2       1       2       1       0       9       1      10
##  7 Leptolyngbya…       2       1       1       1       0       5       2       7
##  8 Nostoc azoll…       2       1       1       1       0       8       1       7
##  9 Dolichosperm…       2       1       1       1       0       7       1       6
## 10 Anabaena cyl…       2       1       1       1       0       8       1       8
## # … with 31 more rows, and 13 more variables: `ace f` <dbl>, `dea d` <dbl>,
## #   `inf b` <dbl>, tig <dbl>, rnr <dbl>, `dna a` <dbl>, `hup b` <dbl>,
## #   `rbf a` <dbl>, `yfl a` <dbl>, pnp <dbl>, csp <dbl>, `ace e` <dbl>,
## #   `des a` <dbl>
```


```r
scaled_clust_map <- clust_map_org 
scaled_clust_map[, c(2:22)] <- scale(scaled_clust_map[, 2:22])

scaled_clust_map
```

```
## # A tibble: 41 × 22
##    organism      `gyr a` `nus a` `inf c` `inf a` `ots a` `dna k` `rec a` `dna j`
##    <chr>           <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
##  1 Nostoc punct…  -0.156  -0.156   4.24   0.0894  -0.448   1.70   -0.425   1.89 
##  2 Oscillatoria…  -0.156   6.25   -0.727  0.0894  -0.448   1.00    3.93   -0.937
##  3 Oscillatoria…  -0.156  -0.156  -0.727  0.0894  -0.448  -0.392   1.75   -0.937
##  4 Nostoc flage…  -0.156  -0.156   0.515 -3.58     2.18   -0.392  -0.425   0.758
##  5 Trichodesmiu…  -0.156  -0.156  -0.727  0.0894  -0.448   2.40   -0.425  -0.937
##  6 Nostoc sphae…  -0.156  -0.156   0.515  0.0894  -0.448   1.00   -0.425   0.758
##  7 Leptolyngbya…  -0.156  -0.156  -0.727  0.0894  -0.448  -1.79    1.75   -0.937
##  8 Nostoc azoll…  -0.156  -0.156  -0.727  0.0894  -0.448   0.306  -0.425  -0.937
##  9 Dolichosperm…  -0.156  -0.156  -0.727  0.0894  -0.448  -0.392  -0.425  -1.50 
## 10 Anabaena cyl…  -0.156  -0.156  -0.727  0.0894  -0.448   0.306  -0.425  -0.372
## # … with 31 more rows, and 13 more variables: `ace f` <dbl>, `dea d` <dbl>,
## #   `inf b` <dbl>, tig <dbl>, rnr <dbl>, `dna a` <dbl>, `hup b` <dbl>,
## #   `rbf a` <dbl>, `yfl a` <dbl>, pnp <dbl>, csp <dbl>, `ace e` <dbl>,
## #   `des a` <dbl>
```




```r
rnames <- scaled_clust_map[,1] # made singe column to use as rownames
rnames
```

```
## # A tibble: 41 × 1
##    organism                               
##    <chr>                                  
##  1 Nostoc punctiforme PCC 73102           
##  2 Oscillatoria acuminata PCC 6304        
##  3 Oscillatoria nigro-viridis PCC 7112    
##  4 Nostoc flagelliforme CCNUN1            
##  5 Trichodesmium erythraeum IMS101        
##  6 Nostoc sphaeroides                     
##  7 Leptolyngbya boryana NIES-2135         
##  8 Nostoc azollae' 0708                   
##  9 Dolichospermum flos-aquae CCAP 1403/13F
## 10 Anabaena cylindrica PCC 7122           
## # … with 31 more rows
```


```r
organisim <- as.list(rnames) # changed to list so I could add to matrix, did not work so I wrote them all out to make a vactor. I think there is a problem with the title "organism"
organisim
```

```
## $organism
##  [1] "Nostoc punctiforme PCC 73102"                     
##  [2] "Oscillatoria acuminata PCC 6304"                  
##  [3] "Oscillatoria nigro-viridis PCC 7112"              
##  [4] "Nostoc flagelliforme CCNUN1"                      
##  [5] "Trichodesmium erythraeum IMS101"                  
##  [6] "Nostoc sphaeroides"                               
##  [7] "Leptolyngbya boryana NIES-2135"                   
##  [8] "Nostoc azollae' 0708"                             
##  [9] "Dolichospermum flos-aquae CCAP 1403/13F"          
## [10] "Anabaena cylindrica PCC 7122"                     
## [11] "Nostoc sp. 'Peltigera membranacea cyanobiont' N6" 
## [12] "Nostoc sp. TCL240-02"                             
## [13] "Nostoc sp. C052"                                  
## [14] "Nostoc edaphicum CCNP1411"                        
## [15] "Nostoc sp. C057"                                  
## [16] "Microcoleus sp. PCC 7113"                         
## [17] "Nostoc sp. ATCC 53789"                            
## [18] "Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont'"
## [19] "Calothrix sp. PCC 7507"                           
## [20] "Anabaena sp. YBS01"                               
## [21] "Nostoc sp. CENA543"                               
## [22] "Nostoc sp. TCL26-01"                              
## [23] "Nostoc sp. PCC 7524"                              
## [24] "Nostoc sphaeroides CCNUC1"                        
## [25] "Nostoc sp. NIES-3756"                             
## [26] "Nostoc sp. PCC 7120 = FACHB-418"                  
## [27] "Leptolyngbya sp. NIES-3755"                       
## [28] "Leptolyngbya boryana dg5"                         
## [29] "Leptolyngbya boryana IAM M-101"                   
## [30] "Anabaena sp. WA102"                               
## [31] "Thermoleptolyngbya sp. PKUAC-SCTA183"             
## [32] "Leptolyngbya sp. O-77"                            
## [33] "Dolichospermum sp. UHCC 0315A"                    
## [34] "Leptolyngbya sp. PCC 7376"                        
## [35] "Pseudanabaena sp. ABRG5-3"                        
## [36] "Pseudanabaena sp. PCC 7367"                       
## [37] "Anabaena sp. 90"                                  
## [38] "Leptolyngbya D1"                                  
## [39] "Phormidesmis"                                     
## [40] "Leptolyngbya D3"                                  
## [41] "Anabaena"
```

### Convert dataframe to matrix with heatmap.2 from gplots

```r
mat_data <- data.matrix(scaled_clust_map[,2:ncol(scaled_clust_map)])
mat_data[is.na(mat_data)] <- 0
mat_data[,colnames(mat_data)!="organism"]
```

```
##            gyr a      nus a      inf c       inf a     ots a      dna k
##  [1,] -0.1561738 -0.1561738  4.2393111  0.08942484 -0.448175  1.7024761
##  [2,] -0.1561738  6.2469505 -0.7267390  0.08942484 -0.448175  1.0044609
##  [3,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -0.3915695
##  [4,] -0.1561738 -0.1561738  0.5147735 -3.57699344  2.176850 -0.3915695
##  [5,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  2.4004912
##  [6,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  1.0044609
##  [7,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -1.7875999
##  [8,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  0.3064457
##  [9,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -0.3915695
## [10,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  0.3064457
## [11,] -0.1561738 -0.1561738  0.5147735  3.75584311 -0.448175 -0.3915695
## [12,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175 -0.3915695
## [13,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175 -1.0895847
## [14,] -0.1561738 -0.1561738  0.5147735  0.08942484  2.176850  1.0044609
## [15,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175 -1.0895847
## [16,] -0.1561738 -0.1561738  1.7562860  0.08942484 -0.448175  0.3064457
## [17,] -0.1561738 -0.1561738  1.7562860  0.08942484 -0.448175  0.3064457
## [18,] -0.1561738 -0.1561738  0.5147735 -3.57699344 -0.448175  0.3064457
## [19,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  0.3064457
## [20,] -0.1561738 -0.1561738  0.5147735  0.08942484  2.176850  1.0044609
## [21,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  1.7024761
## [22,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175 -0.3915695
## [23,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175 -0.3915695
## [24,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  1.7024761
## [25,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  0.3064457
## [26,] -0.1561738 -0.1561738  0.5147735  0.08942484 -0.448175  0.3064457
## [27,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -1.0895847
## [28,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -1.7875999
## [29,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -1.7875999
## [30,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  1.0044609
## [31,] -0.1561738 -0.1561738 -0.7267390  0.08942484  2.176850 -1.0895847
## [32,]  6.2469505 -0.1561738 -0.7267390  0.08942484  2.176850 -0.3915695
## [33,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -0.3915695
## [34,] -0.1561738 -0.1561738 -0.7267390  0.08942484  2.176850 -1.0895847
## [35,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  0.3064457
## [36,] -0.1561738 -0.1561738 -0.7267390  0.08942484  2.176850 -0.3915695
## [37,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  0.3064457
## [38,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175  1.0044609
## [39,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -0.3915695
## [40,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -0.3915695
## [41,] -0.1561738 -0.1561738 -0.7267390  0.08942484 -0.448175 -1.0895847
##            rec a      dna j      ace f       dea d      inf b tig        rnr
##  [1,] -0.4248054  1.8885596  4.1411822  1.14599909 -0.5588104   0 -0.2605158
##  [2,]  3.9294497 -0.9373873 -0.3865103 -1.05217460 -0.5588104   0 -0.2605158
##  [3,]  1.7523222 -0.9373873 -0.3865103  0.59645567  1.2035916   0 -0.2605158
##  [4,] -0.4248054  0.7581809 -0.3865103  2.51985765  1.2035916   0 -0.2605158
##  [5,] -0.4248054 -0.9373873 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
##  [6,] -0.4248054  0.7581809 -0.3865103  0.04691224 -0.5588104   0 -0.2605158
##  [7,]  1.7523222 -0.9373873  1.8773359 -0.77740289 -0.5588104   0 -0.2605158
##  [8,] -0.4248054 -0.9373873 -0.3865103 -1.60171803 -0.5588104   0 -0.2605158
##  [9,] -0.4248054 -1.5025766 -0.3865103 -0.77740289  1.2035916   0 -0.2605158
## [10,] -0.4248054 -0.3721979 -0.3865103  0.32168396 -0.5588104   0 -0.2605158
## [11,] -0.4248054  0.7581809 -0.3865103  2.24508594  1.2035916   0 -0.2605158
## [12,] -0.4248054  0.1929915  1.8773359 -0.22785947 -0.5588104   0 -0.2605158
## [13,] -0.4248054  1.8885596 -0.3865103  1.14599909  1.2035916   0 -0.2605158
## [14,] -0.4248054  0.7581809 -0.3865103 -0.50263118  1.2035916   0 -0.2605158
## [15,] -0.4248054  1.8885596 -0.3865103  1.14599909  1.2035916   0 -0.2605158
## [16,]  1.7523222  0.7581809 -0.3865103  0.32168396  1.2035916   0 -0.2605158
## [17,] -0.4248054  1.3233703  1.8773359  0.32168396 -0.5588104   0 -0.2605158
## [18,] -0.4248054  0.1929915 -0.3865103  0.87122738 -0.5588104   0 -0.2605158
## [19,] -0.4248054  0.1929915 -0.3865103  0.32168396  1.2035916   0 -0.2605158
## [20,] -0.4248054  0.7581809 -0.3865103  0.32168396 -0.5588104   0 -0.2605158
## [21,] -0.4248054  0.1929915 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
## [22,] -0.4248054  0.1929915 -0.3865103 -0.50263118 -0.5588104   0 -0.2605158
## [23,]  1.7523222  1.3233703 -0.3865103 -0.77740289 -0.5588104   0 -0.2605158
## [24,] -0.4248054  0.7581809 -0.3865103  2.51985765 -0.5588104   0 -0.2605158
## [25,] -0.4248054  0.1929915 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
## [26,] -0.4248054  0.7581809 -0.3865103  0.32168396 -0.5588104   0 -0.2605158
## [27,] -0.4248054 -0.3721979 -0.3865103 -0.22785947  2.9659935   0 -0.2605158
## [28,]  1.7523222 -0.9373873  1.8773359 -0.77740289 -0.5588104   0 -0.2605158
## [29,]  1.7523222 -0.9373873  1.8773359 -0.77740289 -0.5588104   0 -0.2605158
## [30,] -0.4248054 -2.0677660 -0.3865103 -0.77740289 -0.5588104   0 -0.2605158
## [31,] -0.4248054 -0.3721979 -0.3865103  0.04691224 -0.5588104   0 -0.2605158
## [32,] -0.4248054 -0.3721979 -0.3865103 -0.22785947 -0.5588104   0  2.4097716
## [33,] -0.4248054 -0.3721979 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
## [34,] -0.4248054  0.1929915 -0.3865103 -1.60171803  2.9659935   0  2.4097716
## [35,] -0.4248054 -1.5025766 -0.3865103  0.87122738 -0.5588104   0 -0.2605158
## [36,] -0.4248054 -1.5025766 -0.3865103 -1.87648974 -0.5588104   0 -0.2605158
## [37,] -0.4248054 -1.5025766 -0.3865103 -1.05217460 -0.5588104   0 -0.2605158
## [38,] -0.4248054  0.7581809 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
## [39,] -0.4248054 -0.3721979 -0.3865103  0.32168396 -0.5588104   0  5.0800590
## [40,] -0.4248054  0.1929915 -0.3865103 -0.22785947 -0.5588104   0 -0.2605158
## [41,] -0.4248054  0.1929915 -0.3865103 -0.50263118 -0.5588104   0 -0.2605158
##             dna a       hup b      rbf a      yfl a        pnp        csp
##  [1,] -0.58034566  1.48022191 -0.3247635 -0.2116563  1.8411606 -0.8583234
##  [2,] -0.03956902 -0.88428842 -0.3247635  2.6809803  1.8411606  0.8174508
##  [3,] -0.03956902 -0.88428842 -0.3247635 -0.2116563 -0.1990444  0.8174508
##  [4,] -0.58034566  0.69205180 -0.3247635 -0.2116563 -0.1990444 -0.8583234
##  [5,]  0.50120762 -0.88428842 -0.3247635  5.5736170 -0.1990444 -0.8583234
##  [6,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
##  [7,]  2.12353753  0.69205180 -0.3247635 -0.2116563 -0.1990444  0.8174508
##  [8,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
##  [9,] -0.03956902 -0.88428842  3.0040623 -0.2116563 -0.1990444  0.8174508
## [10,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [11,] -0.58034566  2.26839202 -0.3247635 -0.2116563  3.8813657 -0.8583234
## [12,] -0.58034566  0.69205180 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [13,] -0.03956902  1.48022191 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [14,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [15,] -0.03956902  1.48022191 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [16,] -0.03956902 -0.88428842 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [17,] -0.58034566  0.69205180 -0.3247635 -0.2116563  1.8411606  0.8174508
## [18,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [19,]  0.50120762 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [20,] -0.03956902 -0.09611831 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [21,] -0.03956902 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [22,] -0.58034566 -0.09611831 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [23,] -0.58034566 -0.09611831 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [24,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [25,] -0.03956902 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [26,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [27,] -0.03956902 -0.09611831 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [28,]  1.58276089  0.69205180 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [29,]  1.58276089  0.69205180 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [30,] -0.58034566 -0.09611831  3.0040623 -0.2116563  1.8411606  2.4932250
## [31,] -0.03956902 -0.88428842 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [32,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [33,] -0.03956902 -0.09611831  3.0040623 -0.2116563 -0.1990444  0.8174508
## [34,] -0.58034566 -0.88428842 -0.3247635 -0.2116563 -2.2392494  0.8174508
## [35,] -0.03956902  0.69205180 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [36,]  4.82742072  0.69205180 -0.3247635 -0.2116563 -2.2392494 -0.8583234
## [37,] -0.03956902 -0.09611831 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [38,] -0.58034566  2.26839202 -0.3247635 -0.2116563 -0.1990444  0.8174508
## [39,]  0.50120762 -0.09611831 -0.3247635 -0.2116563 -0.1990444  2.4932250
## [40,] -0.58034566  2.26839202 -0.3247635 -0.2116563 -0.1990444 -0.8583234
## [41,] -0.03956902 -0.09611831  3.0040623 -0.2116563 -0.1990444  0.8174508
##            ace e       des a
##  [1,] -0.4621465  1.39967730
##  [2,] -2.5674808  0.08547648
##  [3,] -0.4621465 -1.22872435
##  [4,] -0.4621465  1.83774424
##  [5,] -0.4621465 -1.22872435
##  [6,] -0.4621465 -0.79065741
##  [7,] -0.4621465  0.96161036
##  [8,]  1.6431877 -1.66679129
##  [9,] -0.4621465 -0.35259046
## [10,]  1.6431877 -0.79065741
## [11,] -0.4621465 -0.79065741
## [12,]  1.6431877  1.83774424
## [13,] -0.4621465  1.83774424
## [14,] -0.4621465  0.96161036
## [15,] -0.4621465  1.83774424
## [16,] -0.4621465 -0.79065741
## [17,] -0.4621465  0.08547648
## [18,]  1.6431877  0.08547648
## [19,]  1.6431877 -0.79065741
## [20,]  1.6431877 -0.35259046
## [21,]  1.6431877 -0.35259046
## [22,]  1.6431877 -0.79065741
## [23,]  1.6431877 -1.22872435
## [24,] -0.4621465  0.08547648
## [25,] -0.4621465  0.08547648
## [26,]  1.6431877 -0.35259046
## [27,] -0.4621465  0.52354342
## [28,] -0.4621465  0.96161036
## [29,] -0.4621465  0.96161036
## [30,] -0.4621465  0.52354342
## [31,] -0.4621465 -0.79065741
## [32,] -0.4621465 -0.79065741
## [33,] -0.4621465 -0.79065741
## [34,] -0.4621465  0.08547648
## [35,] -0.4621465  0.52354342
## [36,] -0.4621465 -1.66679129
## [37,] -0.4621465 -0.79065741
## [38,] -0.4621465  0.52354342
## [39,] -0.4621465  0.52354342
## [40,] -0.4621465  1.39967730
## [41,] -0.4621465 -0.79065741
```

```r
#rownames(mat_data) <- rnames   # this does not work, gives error: 
```


```r
rownames(mat_data) <- c("Nostoc punctiforme PCC 73102", "Oscillatoria acuminata PCC 6304", "Oscillatoria nigro-viridis PCC 7112", "Nostoc flagelliforme CCNUN1", "Trichodesmium erythraeum IMS101", "Nostoc sphaeroides", "Leptolyngbya boryana NIES-2135", "Nostoc azollae' 0708", "Dolichospermum flos-aquae CCAP 1403/13F", "Anabaena cylindrica PCC 7122", "Nostoc sp. 'Peltigera membranacea cyanobiont' N6", "Nostoc sp. TCL240-02", "Nostoc sp. C052", "Nostoc edaphicum CCNP1411",  "Nostoc sp. C057", "Microcoleus sp. PCC 7113", "Nostoc sp. ATCC 53789",  "Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont'", "Calothrix sp. PCC 7507",  "Anabaena sp. YBS01",  "Nostoc sp. CENA543", "Nostoc sp. TCL26-01",  "Nostoc sp. PCC 7524", "Nostoc sphaeroides CCNUC1", "Nostoc sp. NIES-3756", "Nostoc sp. PCC 7120 = FACHB-418", "Leptolyngbya sp. NIES-3755","Leptolyngbya boryana dg5",  "Leptolyngbya boryana IAM M-101",  "Anabaena sp. WA102", "Thermoleptolyngbya sp. PKUAC-SCTA183", "Leptolyngbya sp. O-77",  "Dolichospermum sp. UHCC 0315A", "Leptolyngbya sp. PCC 7376", "Pseudanabaena sp. ABRG5-3", "Pseudanabaena sp. PCC 7367", "Anabaena sp. 90", "Leptolyngbya D1", "Phormidesmis", "Leptolyngbya D3", "Anabaena")

mat_data
```

```
##                                                        gyr a      nus a
## Nostoc punctiforme PCC 73102                      -0.1561738 -0.1561738
## Oscillatoria acuminata PCC 6304                   -0.1561738  6.2469505
## Oscillatoria nigro-viridis PCC 7112               -0.1561738 -0.1561738
## Nostoc flagelliforme CCNUN1                       -0.1561738 -0.1561738
## Trichodesmium erythraeum IMS101                   -0.1561738 -0.1561738
## Nostoc sphaeroides                                -0.1561738 -0.1561738
## Leptolyngbya boryana NIES-2135                    -0.1561738 -0.1561738
## Nostoc azollae' 0708                              -0.1561738 -0.1561738
## Dolichospermum flos-aquae CCAP 1403/13F           -0.1561738 -0.1561738
## Anabaena cylindrica PCC 7122                      -0.1561738 -0.1561738
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.1561738 -0.1561738
## Nostoc sp. TCL240-02                              -0.1561738 -0.1561738
## Nostoc sp. C052                                   -0.1561738 -0.1561738
## Nostoc edaphicum CCNP1411                         -0.1561738 -0.1561738
## Nostoc sp. C057                                   -0.1561738 -0.1561738
## Microcoleus sp. PCC 7113                          -0.1561738 -0.1561738
## Nostoc sp. ATCC 53789                             -0.1561738 -0.1561738
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.1561738 -0.1561738
## Calothrix sp. PCC 7507                            -0.1561738 -0.1561738
## Anabaena sp. YBS01                                -0.1561738 -0.1561738
## Nostoc sp. CENA543                                -0.1561738 -0.1561738
## Nostoc sp. TCL26-01                               -0.1561738 -0.1561738
## Nostoc sp. PCC 7524                               -0.1561738 -0.1561738
## Nostoc sphaeroides CCNUC1                         -0.1561738 -0.1561738
## Nostoc sp. NIES-3756                              -0.1561738 -0.1561738
## Nostoc sp. PCC 7120 = FACHB-418                   -0.1561738 -0.1561738
## Leptolyngbya sp. NIES-3755                        -0.1561738 -0.1561738
## Leptolyngbya boryana dg5                          -0.1561738 -0.1561738
## Leptolyngbya boryana IAM M-101                    -0.1561738 -0.1561738
## Anabaena sp. WA102                                -0.1561738 -0.1561738
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.1561738 -0.1561738
## Leptolyngbya sp. O-77                              6.2469505 -0.1561738
## Dolichospermum sp. UHCC 0315A                     -0.1561738 -0.1561738
## Leptolyngbya sp. PCC 7376                         -0.1561738 -0.1561738
## Pseudanabaena sp. ABRG5-3                         -0.1561738 -0.1561738
## Pseudanabaena sp. PCC 7367                        -0.1561738 -0.1561738
## Anabaena sp. 90                                   -0.1561738 -0.1561738
## Leptolyngbya D1                                   -0.1561738 -0.1561738
## Phormidesmis                                      -0.1561738 -0.1561738
## Leptolyngbya D3                                   -0.1561738 -0.1561738
## Anabaena                                          -0.1561738 -0.1561738
##                                                        inf c       inf a
## Nostoc punctiforme PCC 73102                       4.2393111  0.08942484
## Oscillatoria acuminata PCC 6304                   -0.7267390  0.08942484
## Oscillatoria nigro-viridis PCC 7112               -0.7267390  0.08942484
## Nostoc flagelliforme CCNUN1                        0.5147735 -3.57699344
## Trichodesmium erythraeum IMS101                   -0.7267390  0.08942484
## Nostoc sphaeroides                                 0.5147735  0.08942484
## Leptolyngbya boryana NIES-2135                    -0.7267390  0.08942484
## Nostoc azollae' 0708                              -0.7267390  0.08942484
## Dolichospermum flos-aquae CCAP 1403/13F           -0.7267390  0.08942484
## Anabaena cylindrica PCC 7122                      -0.7267390  0.08942484
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6   0.5147735  3.75584311
## Nostoc sp. TCL240-02                               0.5147735  0.08942484
## Nostoc sp. C052                                    0.5147735  0.08942484
## Nostoc edaphicum CCNP1411                          0.5147735  0.08942484
## Nostoc sp. C057                                    0.5147735  0.08942484
## Microcoleus sp. PCC 7113                           1.7562860  0.08942484
## Nostoc sp. ATCC 53789                              1.7562860  0.08942484
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont'  0.5147735 -3.57699344
## Calothrix sp. PCC 7507                             0.5147735  0.08942484
## Anabaena sp. YBS01                                 0.5147735  0.08942484
## Nostoc sp. CENA543                                 0.5147735  0.08942484
## Nostoc sp. TCL26-01                                0.5147735  0.08942484
## Nostoc sp. PCC 7524                                0.5147735  0.08942484
## Nostoc sphaeroides CCNUC1                          0.5147735  0.08942484
## Nostoc sp. NIES-3756                               0.5147735  0.08942484
## Nostoc sp. PCC 7120 = FACHB-418                    0.5147735  0.08942484
## Leptolyngbya sp. NIES-3755                        -0.7267390  0.08942484
## Leptolyngbya boryana dg5                          -0.7267390  0.08942484
## Leptolyngbya boryana IAM M-101                    -0.7267390  0.08942484
## Anabaena sp. WA102                                -0.7267390  0.08942484
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.7267390  0.08942484
## Leptolyngbya sp. O-77                             -0.7267390  0.08942484
## Dolichospermum sp. UHCC 0315A                     -0.7267390  0.08942484
## Leptolyngbya sp. PCC 7376                         -0.7267390  0.08942484
## Pseudanabaena sp. ABRG5-3                         -0.7267390  0.08942484
## Pseudanabaena sp. PCC 7367                        -0.7267390  0.08942484
## Anabaena sp. 90                                   -0.7267390  0.08942484
## Leptolyngbya D1                                   -0.7267390  0.08942484
## Phormidesmis                                      -0.7267390  0.08942484
## Leptolyngbya D3                                   -0.7267390  0.08942484
## Anabaena                                          -0.7267390  0.08942484
##                                                       ots a      dna k
## Nostoc punctiforme PCC 73102                      -0.448175  1.7024761
## Oscillatoria acuminata PCC 6304                   -0.448175  1.0044609
## Oscillatoria nigro-viridis PCC 7112               -0.448175 -0.3915695
## Nostoc flagelliforme CCNUN1                        2.176850 -0.3915695
## Trichodesmium erythraeum IMS101                   -0.448175  2.4004912
## Nostoc sphaeroides                                -0.448175  1.0044609
## Leptolyngbya boryana NIES-2135                    -0.448175 -1.7875999
## Nostoc azollae' 0708                              -0.448175  0.3064457
## Dolichospermum flos-aquae CCAP 1403/13F           -0.448175 -0.3915695
## Anabaena cylindrica PCC 7122                      -0.448175  0.3064457
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.448175 -0.3915695
## Nostoc sp. TCL240-02                              -0.448175 -0.3915695
## Nostoc sp. C052                                   -0.448175 -1.0895847
## Nostoc edaphicum CCNP1411                          2.176850  1.0044609
## Nostoc sp. C057                                   -0.448175 -1.0895847
## Microcoleus sp. PCC 7113                          -0.448175  0.3064457
## Nostoc sp. ATCC 53789                             -0.448175  0.3064457
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.448175  0.3064457
## Calothrix sp. PCC 7507                            -0.448175  0.3064457
## Anabaena sp. YBS01                                 2.176850  1.0044609
## Nostoc sp. CENA543                                -0.448175  1.7024761
## Nostoc sp. TCL26-01                               -0.448175 -0.3915695
## Nostoc sp. PCC 7524                               -0.448175 -0.3915695
## Nostoc sphaeroides CCNUC1                         -0.448175  1.7024761
## Nostoc sp. NIES-3756                              -0.448175  0.3064457
## Nostoc sp. PCC 7120 = FACHB-418                   -0.448175  0.3064457
## Leptolyngbya sp. NIES-3755                        -0.448175 -1.0895847
## Leptolyngbya boryana dg5                          -0.448175 -1.7875999
## Leptolyngbya boryana IAM M-101                    -0.448175 -1.7875999
## Anabaena sp. WA102                                -0.448175  1.0044609
## Thermoleptolyngbya sp. PKUAC-SCTA183               2.176850 -1.0895847
## Leptolyngbya sp. O-77                              2.176850 -0.3915695
## Dolichospermum sp. UHCC 0315A                     -0.448175 -0.3915695
## Leptolyngbya sp. PCC 7376                          2.176850 -1.0895847
## Pseudanabaena sp. ABRG5-3                         -0.448175  0.3064457
## Pseudanabaena sp. PCC 7367                         2.176850 -0.3915695
## Anabaena sp. 90                                   -0.448175  0.3064457
## Leptolyngbya D1                                   -0.448175  1.0044609
## Phormidesmis                                      -0.448175 -0.3915695
## Leptolyngbya D3                                   -0.448175 -0.3915695
## Anabaena                                          -0.448175 -1.0895847
##                                                        rec a      dna j
## Nostoc punctiforme PCC 73102                      -0.4248054  1.8885596
## Oscillatoria acuminata PCC 6304                    3.9294497 -0.9373873
## Oscillatoria nigro-viridis PCC 7112                1.7523222 -0.9373873
## Nostoc flagelliforme CCNUN1                       -0.4248054  0.7581809
## Trichodesmium erythraeum IMS101                   -0.4248054 -0.9373873
## Nostoc sphaeroides                                -0.4248054  0.7581809
## Leptolyngbya boryana NIES-2135                     1.7523222 -0.9373873
## Nostoc azollae' 0708                              -0.4248054 -0.9373873
## Dolichospermum flos-aquae CCAP 1403/13F           -0.4248054 -1.5025766
## Anabaena cylindrica PCC 7122                      -0.4248054 -0.3721979
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.4248054  0.7581809
## Nostoc sp. TCL240-02                              -0.4248054  0.1929915
## Nostoc sp. C052                                   -0.4248054  1.8885596
## Nostoc edaphicum CCNP1411                         -0.4248054  0.7581809
## Nostoc sp. C057                                   -0.4248054  1.8885596
## Microcoleus sp. PCC 7113                           1.7523222  0.7581809
## Nostoc sp. ATCC 53789                             -0.4248054  1.3233703
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.4248054  0.1929915
## Calothrix sp. PCC 7507                            -0.4248054  0.1929915
## Anabaena sp. YBS01                                -0.4248054  0.7581809
## Nostoc sp. CENA543                                -0.4248054  0.1929915
## Nostoc sp. TCL26-01                               -0.4248054  0.1929915
## Nostoc sp. PCC 7524                                1.7523222  1.3233703
## Nostoc sphaeroides CCNUC1                         -0.4248054  0.7581809
## Nostoc sp. NIES-3756                              -0.4248054  0.1929915
## Nostoc sp. PCC 7120 = FACHB-418                   -0.4248054  0.7581809
## Leptolyngbya sp. NIES-3755                        -0.4248054 -0.3721979
## Leptolyngbya boryana dg5                           1.7523222 -0.9373873
## Leptolyngbya boryana IAM M-101                     1.7523222 -0.9373873
## Anabaena sp. WA102                                -0.4248054 -2.0677660
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.4248054 -0.3721979
## Leptolyngbya sp. O-77                             -0.4248054 -0.3721979
## Dolichospermum sp. UHCC 0315A                     -0.4248054 -0.3721979
## Leptolyngbya sp. PCC 7376                         -0.4248054  0.1929915
## Pseudanabaena sp. ABRG5-3                         -0.4248054 -1.5025766
## Pseudanabaena sp. PCC 7367                        -0.4248054 -1.5025766
## Anabaena sp. 90                                   -0.4248054 -1.5025766
## Leptolyngbya D1                                   -0.4248054  0.7581809
## Phormidesmis                                      -0.4248054 -0.3721979
## Leptolyngbya D3                                   -0.4248054  0.1929915
## Anabaena                                          -0.4248054  0.1929915
##                                                        ace f       dea d
## Nostoc punctiforme PCC 73102                       4.1411822  1.14599909
## Oscillatoria acuminata PCC 6304                   -0.3865103 -1.05217460
## Oscillatoria nigro-viridis PCC 7112               -0.3865103  0.59645567
## Nostoc flagelliforme CCNUN1                       -0.3865103  2.51985765
## Trichodesmium erythraeum IMS101                   -0.3865103 -0.22785947
## Nostoc sphaeroides                                -0.3865103  0.04691224
## Leptolyngbya boryana NIES-2135                     1.8773359 -0.77740289
## Nostoc azollae' 0708                              -0.3865103 -1.60171803
## Dolichospermum flos-aquae CCAP 1403/13F           -0.3865103 -0.77740289
## Anabaena cylindrica PCC 7122                      -0.3865103  0.32168396
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.3865103  2.24508594
## Nostoc sp. TCL240-02                               1.8773359 -0.22785947
## Nostoc sp. C052                                   -0.3865103  1.14599909
## Nostoc edaphicum CCNP1411                         -0.3865103 -0.50263118
## Nostoc sp. C057                                   -0.3865103  1.14599909
## Microcoleus sp. PCC 7113                          -0.3865103  0.32168396
## Nostoc sp. ATCC 53789                              1.8773359  0.32168396
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.3865103  0.87122738
## Calothrix sp. PCC 7507                            -0.3865103  0.32168396
## Anabaena sp. YBS01                                -0.3865103  0.32168396
## Nostoc sp. CENA543                                -0.3865103 -0.22785947
## Nostoc sp. TCL26-01                               -0.3865103 -0.50263118
## Nostoc sp. PCC 7524                               -0.3865103 -0.77740289
## Nostoc sphaeroides CCNUC1                         -0.3865103  2.51985765
## Nostoc sp. NIES-3756                              -0.3865103 -0.22785947
## Nostoc sp. PCC 7120 = FACHB-418                   -0.3865103  0.32168396
## Leptolyngbya sp. NIES-3755                        -0.3865103 -0.22785947
## Leptolyngbya boryana dg5                           1.8773359 -0.77740289
## Leptolyngbya boryana IAM M-101                     1.8773359 -0.77740289
## Anabaena sp. WA102                                -0.3865103 -0.77740289
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.3865103  0.04691224
## Leptolyngbya sp. O-77                             -0.3865103 -0.22785947
## Dolichospermum sp. UHCC 0315A                     -0.3865103 -0.22785947
## Leptolyngbya sp. PCC 7376                         -0.3865103 -1.60171803
## Pseudanabaena sp. ABRG5-3                         -0.3865103  0.87122738
## Pseudanabaena sp. PCC 7367                        -0.3865103 -1.87648974
## Anabaena sp. 90                                   -0.3865103 -1.05217460
## Leptolyngbya D1                                   -0.3865103 -0.22785947
## Phormidesmis                                      -0.3865103  0.32168396
## Leptolyngbya D3                                   -0.3865103 -0.22785947
## Anabaena                                          -0.3865103 -0.50263118
##                                                        inf b tig        rnr
## Nostoc punctiforme PCC 73102                      -0.5588104   0 -0.2605158
## Oscillatoria acuminata PCC 6304                   -0.5588104   0 -0.2605158
## Oscillatoria nigro-viridis PCC 7112                1.2035916   0 -0.2605158
## Nostoc flagelliforme CCNUN1                        1.2035916   0 -0.2605158
## Trichodesmium erythraeum IMS101                   -0.5588104   0 -0.2605158
## Nostoc sphaeroides                                -0.5588104   0 -0.2605158
## Leptolyngbya boryana NIES-2135                    -0.5588104   0 -0.2605158
## Nostoc azollae' 0708                              -0.5588104   0 -0.2605158
## Dolichospermum flos-aquae CCAP 1403/13F            1.2035916   0 -0.2605158
## Anabaena cylindrica PCC 7122                      -0.5588104   0 -0.2605158
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6   1.2035916   0 -0.2605158
## Nostoc sp. TCL240-02                              -0.5588104   0 -0.2605158
## Nostoc sp. C052                                    1.2035916   0 -0.2605158
## Nostoc edaphicum CCNP1411                          1.2035916   0 -0.2605158
## Nostoc sp. C057                                    1.2035916   0 -0.2605158
## Microcoleus sp. PCC 7113                           1.2035916   0 -0.2605158
## Nostoc sp. ATCC 53789                             -0.5588104   0 -0.2605158
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.5588104   0 -0.2605158
## Calothrix sp. PCC 7507                             1.2035916   0 -0.2605158
## Anabaena sp. YBS01                                -0.5588104   0 -0.2605158
## Nostoc sp. CENA543                                -0.5588104   0 -0.2605158
## Nostoc sp. TCL26-01                               -0.5588104   0 -0.2605158
## Nostoc sp. PCC 7524                               -0.5588104   0 -0.2605158
## Nostoc sphaeroides CCNUC1                         -0.5588104   0 -0.2605158
## Nostoc sp. NIES-3756                              -0.5588104   0 -0.2605158
## Nostoc sp. PCC 7120 = FACHB-418                   -0.5588104   0 -0.2605158
## Leptolyngbya sp. NIES-3755                         2.9659935   0 -0.2605158
## Leptolyngbya boryana dg5                          -0.5588104   0 -0.2605158
## Leptolyngbya boryana IAM M-101                    -0.5588104   0 -0.2605158
## Anabaena sp. WA102                                -0.5588104   0 -0.2605158
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.5588104   0 -0.2605158
## Leptolyngbya sp. O-77                             -0.5588104   0  2.4097716
## Dolichospermum sp. UHCC 0315A                     -0.5588104   0 -0.2605158
## Leptolyngbya sp. PCC 7376                          2.9659935   0  2.4097716
## Pseudanabaena sp. ABRG5-3                         -0.5588104   0 -0.2605158
## Pseudanabaena sp. PCC 7367                        -0.5588104   0 -0.2605158
## Anabaena sp. 90                                   -0.5588104   0 -0.2605158
## Leptolyngbya D1                                   -0.5588104   0 -0.2605158
## Phormidesmis                                      -0.5588104   0  5.0800590
## Leptolyngbya D3                                   -0.5588104   0 -0.2605158
## Anabaena                                          -0.5588104   0 -0.2605158
##                                                         dna a       hup b
## Nostoc punctiforme PCC 73102                      -0.58034566  1.48022191
## Oscillatoria acuminata PCC 6304                   -0.03956902 -0.88428842
## Oscillatoria nigro-viridis PCC 7112               -0.03956902 -0.88428842
## Nostoc flagelliforme CCNUN1                       -0.58034566  0.69205180
## Trichodesmium erythraeum IMS101                    0.50120762 -0.88428842
## Nostoc sphaeroides                                -0.58034566 -0.88428842
## Leptolyngbya boryana NIES-2135                     2.12353753  0.69205180
## Nostoc azollae' 0708                              -0.58034566 -0.88428842
## Dolichospermum flos-aquae CCAP 1403/13F           -0.03956902 -0.88428842
## Anabaena cylindrica PCC 7122                      -0.58034566 -0.88428842
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.58034566  2.26839202
## Nostoc sp. TCL240-02                              -0.58034566  0.69205180
## Nostoc sp. C052                                   -0.03956902  1.48022191
## Nostoc edaphicum CCNP1411                         -0.58034566 -0.88428842
## Nostoc sp. C057                                   -0.03956902  1.48022191
## Microcoleus sp. PCC 7113                          -0.03956902 -0.88428842
## Nostoc sp. ATCC 53789                             -0.58034566  0.69205180
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.58034566 -0.88428842
## Calothrix sp. PCC 7507                             0.50120762 -0.88428842
## Anabaena sp. YBS01                                -0.03956902 -0.09611831
## Nostoc sp. CENA543                                -0.03956902 -0.88428842
## Nostoc sp. TCL26-01                               -0.58034566 -0.09611831
## Nostoc sp. PCC 7524                               -0.58034566 -0.09611831
## Nostoc sphaeroides CCNUC1                         -0.58034566 -0.88428842
## Nostoc sp. NIES-3756                              -0.03956902 -0.88428842
## Nostoc sp. PCC 7120 = FACHB-418                   -0.58034566 -0.88428842
## Leptolyngbya sp. NIES-3755                        -0.03956902 -0.09611831
## Leptolyngbya boryana dg5                           1.58276089  0.69205180
## Leptolyngbya boryana IAM M-101                     1.58276089  0.69205180
## Anabaena sp. WA102                                -0.58034566 -0.09611831
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.03956902 -0.88428842
## Leptolyngbya sp. O-77                             -0.58034566 -0.88428842
## Dolichospermum sp. UHCC 0315A                     -0.03956902 -0.09611831
## Leptolyngbya sp. PCC 7376                         -0.58034566 -0.88428842
## Pseudanabaena sp. ABRG5-3                         -0.03956902  0.69205180
## Pseudanabaena sp. PCC 7367                         4.82742072  0.69205180
## Anabaena sp. 90                                   -0.03956902 -0.09611831
## Leptolyngbya D1                                   -0.58034566  2.26839202
## Phormidesmis                                       0.50120762 -0.09611831
## Leptolyngbya D3                                   -0.58034566  2.26839202
## Anabaena                                          -0.03956902 -0.09611831
##                                                        rbf a      yfl a
## Nostoc punctiforme PCC 73102                      -0.3247635 -0.2116563
## Oscillatoria acuminata PCC 6304                   -0.3247635  2.6809803
## Oscillatoria nigro-viridis PCC 7112               -0.3247635 -0.2116563
## Nostoc flagelliforme CCNUN1                       -0.3247635 -0.2116563
## Trichodesmium erythraeum IMS101                   -0.3247635  5.5736170
## Nostoc sphaeroides                                -0.3247635 -0.2116563
## Leptolyngbya boryana NIES-2135                    -0.3247635 -0.2116563
## Nostoc azollae' 0708                              -0.3247635 -0.2116563
## Dolichospermum flos-aquae CCAP 1403/13F            3.0040623 -0.2116563
## Anabaena cylindrica PCC 7122                      -0.3247635 -0.2116563
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.3247635 -0.2116563
## Nostoc sp. TCL240-02                              -0.3247635 -0.2116563
## Nostoc sp. C052                                   -0.3247635 -0.2116563
## Nostoc edaphicum CCNP1411                         -0.3247635 -0.2116563
## Nostoc sp. C057                                   -0.3247635 -0.2116563
## Microcoleus sp. PCC 7113                          -0.3247635 -0.2116563
## Nostoc sp. ATCC 53789                             -0.3247635 -0.2116563
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.3247635 -0.2116563
## Calothrix sp. PCC 7507                            -0.3247635 -0.2116563
## Anabaena sp. YBS01                                -0.3247635 -0.2116563
## Nostoc sp. CENA543                                -0.3247635 -0.2116563
## Nostoc sp. TCL26-01                               -0.3247635 -0.2116563
## Nostoc sp. PCC 7524                               -0.3247635 -0.2116563
## Nostoc sphaeroides CCNUC1                         -0.3247635 -0.2116563
## Nostoc sp. NIES-3756                              -0.3247635 -0.2116563
## Nostoc sp. PCC 7120 = FACHB-418                   -0.3247635 -0.2116563
## Leptolyngbya sp. NIES-3755                        -0.3247635 -0.2116563
## Leptolyngbya boryana dg5                          -0.3247635 -0.2116563
## Leptolyngbya boryana IAM M-101                    -0.3247635 -0.2116563
## Anabaena sp. WA102                                 3.0040623 -0.2116563
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.3247635 -0.2116563
## Leptolyngbya sp. O-77                             -0.3247635 -0.2116563
## Dolichospermum sp. UHCC 0315A                      3.0040623 -0.2116563
## Leptolyngbya sp. PCC 7376                         -0.3247635 -0.2116563
## Pseudanabaena sp. ABRG5-3                         -0.3247635 -0.2116563
## Pseudanabaena sp. PCC 7367                        -0.3247635 -0.2116563
## Anabaena sp. 90                                   -0.3247635 -0.2116563
## Leptolyngbya D1                                   -0.3247635 -0.2116563
## Phormidesmis                                      -0.3247635 -0.2116563
## Leptolyngbya D3                                   -0.3247635 -0.2116563
## Anabaena                                           3.0040623 -0.2116563
##                                                          pnp        csp
## Nostoc punctiforme PCC 73102                       1.8411606 -0.8583234
## Oscillatoria acuminata PCC 6304                    1.8411606  0.8174508
## Oscillatoria nigro-viridis PCC 7112               -0.1990444  0.8174508
## Nostoc flagelliforme CCNUN1                       -0.1990444 -0.8583234
## Trichodesmium erythraeum IMS101                   -0.1990444 -0.8583234
## Nostoc sphaeroides                                -0.1990444 -0.8583234
## Leptolyngbya boryana NIES-2135                    -0.1990444  0.8174508
## Nostoc azollae' 0708                              -0.1990444 -0.8583234
## Dolichospermum flos-aquae CCAP 1403/13F           -0.1990444  0.8174508
## Anabaena cylindrica PCC 7122                      -0.1990444 -0.8583234
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6   3.8813657 -0.8583234
## Nostoc sp. TCL240-02                              -0.1990444 -0.8583234
## Nostoc sp. C052                                   -0.1990444 -0.8583234
## Nostoc edaphicum CCNP1411                         -0.1990444 -0.8583234
## Nostoc sp. C057                                   -0.1990444 -0.8583234
## Microcoleus sp. PCC 7113                          -0.1990444  0.8174508
## Nostoc sp. ATCC 53789                              1.8411606  0.8174508
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont' -0.1990444  0.8174508
## Calothrix sp. PCC 7507                            -0.1990444 -0.8583234
## Anabaena sp. YBS01                                -0.1990444  0.8174508
## Nostoc sp. CENA543                                -0.1990444 -0.8583234
## Nostoc sp. TCL26-01                               -0.1990444 -0.8583234
## Nostoc sp. PCC 7524                               -0.1990444  0.8174508
## Nostoc sphaeroides CCNUC1                         -0.1990444 -0.8583234
## Nostoc sp. NIES-3756                              -0.1990444 -0.8583234
## Nostoc sp. PCC 7120 = FACHB-418                   -0.1990444 -0.8583234
## Leptolyngbya sp. NIES-3755                        -0.1990444 -0.8583234
## Leptolyngbya boryana dg5                          -0.1990444  0.8174508
## Leptolyngbya boryana IAM M-101                    -0.1990444  0.8174508
## Anabaena sp. WA102                                 1.8411606  2.4932250
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.1990444 -0.8583234
## Leptolyngbya sp. O-77                             -0.1990444  0.8174508
## Dolichospermum sp. UHCC 0315A                     -0.1990444  0.8174508
## Leptolyngbya sp. PCC 7376                         -2.2392494  0.8174508
## Pseudanabaena sp. ABRG5-3                         -0.1990444 -0.8583234
## Pseudanabaena sp. PCC 7367                        -2.2392494 -0.8583234
## Anabaena sp. 90                                   -0.1990444  0.8174508
## Leptolyngbya D1                                   -0.1990444  0.8174508
## Phormidesmis                                      -0.1990444  2.4932250
## Leptolyngbya D3                                   -0.1990444 -0.8583234
## Anabaena                                          -0.1990444  0.8174508
##                                                        ace e       des a
## Nostoc punctiforme PCC 73102                      -0.4621465  1.39967730
## Oscillatoria acuminata PCC 6304                   -2.5674808  0.08547648
## Oscillatoria nigro-viridis PCC 7112               -0.4621465 -1.22872435
## Nostoc flagelliforme CCNUN1                       -0.4621465  1.83774424
## Trichodesmium erythraeum IMS101                   -0.4621465 -1.22872435
## Nostoc sphaeroides                                -0.4621465 -0.79065741
## Leptolyngbya boryana NIES-2135                    -0.4621465  0.96161036
## Nostoc azollae' 0708                               1.6431877 -1.66679129
## Dolichospermum flos-aquae CCAP 1403/13F           -0.4621465 -0.35259046
## Anabaena cylindrica PCC 7122                       1.6431877 -0.79065741
## Nostoc sp. 'Peltigera membranacea cyanobiont' N6  -0.4621465 -0.79065741
## Nostoc sp. TCL240-02                               1.6431877  1.83774424
## Nostoc sp. C052                                   -0.4621465  1.83774424
## Nostoc edaphicum CCNP1411                         -0.4621465  0.96161036
## Nostoc sp. C057                                   -0.4621465  1.83774424
## Microcoleus sp. PCC 7113                          -0.4621465 -0.79065741
## Nostoc sp. ATCC 53789                             -0.4621465  0.08547648
## Nostoc sp. 'Lobaria pulmonaria (5183) cyanobiont'  1.6431877  0.08547648
## Calothrix sp. PCC 7507                             1.6431877 -0.79065741
## Anabaena sp. YBS01                                 1.6431877 -0.35259046
## Nostoc sp. CENA543                                 1.6431877 -0.35259046
## Nostoc sp. TCL26-01                                1.6431877 -0.79065741
## Nostoc sp. PCC 7524                                1.6431877 -1.22872435
## Nostoc sphaeroides CCNUC1                         -0.4621465  0.08547648
## Nostoc sp. NIES-3756                              -0.4621465  0.08547648
## Nostoc sp. PCC 7120 = FACHB-418                    1.6431877 -0.35259046
## Leptolyngbya sp. NIES-3755                        -0.4621465  0.52354342
## Leptolyngbya boryana dg5                          -0.4621465  0.96161036
## Leptolyngbya boryana IAM M-101                    -0.4621465  0.96161036
## Anabaena sp. WA102                                -0.4621465  0.52354342
## Thermoleptolyngbya sp. PKUAC-SCTA183              -0.4621465 -0.79065741
## Leptolyngbya sp. O-77                             -0.4621465 -0.79065741
## Dolichospermum sp. UHCC 0315A                     -0.4621465 -0.79065741
## Leptolyngbya sp. PCC 7376                         -0.4621465  0.08547648
## Pseudanabaena sp. ABRG5-3                         -0.4621465  0.52354342
## Pseudanabaena sp. PCC 7367                        -0.4621465 -1.66679129
## Anabaena sp. 90                                   -0.4621465 -0.79065741
## Leptolyngbya D1                                   -0.4621465  0.52354342
## Phormidesmis                                      -0.4621465  0.52354342
## Leptolyngbya D3                                   -0.4621465  1.39967730
## Anabaena                                          -0.4621465 -0.79065741
```






```r
png(file="fil_matrix_heatmap.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*275,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size)

row_distance = dist(mat_data, method = "manhattan") #I dont understand distance so I'm not going to use it yet.
row_cluster = hclust(row_distance, method = "ward.D")
col_distance = dist(t(mat_data), method = "manhattan")
col_cluster = hclust(col_distance, method = "ward.D")

Colors=brewer.pal(11,"PuBuGn")  
```

```
## Warning in brewer.pal(11, "PuBuGn"): n too large, allowed maximum for palette PuBuGn is 9
## Returning the palette you asked for with that many colors
```

```r
#Colors=c("white","blue", "green")
#Colors=colorRampPalette(Colors)(100)
heatmap.2(mat_data, 
          margins = c(7, 19), # words aren't cut off in png output
          trace = "none", # trace is the cyan histogram, get rid of it
          density.info="density", # density plot, not histogram
          col=Colors,
          Rowv = as.dendrogram(row_cluster), # changing the cluster and distance methods from the default
          Colv = as.dendrogram(col_cluster)) 
dev.off()
```

```
## quartz_off_screen 
##                 2
```

<img src="fil_matrix_heatmap.png" width="1500" />

