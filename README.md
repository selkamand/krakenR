
<!-- README.md is generated from README.Rmd. Please edit that file -->

# krakenR <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of krakenR is to simplify cohort-analysis of sample analysed
with the kraken2 taxonomic classification system

## Installation

You can install the development version of krakenR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/krakenR")
```

## Quick Start

``` r
library(krakenR)

# Read in all kraken reports from a directory into a data.frame
kreport_dir <- system.file(package="krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")
df_kreports <- kraken_reports_parse(kreport_dir)

# Print Data
df_kreports
#>                          Filename PercentReadsCoveredByCladeLowResolution
#>                            <char>                                   <num>
#>      1:          e_coli_1.kreport                                    0.23
#>      2:          e_coli_1.kreport                                   99.77
#>      3:          e_coli_1.kreport                                   95.07
#>      4:          e_coli_1.kreport                                   95.07
#>      5:          e_coli_1.kreport                                   95.00
#>     ---                                                                  
#> 422236: l_monocytogenes_9.kreport                                    0.00
#> 422237: l_monocytogenes_9.kreport                                    0.00
#> 422238: l_monocytogenes_9.kreport                                    0.00
#> 422239: l_monocytogenes_9.kreport                                    0.00
#> 422240: l_monocytogenes_9.kreport                                    0.00
#>         ReadsCoveredByClade ReadsDirectlyAssigned   Rank TaxonomyID
#>                       <int>                 <int> <char>      <int>
#>      1:                   7                     7      U          0
#>      2:                2993                   141      R          1
#>      3:                2852                     0     R1     131567
#>      4:                2852                     2      D          2
#>      5:                2850                     1      P       1224
#>     ---                                                            
#> 422236:                   0                     0      S     758918
#> 422237:                   0                     0      S    1423421
#> 422238:                   0                     0      S    1737345
#> 422239:                   0                     0      S    1737346
#> 422240:                   0                     0     R1      28384
#>                                              ScientificName Level RankSimple
#>                                                      <char> <int>     <char>
#>      1:                                        unclassified     0          U
#>      2:                                                root     0          R
#>      3:                                  cellular organisms     1          R
#>      4:                                            Bacteria     2          D
#>      5:                                      Proteobacteria     3          P
#>     ---                                                                     
#> 422236:                            Rubber viroid India/2009     3          S
#> 422237: Cherry leaf scorch small circular viroid-like RNA 1     3          S
#> 422238:                             Grapevine latent viroid     3          S
#> 422239:           Apple hammerhead viroid-like circular RNA     3          S
#> 422240:                                     other sequences     1          R
#>                  SampleID TotalReadsInSample        RPM
#>                    <char>              <int>      <num>
#>      1:          e_coli_1               3000   2333.333
#>      2:          e_coli_1               3000 997666.667
#>      3:          e_coli_1               3000 950666.667
#>      4:          e_coli_1               3000 950666.667
#>      5:          e_coli_1               3000 950000.000
#>     ---                                                
#> 422236: l_monocytogenes_9               1000      0.000
#> 422237: l_monocytogenes_9               1000      0.000
#> 422238: l_monocytogenes_9               1000      0.000
#> 422239: l_monocytogenes_9               1000      0.000
#> 422240: l_monocytogenes_9               1000      0.000
```

## Visualising Results

``` r
library(plotly)
#> Loading required package: ggplot2
#> 
#> Attaching package: 'plotly'
#> The following object is masked from 'package:ggplot2':
#> 
#>     last_plot
#> The following object is masked from 'package:stats':
#> 
#>     filter
#> The following object is masked from 'package:graphics':
#> 
#>     layout

sunburst <- kraken_visualise_sunburst(df_kreports, sample = "e_coli_1", ranks = c("F", "G", "S"), ancestor = "E. coli\nSample", insidetextorientation = "tangential")

print(sunburst)
```

### Custom Visualisations

``` r
# If building more custom visualisations you may need to annotate the kraken datasets with their 'parent' species. You can do this using the kraken_annotate_parents commands
df_kreports_annotated <- df_kreports |>
  filter(Rank %in% c("F", "G", "S")) |>
  kraken_annotate_parents(ancestor = "Ancestor")
```

## Storing as a database

``` r
# Path to folder containing kraken reports 
# (reports must be named with sample_id as prefix)
kreport_dir <- system.file(package="krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")

# Create kraken sqlite database
kraken_db <- kraken_reports_parse_to_sqlite_db(kreport_dir)
```
