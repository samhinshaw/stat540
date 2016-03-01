# Sam Hinshaw Homework for STAT540



Please view the HTML version of this document, not the GitHub rendering of the `.md` file for proper Table of Contents. 

## Setup

```r
# library(BiocInstaller) # or source("https://bioconductor.org/biocLite.R")
suppressPackageStartupMessages({
	library(limma) # biocLite("limma")
	library(edgeR) # biocLite("edgeR")
	library(DESeq2) # biocLite("DESeq2")
	library(dplyr)
	library(magrittr)
	library(devtools)
	library(broom)
	library(readr)
})
```

Let's check our R version & Packages are up-to-date
We want `R >= 3.2.3`, `limma >= 3.26.7`,  & `edgeR >= 3.12.0`

```r
devtools::session_info()
```

```
## Session info --------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.2.3 (2015-12-10)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language en_CA:en                    
##  collate  en_CA.UTF-8                 
##  tz       <NA>                        
##  date     2016-03-01
```

```
## Packages ------------------------------------------------------------------
```

```
##  package              * version     date       source        
##  acepack                1.3-3.3     2014-11-24 CRAN (R 3.2.2)
##  annotate               1.48.0      2015-12-17 Bioconductor  
##  AnnotationDbi          1.32.3      2016-01-06 Bioconductor  
##  assertthat             0.1         2013-12-06 CRAN (R 3.2.2)
##  Biobase              * 2.30.0      2016-01-19 Bioconductor  
##  BiocGenerics         * 0.16.1      2015-12-17 Bioconductor  
##  BiocParallel           1.4.3       2015-12-17 Bioconductor  
##  broom                * 0.4.0       2015-11-30 CRAN (R 3.2.3)
##  cluster                2.0.3       2015-07-21 CRAN (R 3.2.2)
##  colorspace             1.2-6       2015-03-11 CRAN (R 3.2.2)
##  DBI                    0.3.1       2014-09-24 CRAN (R 3.2.2)
##  DESeq2               * 1.10.1      2016-01-06 Bioconductor  
##  devtools             * 1.10.0      2016-01-23 CRAN (R 3.2.3)
##  digest                 0.6.9       2016-01-08 CRAN (R 3.2.3)
##  dplyr                * 0.4.3       2015-09-01 CRAN (R 3.2.2)
##  edgeR                * 3.12.0      2015-12-17 Bioconductor  
##  evaluate               0.8         2015-09-18 CRAN (R 3.2.2)
##  foreign                0.8-66      2015-08-19 CRAN (R 3.2.2)
##  formatR                1.2.1       2015-09-18 CRAN (R 3.2.2)
##  Formula                1.2-1       2015-04-07 CRAN (R 3.2.2)
##  futile.logger          1.4.1       2015-04-20 CRAN (R 3.2.3)
##  futile.options         1.0.0       2010-04-06 CRAN (R 3.2.3)
##  genefilter             1.52.1      2016-01-30 Bioconductor  
##  geneplotter            1.48.0      2015-12-17 Bioconductor  
##  GenomeInfoDb         * 1.6.3       2016-01-30 Bioconductor  
##  GenomicRanges        * 1.22.4      2016-02-02 Bioconductor  
##  ggplot2                2.0.0       2015-12-18 CRAN (R 3.2.3)
##  gridExtra              2.0.0       2015-07-14 CRAN (R 3.2.2)
##  gtable                 0.1.2       2012-12-05 CRAN (R 3.2.2)
##  Hmisc                  3.17-2      2016-02-21 CRAN (R 3.2.3)
##  htmltools              0.3         2015-12-29 CRAN (R 3.2.3)
##  IRanges              * 2.4.6       2015-12-17 Bioconductor  
##  knitr                  1.12.3      2016-01-22 CRAN (R 3.2.3)
##  lambda.r               1.1.7       2015-03-20 CRAN (R 3.2.3)
##  lattice                0.20-33     2015-07-14 CRAN (R 3.2.1)
##  latticeExtra           0.6-28      2016-02-09 CRAN (R 3.2.3)
##  limma                * 3.26.7      2016-01-30 Bioconductor  
##  locfit                 1.5-9.1     2013-04-20 CRAN (R 3.2.3)
##  magrittr             * 1.5         2014-11-22 CRAN (R 3.2.2)
##  memoise                1.0.0       2016-01-29 CRAN (R 3.2.3)
##  mnormt                 1.5-3       2015-05-25 CRAN (R 3.2.2)
##  munsell                0.4.3       2016-02-13 CRAN (R 3.2.3)
##  nlme                   3.1-124     2016-01-20 CRAN (R 3.2.3)
##  nnet                   7.3-12      2016-02-02 CRAN (R 3.2.3)
##  plyr                   1.8.3       2015-06-12 CRAN (R 3.2.2)
##  psych                  1.5.8       2015-08-30 CRAN (R 3.2.2)
##  R6                     2.1.2       2016-01-26 CRAN (R 3.2.3)
##  RColorBrewer           1.1-2       2014-12-07 CRAN (R 3.2.2)
##  Rcpp                 * 0.12.3      2016-01-10 CRAN (R 3.2.3)
##  RcppArmadillo        * 0.6.500.4.0 2016-01-27 CRAN (R 3.2.3)
##  readr                * 0.2.2       2015-10-22 CRAN (R 3.2.3)
##  reshape2               1.4.1       2014-12-06 CRAN (R 3.2.2)
##  rmarkdown              0.9.5       2016-02-22 CRAN (R 3.2.3)
##  rpart                  4.1-10      2015-06-29 CRAN (R 3.2.1)
##  RSQLite                1.0.0       2014-10-25 CRAN (R 3.2.2)
##  S4Vectors            * 0.8.11      2016-01-30 Bioconductor  
##  scales                 0.3.0       2015-08-25 CRAN (R 3.2.2)
##  stringi                1.0-1       2015-10-22 CRAN (R 3.2.3)
##  stringr                1.0.0       2015-04-30 CRAN (R 3.2.3)
##  SummarizedExperiment * 1.0.2       2016-01-06 Bioconductor  
##  survival               2.38-3      2015-07-02 CRAN (R 3.2.1)
##  tidyr                  0.4.1       2016-02-05 CRAN (R 3.2.3)
##  XML                    3.98-1.3    2015-06-30 CRAN (R 3.2.2)
##  xtable                 1.8-2       2016-02-05 CRAN (R 3.2.3)
##  XVector                0.10.0      2015-12-17 Bioconductor  
##  yaml                   2.1.13      2014-06-12 CRAN (R 3.2.2)
##  zlibbioc               1.16.0      2015-12-17 Bioconductor
```

## Data Inspection

### Download & Import Data

For this, I am using a makefile to avoid downloading the datasets more than once. 
See external file in this directory, `Makefile`. The contents are pasted below:
```
all: homework.html

clean: 
	rm -rf data.txt.gz design.txt SamHinshawHomework.md SamHinshawHomework.html
	sed -i '/data.txt.gz/d' ../.gitignore
	sed -i '/design.txt/d' ../.gitignore

data.txt: 
	curl -o data.txt.gz 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_transcriptome_data.txt.gz?raw=true'

design.txt: 
	curl -o design.txt 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_design.txt?raw=true'

gitignore: 
	grep -q -F 'data.txt.gz' ../.gitignore || echo 'data.txt.gz' >> ../.gitignore
	grep -q -F 'design.txt' ../.gitignore || echo 'design.txt' >> ../.gitignore

homework.html:  SamHinshawHomework.Rmd data.txt design.txt gitignore
	Rscript -e "rmarkdown::render('SamHinshawHomework.Rmd')"
```
Now to Import

```r
data <- read_delim("data.txt", delim = "\t")
```

```
## Warning: 22737 parsing failures.
## row col   expected     actual
##   1  -- 23 columns 24 columns
##   2  -- 23 columns 24 columns
##   3  -- 23 columns 24 columns
##   4  -- 23 columns 24 columns
##   5  -- 23 columns 24 columns
## ... ... .......... ..........
## .See problems(...) for more details.
```

Huh, thanks to `readr`, we're getting an error. Something is wrong with our column count. Let's head to the next step and figure out how to handle this error.  

### Inspect Data

```r
class(data) # First let's make sure this is expected. We should have a data.frame, specifically a tbl_df. 
```

```
## [1] "tbl_df"     "tbl"        "data.frame"
```

```r
str(data) # Huh, it looks like we may have an error with our header, as we've got probe IDs labeled as a sample. 
```

```
## Classes 'tbl_df', 'tbl' and 'data.frame':	22737 obs. of  23 variables:
##  $ GSE10718_Biomat_1 : chr  "1294_at" "1316_at" "1320_at" "1431_at" ...
##  $ GSE10718_Biomat_10: num  7.9 6.13 6.82 6.63 6.95 ...
##  $ GSE10718_Biomat_11: num  7.62 5.85 7.62 5.7 7.51 ...
##  $ GSE10718_Biomat_12: num  8.22 6.3 7.22 6.22 7.29 ...
##  $ GSE10718_Biomat_13: num  8.01 5.91 6.17 6.06 7.95 ...
##  $ GSE10718_Biomat_14: num  7.61 6.62 7.2 6.89 7.87 ...
##  $ GSE10718_Biomat_15: num  7.45 6.82 7.4 7.02 7.55 ...
##  $ GSE10718_Biomat_16: num  7.41 6.88 6.78 6.59 7.27 ...
##  $ GSE10718_Biomat_17: num  7.67 6.83 6.99 6.45 7.56 ...
##  $ GSE10718_Biomat_19: num  7.89 6.45 5.85 6.18 7.24 ...
##  $ GSE10718_Biomat_2 : num  7.8 6.66 6.83 6.19 6.96 ...
##  $ GSE10718_Biomat_20: num  7.72 6.57 6.67 6.74 7.05 ...
##  $ GSE10718_Biomat_21: num  7.55 5.81 6.28 5.57 8.04 ...
##  $ GSE10718_Biomat_22: num  6.34 7.06 7.23 5.77 7.62 ...
##  $ GSE10718_Biomat_23: num  7.8 6.49 6.75 5.59 7.32 ...
##  $ GSE10718_Biomat_24: num  7.95 7.08 7.11 7.14 7.83 ...
##  $ GSE10718_Biomat_3 : num  7.88 6.24 6.97 6.34 7.03 ...
##  $ GSE10718_Biomat_4 : num  7.51 6.38 6.98 6.73 7.06 ...
##  $ GSE10718_Biomat_5 : num  7.31 6.89 7.62 6.48 8.46 ...
##  $ GSE10718_Biomat_6 : num  8 7.1 7.29 5.72 7.73 ...
##  $ GSE10718_Biomat_7 : num  7.5 6.42 7.68 6.28 7.54 ...
##  $ GSE10718_Biomat_8 : num  7.5 6.27 7.12 6.63 8.17 ...
##  $ GSE10718_Biomat_9 : num  7.51 6.53 6.61 6.52 7.67 ...
##  - attr(*, "problems")=Classes 'tbl_df', 'tbl' and 'data.frame':	22737 obs. of  4 variables:
##   ..$ row     : int  1 2 3 4 5 6 7 8 9 10 ...
##   ..$ col     : chr  NA NA NA NA ...
##   ..$ expected: chr  "23 columns" "23 columns" "23 columns" "23 columns" ...
##   ..$ actual  : chr  "24 columns" "24 columns" "24 columns" "24 columns" ...
```

Could this be an artifact of the compression? Let's try something else real quick. 

********
This page was last updated on  Tuesday, March 01, 2016 at 02:38PM
