# Homework for STAT540
Sam Hinshaw  



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
	library(ggplot2)
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
##  ggplot2              * 2.0.0       2015-12-18 CRAN (R 3.2.3)
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
all: homework

clean: 
	rm -rf SamHinshawHomework.md SamHinshawHomework.html data_fixed.txt design_fixed.txt

data: 
	curl -o data.txt.gz 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_transcriptome_data.txt.gz?raw=true'
	gunzip -kq data.txt.gz

design: 
	curl -o design.txt 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_design.txt?raw=true'

gitignore: data.txt.gz design.txt data.txt
	grep -q -F 'data.txt' ../.gitignore || echo 'data.txt' >> ../.gitignore
	grep -q -F 'data.txt.gz' ../.gitignore || echo 'data.txt.gz' >> ../.gitignore
	grep -q -F 'design.txt' ../.gitignore || echo 'design.txt' >> ../.gitignore
	grep -q -F 'design_fixed.txt' ../.gitignore || echo 'design_fixed.txt' >> ../.gitignore
	grep -q -F 'data_fixed.txt' ../.gitignore || echo 'data_fixed.txt' >> ../.gitignore

data_fixed: data.txt
	cp data.txt data_fixed.txt
	sed -i -e 's/^"GSE10718_Biomat_1"/"ProbeID"	"GSE10718_Biomat_1"/' data_fixed.txt
	
design_fixed: design.txt
	cp design.txt design_fixed.txt
	sed -i -e 's/^"ExternalID"/"InternalID"	"ExternalID"/' design_fixed.txt
	
homework:  SamHinshawHomework.Rmd gitignore data_fixed design_fixed
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
glimpse(data)
```

```
## Observations: 22,737
## Variables: 23
## $ GSE10718_Biomat_1  (chr) "1294_at", "1316_at", "1320_at", "1431_at",...
## $ GSE10718_Biomat_10 (dbl) 7.900022, 6.126008, 6.822491, 6.633596, 6.9...
## $ GSE10718_Biomat_11 (dbl) 7.623666, 5.846865, 7.616309, 5.698143, 7.5...
## $ GSE10718_Biomat_12 (dbl) 8.224398, 6.299318, 7.217153, 6.221293, 7.2...
## $ GSE10718_Biomat_13 (dbl) 8.009710, 5.910947, 6.172960, 6.055411, 7.9...
## $ GSE10718_Biomat_14 (dbl) 7.607968, 6.623939, 7.201396, 6.885087, 7.8...
## $ GSE10718_Biomat_15 (dbl) 7.448574, 6.824437, 7.397821, 7.023587, 7.5...
## $ GSE10718_Biomat_16 (dbl) 7.414247, 6.884227, 6.782186, 6.591812, 7.2...
## $ GSE10718_Biomat_17 (dbl) 7.671862, 6.830394, 6.987745, 6.454693, 7.5...
## $ GSE10718_Biomat_19 (dbl) 7.892385, 6.453085, 5.845809, 6.176024, 7.2...
## $ GSE10718_Biomat_2  (dbl) 7.804841, 6.658476, 6.831755, 6.192702, 6.9...
## $ GSE10718_Biomat_20 (dbl) 7.721433, 6.569794, 6.674275, 6.741030, 7.0...
## $ GSE10718_Biomat_21 (dbl) 7.548707, 5.806634, 6.282316, 5.568127, 8.0...
## $ GSE10718_Biomat_22 (dbl) 6.343507, 7.064967, 7.234590, 5.768635, 7.6...
## $ GSE10718_Biomat_23 (dbl) 7.797584, 6.485319, 6.750336, 5.585600, 7.3...
## $ GSE10718_Biomat_24 (dbl) 7.951058, 7.084865, 7.111929, 7.135406, 7.8...
## $ GSE10718_Biomat_3  (dbl) 7.877641, 6.237747, 6.972026, 6.343507, 7.0...
## $ GSE10718_Biomat_4  (dbl) 7.510909, 6.384493, 6.979904, 6.727300, 7.0...
## $ GSE10718_Biomat_5  (dbl) 7.309565, 6.885804, 7.621911, 6.481674, 8.4...
## $ GSE10718_Biomat_6  (dbl) 8.004972, 7.100588, 7.294042, 5.722269, 7.7...
## $ GSE10718_Biomat_7  (dbl) 7.502813, 6.421699, 7.684150, 6.281308, 7.5...
## $ GSE10718_Biomat_8  (dbl) 7.500043, 6.269029, 7.116839, 6.625745, 8.1...
## $ GSE10718_Biomat_9  (dbl) 7.509491, 6.534320, 6.613354, 6.523722, 7.6...
```

It looks like we may have an error with our header, as we've got probe IDs labeled as a sample. Could this be an artifact of the compression? Let's try something else real quick. Also, now we know what dataset we're using! GSE10718.

```r
data_method2 <- gzfile("data.txt.gz") %>% read_delim(delim = "\t")
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

Nope, looks like it's just a problem with the source file. 

```r
rm(data_method2)
```

Well, since we're not supposed to edit our souce file, I've copied the file do `data_fixed.txt` and appended `"probeIDs	"` to the beginning of the first line with `sed`.  Now we should be ready to rock!

```r
data <- read_delim("data_fixed.txt", delim = "\t")
```

Sweet! No errors. Let's take a closer look.

```r
glimpse(data)
```

```
## Observations: 22,737
## Variables: 24
## $ ProbeID            (chr) "1294_at", "1316_at", "1320_at", "1431_at",...
## $ GSE10718_Biomat_1  (dbl) 7.900022, 6.126008, 6.822491, 6.633596, 6.9...
## $ GSE10718_Biomat_10 (dbl) 7.623666, 5.846865, 7.616309, 5.698143, 7.5...
## $ GSE10718_Biomat_11 (dbl) 8.224398, 6.299318, 7.217153, 6.221293, 7.2...
## $ GSE10718_Biomat_12 (dbl) 8.009710, 5.910947, 6.172960, 6.055411, 7.9...
## $ GSE10718_Biomat_13 (dbl) 7.607968, 6.623939, 7.201396, 6.885087, 7.8...
## $ GSE10718_Biomat_14 (dbl) 7.448574, 6.824437, 7.397821, 7.023587, 7.5...
## $ GSE10718_Biomat_15 (dbl) 7.414247, 6.884227, 6.782186, 6.591812, 7.2...
## $ GSE10718_Biomat_16 (dbl) 7.671862, 6.830394, 6.987745, 6.454693, 7.5...
## $ GSE10718_Biomat_17 (dbl) 7.892385, 6.453085, 5.845809, 6.176024, 7.2...
## $ GSE10718_Biomat_19 (dbl) 7.804841, 6.658476, 6.831755, 6.192702, 6.9...
## $ GSE10718_Biomat_2  (dbl) 7.721433, 6.569794, 6.674275, 6.741030, 7.0...
## $ GSE10718_Biomat_20 (dbl) 7.548707, 5.806634, 6.282316, 5.568127, 8.0...
## $ GSE10718_Biomat_21 (dbl) 6.343507, 7.064967, 7.234590, 5.768635, 7.6...
## $ GSE10718_Biomat_22 (dbl) 7.797584, 6.485319, 6.750336, 5.585600, 7.3...
## $ GSE10718_Biomat_23 (dbl) 7.951058, 7.084865, 7.111929, 7.135406, 7.8...
## $ GSE10718_Biomat_24 (dbl) 7.877641, 6.237747, 6.972026, 6.343507, 7.0...
## $ GSE10718_Biomat_3  (dbl) 7.510909, 6.384493, 6.979904, 6.727300, 7.0...
## $ GSE10718_Biomat_4  (dbl) 7.309565, 6.885804, 7.621911, 6.481674, 8.4...
## $ GSE10718_Biomat_5  (dbl) 8.004972, 7.100588, 7.294042, 5.722269, 7.7...
## $ GSE10718_Biomat_6  (dbl) 7.502813, 6.421699, 7.684150, 6.281308, 7.5...
## $ GSE10718_Biomat_7  (dbl) 7.500043, 6.269029, 7.116839, 6.625745, 8.1...
## $ GSE10718_Biomat_8  (dbl) 7.509491, 6.534320, 6.613354, 6.523722, 7.6...
## $ GSE10718_Biomat_9  (dbl) 7.735819, 6.888977, 6.708225, 7.063494, 7.4...
```

Looking good. These seem to be log2 transformed microarray intensity values.  How about we examine the metadata and then do a quick intensity plot. 

```r
design <- read_delim("design.txt", delim = "\t")
```

```
## Warning: 23 parsing failures.
## row col  expected    actual
##   1  -- 3 columns 4 columns
##   2  -- 3 columns 4 columns
##   3  -- 3 columns 4 columns
##   4  -- 3 columns 4 columns
##   5  -- 3 columns 4 columns
## ... ... ......... .........
## .See problems(...) for more details.
```

Ugh, same problem again? C'mon guys. Back to `sed`. 

```r
design <- read_delim("design_fixed.txt", delim = "\t")
```

********
This page was last updated on  Tuesday, March 01, 2016 at 03:41PM
