# Homework for STAT540
Sam Hinshaw  



Please view the HTML version of this document, not the GitHub rendering of the `.md` file for proper Table of Contents. 

# Setup

```r
# library(BiocInstaller) # or source("https://bioconductor.org/biocLite.R")
suppressPackageStartupMessages({
	library(limma) # biocLite("limma")
	library(edgeR) # biocLite("edgeR")
	library(DESeq2) # biocLite("DESeq2")
	library(biomaRt)
	library(magrittr)
	library(devtools)
	library(broom)
	library(ggplot2)
	library(reshape2)
	library(readr)
	library(knitr)
	library(xtable)
	library(pander)
	library(tidyr)
	library(RColorBrewer)
	library(dplyr)
})
```

Let's check our R version & Packages are up-to-date
We want `R >= 3.2.3`, `limma >= 3.26.7`,  & `edgeR >= 3.12.0`

```r
sessioninfo <- devtools::session_info()
packages <- sessioninfo$packages
packages %>% 
	filter(package %in% c("limma", "edgeR"))
```

```
##   package * version       date       source
## 1   edgeR *  3.12.0 2015-12-17 Bioconductor
## 2   limma *  3.26.7 2016-01-30 Bioconductor
```

```r
sessioninfo$platform$version
```

```
## [1] "R version 3.2.3 (2015-12-10)"
```

Looks good, let's continue

# Data Inspection

## Download Data

For this, I am using a makefile to avoid downloading the datasets more than once. 
See [Makefile](./Makefile) in this directory. Contents pasted here:
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

## Experimental Design

Before we do anything, let's see what's what in our experiment, [GSE10718](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10718)
<blockquote>
Gene expression patterns were assessed in normal human bronchial epithelial (NHBE) cells exposed to cigarette smoke (CS) from a typical "full flavor" American brand of cigarettes in order to develop a better understanding of the genomic impact of tobacco exposure, which can ultimately define biomarkers that discriminate tobacco-related effects and outcomes in a clinical setting. NHBE cells were treated with CS for 15 minutes and alterations to the transcriptome assessed at 1,2,4 and 24 hours post-CS-exposure using high-density oligonucleotide microarrays.
</blockquote>

## Import Data

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

## Inspect Data

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

### Fix Import Errors

Well, since we're not supposed to edit our souce file, I've copied the file do `data_fixed.txt` and appended `"probeIDs	"` to the beginning of the first line with `sed`.  Now we should be ready to rock!

```r
data <- read_delim("data_fixed.txt", delim = "\t")
fitdata <- data
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
glimpse(design)
```

```
## Observations: 23
## Variables: 4
## $ InternalID (chr) "GSE10718_Biomat_1", "GSE10718_Biomat_10", "GSE1071...
## $ ExternalID (chr) "GSM270883", "GSM270874", "GSM270873", "GSM270872",...
## $ Treatment  (chr) "control", "control", "control", "control", "contro...
## $ time       (chr) "24_h", "1_h", "1_h", "1_h", "4_h", "4_h", "4_h", "...
```

### Inspect Data

There we go! So we've got 23 samples. How many different treatment groups, and how many different timepoints?

```r
design %>% 
	group_by(Treatment) %>% 
	tally()
```

```
## Source: local data frame [2 x 2]
## 
##         Treatment     n
##             (chr) (int)
## 1 cigarette_smoke    11
## 2         control    12
```

```r
design %>% 
	group_by(time) %>% 
	tally()
```

```
## Source: local data frame [4 x 2]
## 
##    time     n
##   (chr) (int)
## 1   1_h     5
## 2  24_h     6
## 3   2_h     6
## 4   4_h     6
```

```r
design %>% 
	group_by(Treatment, time) %>% 
	tally()
```

```
## Source: local data frame [8 x 3]
## Groups: Treatment [?]
## 
##         Treatment  time     n
##             (chr) (chr) (int)
## 1 cigarette_smoke   1_h     2
## 2 cigarette_smoke  24_h     3
## 3 cigarette_smoke   2_h     3
## 4 cigarette_smoke   4_h     3
## 5         control   1_h     3
## 6         control  24_h     3
## 7         control   2_h     3
## 8         control   4_h     3
```

That's a pretty good summary of what's going on, I'd say. Now back to our data file. How many different genes do we have reported?

```r
nrow(data)
```

```
## [1] 22737
```

22,737 probes.  Do we know any microarrays that have exactly that many genes?  Fortunately, we don't need to wonder. This study, [GSE10718](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10718) used [GPL570](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570), the well known Affymetrix HG U133+ 2.0 Array. Later we can use biomaRt to convert our probeIDs to ensembl IDs (or another equivalent), given the platform information. 

## Basic Data Manipulation
Before we continue, let's create a new column in our dataset that holds our time post-treatment as hours.


```r
design %<>% 
	mutate(hours = gsub("\\_h$","", time) %>% as.numeric())
glimpse(design)
```

```
## Observations: 23
## Variables: 5
## $ InternalID (chr) "GSE10718_Biomat_1", "GSE10718_Biomat_10", "GSE1071...
## $ ExternalID (chr) "GSM270883", "GSM270874", "GSM270873", "GSM270872",...
## $ Treatment  (chr) "control", "control", "control", "control", "contro...
## $ time       (chr) "24_h", "1_h", "1_h", "1_h", "4_h", "4_h", "4_h", "...
## $ hours      (dbl) 24, 1, 1, 1, 4, 4, 4, 1, 1, 2, 24, 2, 2, 4, 4, 4, 2...
```


## Basic Graphing

However, for now we can just plot some intensity distributions. To help our plot colors out, we can convert our Treatment and time variables to factors. This is more for my own benefit (to see intensity distributions) than it is to follow the homework to the letter. Homework progress will continue below. 

```r
design$Treatment <- as.factor(design$Treatment)
design$time <- as.factor(design$time)
```

### Intensity Distributions

Let's get biomaRt out of the way so we can move on.  

```r
ensembl <- useMart("ensembl") # set up biomaRt to use the ensembl database
datasets <- listDatasets(ensembl) # store our datasets so we can see which to use
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl) # use the human data set!
filters <- listFilters(ensembl) # find what filters we can use and cross fingers that it's here
attributes <- listAttributes(ensembl) # ditto with the attributes
affy_probe_IDs <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"), filters = "affy_hg_u133_plus_2", values = data$ProbeID, mart = ensembl) # query the probeIDs
head(affy_probe_IDs)
```

```
##   ensembl_gene_id affy_hg_u133_plus_2
## 1 ENSG00000198888        1553551_s_at
## 2 ENSG00000210100        1553551_s_at
## 3 ENSG00000210112        1553551_s_at
## 4 ENSG00000198763        1553551_s_at
## 5 ENSG00000198804          1553569_at
## 6 ENSG00000198804        1553538_s_at
```

```r
colnames(affy_probe_IDs) <- c("ensembl_gene_id", "ProbeID")
data <- inner_join(data, affy_probe_IDs, by = "ProbeID")
glimpse(data)
```

```
## Observations: 21,323
## Variables: 25
## $ ProbeID            (chr) "1294_at", "1294_at", "1316_at", "1320_at",...
## $ GSE10718_Biomat_1  (dbl) 7.900022, 7.900022, 6.126008, 6.822491, 6.6...
## $ GSE10718_Biomat_10 (dbl) 7.623666, 7.623666, 5.846865, 7.616309, 5.6...
## $ GSE10718_Biomat_11 (dbl) 8.224398, 8.224398, 6.299318, 7.217153, 6.2...
## $ GSE10718_Biomat_12 (dbl) 8.009710, 8.009710, 5.910947, 6.172960, 6.0...
## $ GSE10718_Biomat_13 (dbl) 7.607968, 7.607968, 6.623939, 7.201396, 6.8...
## $ GSE10718_Biomat_14 (dbl) 7.448574, 7.448574, 6.824437, 7.397821, 7.0...
## $ GSE10718_Biomat_15 (dbl) 7.414247, 7.414247, 6.884227, 6.782186, 6.5...
## $ GSE10718_Biomat_16 (dbl) 7.671862, 7.671862, 6.830394, 6.987745, 6.4...
## $ GSE10718_Biomat_17 (dbl) 7.892385, 7.892385, 6.453085, 5.845809, 6.1...
## $ GSE10718_Biomat_19 (dbl) 7.804841, 7.804841, 6.658476, 6.831755, 6.1...
## $ GSE10718_Biomat_2  (dbl) 7.721433, 7.721433, 6.569794, 6.674275, 6.7...
## $ GSE10718_Biomat_20 (dbl) 7.548707, 7.548707, 5.806634, 6.282316, 5.5...
## $ GSE10718_Biomat_21 (dbl) 6.343507, 6.343507, 7.064967, 7.234590, 5.7...
## $ GSE10718_Biomat_22 (dbl) 7.797584, 7.797584, 6.485319, 6.750336, 5.5...
## $ GSE10718_Biomat_23 (dbl) 7.951058, 7.951058, 7.084865, 7.111929, 7.1...
## $ GSE10718_Biomat_24 (dbl) 7.877641, 7.877641, 6.237747, 6.972026, 6.3...
## $ GSE10718_Biomat_3  (dbl) 7.510909, 7.510909, 6.384493, 6.979904, 6.7...
## $ GSE10718_Biomat_4  (dbl) 7.309565, 7.309565, 6.885804, 7.621911, 6.4...
## $ GSE10718_Biomat_5  (dbl) 8.004972, 8.004972, 7.100588, 7.294042, 5.7...
## $ GSE10718_Biomat_6  (dbl) 7.502813, 7.502813, 6.421699, 7.684150, 6.2...
## $ GSE10718_Biomat_7  (dbl) 7.500043, 7.500043, 6.269029, 7.116839, 6.6...
## $ GSE10718_Biomat_8  (dbl) 7.509491, 7.509491, 6.534320, 6.613354, 6.5...
## $ GSE10718_Biomat_9  (dbl) 7.735819, 7.735819, 6.888977, 6.708225, 7.0...
## $ ensembl_gene_id    (chr) "ENSG00000263506", "ENSG00000182179", "ENSG...
```
Now we can `gather()` our data for easy plotting, and even join our design data as well. 

```r
gathered_data <- data %>% 
	gather("InternalID", "intensity", 2:24)
glimpse(gathered_data)
```

```
## Observations: 490,429
## Variables: 4
## $ ProbeID         (chr) "1294_at", "1294_at", "1316_at", "1320_at", "1...
## $ ensembl_gene_id (chr) "ENSG00000263506", "ENSG00000182179", "ENSG000...
## $ InternalID      (chr) "GSE10718_Biomat_1", "GSE10718_Biomat_1", "GSE...
## $ intensity       (dbl) 7.900022, 7.900022, 6.126008, 6.822491, 6.6335...
```

```r
gathered_data <- inner_join(gathered_data, design, by = "InternalID")
glimpse(gathered_data)
```

```
## Observations: 490,429
## Variables: 8
## $ ProbeID         (chr) "1294_at", "1294_at", "1316_at", "1320_at", "1...
## $ ensembl_gene_id (chr) "ENSG00000263506", "ENSG00000182179", "ENSG000...
## $ InternalID      (chr) "GSE10718_Biomat_1", "GSE10718_Biomat_1", "GSE...
## $ intensity       (dbl) 7.900022, 7.900022, 6.126008, 6.822491, 6.6335...
## $ ExternalID      (chr) "GSM270883", "GSM270883", "GSM270883", "GSM270...
## $ Treatment       (fctr) control, control, control, control, control, ...
## $ time            (fctr) 24_h, 24_h, 24_h, 24_h, 24_h, 24_h, 24_h, 24_...
## $ hours           (dbl) 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24...
```

Before we plot, let's set up our color palette.

```r
display.brewer.all()
```

![](SamHinshawHomework_files/figure-html/color palette 1-1.png)

```r
timepoints <- brewer.pal(4, "Set1")
names(timepoints) <- unique(design$time)
treatments <- c("#7FC97F", "#FDC086")
names(treatments) <- unique(design$Treatment)
```

And now, the plot.  First with time-coding...

```r
p <- ggplot(gathered_data, aes(x = intensity))
p + geom_density(aes(fill = time)) + 
	scale_fill_manual(values = timepoints) + xlab("log2 intensity") + ylab("Frequency") + 
	facet_wrap( ~ InternalID, ncol = 5) + ggtitle("Log2 Intensities of Samples in GSE10718")
```

![](SamHinshawHomework_files/figure-html/intensity plot faceted time-coded-1.png)

...and then with treatment coding.

```r
p <- ggplot(gathered_data, aes(x = intensity))
p + geom_density(aes(fill = Treatment)) + 
	scale_fill_manual(values = treatments) + xlab("log2 intensity") + ylab("Frequency") + 
	facet_wrap( ~ InternalID, ncol = 5) + ggtitle("Log2 Intensities of Samples in GSE10718")
```

![](SamHinshawHomework_files/figure-html/intensity plot faceted Treatment-coded-1.png)

### Single Probe Intensity

As per the homework, let's plot the gene expression data for a single probe. Let's pick one at pseudo-random!

```r
set.seed(20)
data[sample(nrow(data), 1), 1:6] %>% kable("markdown")
```



|ProbeID   | GSE10718_Biomat_1| GSE10718_Biomat_10| GSE10718_Biomat_11| GSE10718_Biomat_12| GSE10718_Biomat_13|
|:---------|-----------------:|------------------:|------------------:|------------------:|------------------:|
|231836_at |          8.371558|           7.766896|           8.000101|           8.441531|           8.507032|

```r
singleprobe <- gathered_data %>% 
	filter(ProbeID == "231836_at")
```


```r
ggplot(data = singleprobe, aes(y = intensity, x = ExternalID)) + 
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle = 65, hjust = 1))
```

![](SamHinshawHomework_files/figure-html/plot single probe-1.png)

Okay, so we see we've got pretty similar intensities throughout the samples.  But let's group them anyways! First by treatment.

```r
ggplot(data = singleprobe, aes(y = intensity, x = Treatment)) + 
	geom_bar(stat = "identity", width = 0.5) +
	theme(axis.text.x = element_text(angle = 65, hjust = 1))
```

![](SamHinshawHomework_files/figure-html/plot single probe treatments-1.png)

...then by time

```r
ggplot(data = singleprobe, aes(y = intensity, x = time)) + 
	geom_bar(stat = "identity", width = 0.5) +
	theme(axis.text.x = element_text(angle = 65, hjust = 1))
```

![](SamHinshawHomework_files/figure-html/plot single probe times-1.png)
Ah-hah! That's rather interesting...
Let's do all permutations of combinations. 

```r
gathered_data %<>% mutate(perm = paste(as.character(Treatment), as.character(time), sep = "_"))
singleprobe <- gathered_data %>% 
	filter(ProbeID == "231836_at")
```

And now to plot...

```r
ggplot(data = singleprobe, aes(y = intensity, x = perm)) + 
	geom_bar(stat = "identity", width = 0.5) +
	theme(axis.text.x = element_text(angle = 65, hjust = 1))
```

![](SamHinshawHomework_files/figure-html/plot single probe perms-1.png)

Note that instead of creating new groups to chart, we could have simply color-coded the original bar chart. For example, using varying intensities of colors for time and complementary colors for the treatment. 

```r
smoke.colors <- brewer.pal(4, "Reds")
control.colors <- brewer.pal(4, "Blues")
names(smoke.colors) <- c("cigarette_smoke_1_h", "cigarette_smoke_2_h", "cigarette_smoke_4_h", "cigarette_smoke_24_h")
names(control.colors) <- c("control_1_h", "control_2_h", "control_4_h", "control_24_h")
perm_palette <- c(control.colors, smoke.colors)
```


```r
singleprobe %<>% 
    mutate(IDnum = gsub("^GSM","", ExternalID) %>% as.numeric())
singleprobe$perm <- as.factor(singleprobe$perm)
singleprobe %<>% 
    mutate(perm = reorder(perm, IDnum, max))
```


```r
ggplot(data = singleprobe, aes(y = intensity, x = ExternalID)) + 
	geom_bar(stat = "identity", aes(fill = perm)) +
	scale_fill_manual(values = perm_palette) + labs(fill = "Conditions") +
	theme(axis.text.x = element_text(angle = 65, hjust = 1))
```

![](SamHinshawHomework_files/figure-html/plot color coded-1.png)


# Assessing Data Quality

Moving on, let's start looking at our data quality, not just visualizing it. 
## Heatmaps


```r
datacor <- data %>% 
	dplyr::select(-ensembl_gene_id, -ProbeID) %>% 
	cor()
melted_corr <- melt(datacor) %>% tbl_df()
head(melted_corr)
```

```
## Source: local data frame [6 x 3]
## 
##                 Var1              Var2     value
##               (fctr)            (fctr)     (dbl)
## 1  GSE10718_Biomat_1 GSE10718_Biomat_1 1.0000000
## 2 GSE10718_Biomat_10 GSE10718_Biomat_1 0.9235724
## 3 GSE10718_Biomat_11 GSE10718_Biomat_1 0.9330836
## 4 GSE10718_Biomat_12 GSE10718_Biomat_1 0.9315279
## 5 GSE10718_Biomat_13 GSE10718_Biomat_1 0.9376926
## 6 GSE10718_Biomat_14 GSE10718_Biomat_1 0.9375481
```

Let's try reordering these factors with `arrange()` to reveal discrepancies. 

```r
design$qualitative <- NULL
design %<>% 
	mutate(qualitative = paste0(Treatment, "_", as.character(hours), "_", ExternalID))
colnames(melted_corr) <- c("InternalID", "Var2", "value")
design_for_corr <- design %>% 
	select(-ExternalID)
joined_corr <- melted_corr %>% 
	inner_join(design_for_corr, by = "InternalID")
```

```
## Warning in inner_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
##
design_for_corr2 <- design_for_corr %>% 
	select(InternalID, qualitative, hours, Treatment)
colnames(design_for_corr2) <- c("Var2", "qualitative2", "hours2", "Treatment2")
joined_corr %<>% 
	inner_join(design_for_corr2, by = "Var2")
```

```
## Warning in inner_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
## Rearrange & reorder by time & then treatment
by_time <- joined_corr %>% 
	arrange(-hours, Treatment)
by_time$factororder <- 1:nrow(by_time)
by_time %<>% 
	mutate(qualitative = reorder(qualitative, factororder, max))
by_time %<>% 
	arrange(hours2, Treatment2)
by_time$factororder2 <- 1:nrow(by_time)
by_time %<>% 
	mutate(qualitative2 = reorder(qualitative2, factororder2, max))
## Rearrange & reorder by treatment & then time
by_treatment <- joined_corr %>% 
	arrange(Treatment, -hours)
by_treatment$factororder <- 1:nrow(by_treatment)
by_treatment %<>% 
	mutate(qualitative = reorder(qualitative, factororder, max))
by_treatment %<>% 
	arrange(Treatment2, hours2)
by_treatment$factororder2 <- 1:nrow(by_time)
by_treatment %<>% 
	mutate(qualitative2 = reorder(qualitative2, factororder2, max))
```


```r
tilegradient <- brewer.pal(11, "Spectral")
```

Now to plot!

```r
ggplot(joined_corr, aes(x = qualitative2, y = qualitative, fill = value)) + 
	geom_tile() +
	theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
	scale_fill_gradientn(colors = tilegradient)
```

![](SamHinshawHomework_files/figure-html/unnamed-chunk-4-1.png)

```r
ggplot(by_time, aes(x = qualitative2, y = qualitative, fill = value)) + 
	geom_tile() +
	theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
	scale_fill_gradientn(colors = tilegradient)
```

![](SamHinshawHomework_files/figure-html/unnamed-chunk-4-2.png)

```r
ggplot(by_treatment, aes(x = qualitative2, y = qualitative, fill = value)) + 
	geom_tile() +
	theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
	scale_fill_gradientn(colors = tilegradient)
```

![](SamHinshawHomework_files/figure-html/unnamed-chunk-4-3.png)

## Outliers

Okay. That's an absolutely horrendous amount of code just to reorder a damn heatmap. I'm never doing that again, I'll stick with my unordered heatmap, which told me what I needed to know. GSM270874 or "GSE10718-Biomat-10" is an outlier. No significant batch effects noticeable. 


```r
outliers <- joined_corr %>% 
	filter(InternalID == "GSE10718_Biomat_10")
outliers %<>% 
	mutate(qualitative2 = reorder(qualitative2, value, max))
## Mean Correlation for this sample
outliers %>% 
	filter(Var2 != "GSE10718_Biomat_10") %>% 
	summarize(MeanCorr = mean(value))
```

```
## Source: local data frame [1 x 1]
## 
##    MeanCorr
##       (dbl)
## 1 0.9215407
```

Okay, looks like an outlier, but how can we be sure? What's the average correlation?

```r
joined_corr %>% 
	filter(value != 1.0) %>% 
	group_by(InternalID) %>% 
	summarize(MeanCorr = mean(value)) %>% 
	arrange(-MeanCorr) %>% 
	kable("markdown")
```



|InternalID         |  MeanCorr|
|:------------------|---------:|
|GSE10718_Biomat_15 | 0.9367632|
|GSE10718_Biomat_16 | 0.9365047|
|GSE10718_Biomat_17 | 0.9348710|
|GSE10718_Biomat_13 | 0.9346969|
|GSE10718_Biomat_14 | 0.9345540|
|GSE10718_Biomat_6  | 0.9331644|
|GSE10718_Biomat_19 | 0.9331202|
|GSE10718_Biomat_3  | 0.9328543|
|GSE10718_Biomat_1  | 0.9320065|
|GSE10718_Biomat_11 | 0.9312138|
|GSE10718_Biomat_8  | 0.9309811|
|GSE10718_Biomat_2  | 0.9302614|
|GSE10718_Biomat_7  | 0.9300924|
|GSE10718_Biomat_20 | 0.9295662|
|GSE10718_Biomat_12 | 0.9291928|
|GSE10718_Biomat_22 | 0.9283850|
|GSE10718_Biomat_9  | 0.9275301|
|GSE10718_Biomat_21 | 0.9258038|
|GSE10718_Biomat_24 | 0.9254918|
|GSE10718_Biomat_5  | 0.9253976|
|GSE10718_Biomat_4  | 0.9241465|
|GSE10718_Biomat_23 | 0.9239977|
|GSE10718_Biomat_10 | 0.9215407|

```r
## What about within groups?

joined_corr %>% 
	filter(value != 1.0) %>% 
	group_by(Treatment, hours) %>% 
	summarize(MeanCorr = mean(value)) %>% 
	arrange(-MeanCorr) %>% 
	kable("markdown")
```



|Treatment       | hours|  MeanCorr|
|:---------------|-----:|---------:|
|cigarette_smoke |     1| 0.9356879|
|control         |     4| 0.9353380|
|control         |    24| 0.9317074|
|control         |     2| 0.9295345|
|cigarette_smoke |     2| 0.9294967|
|cigarette_smoke |    24| 0.9275695|
|control         |     1| 0.9273158|
|cigarette_smoke |     4| 0.9259581|

Lastly, how does this outlier sample compare to specific treatments (This is a 1hr control sample). For this I'll filter out correlation with itself.

```r
outliers %>% 
	filter(Var2 != "GSE10718_Biomat_10") %>% 
	group_by(Treatment2, hours2) %>% 
	summarize(MeanCorr = mean(value)) %>% 
	arrange(-MeanCorr) %>% 
	kable("markdown")
```



|Treatment2      | hours2|  MeanCorr|
|:---------------|------:|---------:|
|control         |      1| 0.9274908|
|cigarette_smoke |      1| 0.9264554|
|control         |     24| 0.9255085|
|control         |      4| 0.9252802|
|control         |      2| 0.9236324|
|cigarette_smoke |      2| 0.9197050|
|cigarette_smoke |      4| 0.9148929|
|cigarette_smoke |     24| 0.9129820|

At last, something meaningful!! This sample correlates MOST those samples within its own group, but also highly with the 1hr cigarette smoke group.  

# Differential Expression with Respect to Treatment

## Linear Model

```r
combn <- factor(paste(design$Treatment,
					  design$hours, sep = "_"))
```


```r
design.matrix <- model.matrix(~combn)
row.names(fitdata) <- fitdata$ProbeID
fitdata$ProbeID <- NULL
fit <- lmFit(fitdata, design.matrix)
efit <- eBayes(fit)
df <- topTable(efit, coef=2)
df %>% kable("markdown")
```



|            |    logFC|   AveExpr|        t| P.Value| adj.P.Val|        B|
|:-----------|--------:|---------:|--------:|-------:|---------:|--------:|
|204475_at   | 2.896548| 12.766287| 21.83906|       0|     0e+00| 19.00288|
|203665_at   | 4.308239| 10.745720| 21.27219|       0|     0e+00| 18.72882|
|241418_at   | 2.942141|  9.167353| 19.73735|       0|     0e+00| 17.92296|
|226650_at   | 3.495259| 11.589417| 17.97868|       0|     0e+00| 16.86965|
|203810_at   | 2.870408| 10.684728| 17.43976|       0|     0e+00| 16.51543|
|203811_s_at | 3.391046| 10.276328| 16.02841|       0|     1e-07| 15.50864|
|225061_at   | 2.636445| 10.503353| 15.83022|       0|     1e-07| 15.35740|
|214696_at   | 1.864804| 11.730937| 14.83136|       0|     2e-07| 14.55469|
|226541_at   | 2.197414| 11.269388| 14.41075|       0|     3e-07| 14.19530|
|225252_at   | 1.603954| 12.195739| 14.20703|       0|     3e-07| 14.01641|

That's good for now. 
********
This page was last updated on  Tuesday, March 08, 2016 at 05:23PM
