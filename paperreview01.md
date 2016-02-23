### [Phylogeographical analysis of the dominant multidrug-resistant H58 clade of Salmonella Typhi identifies inter- and intracontinental transmission events](http://www.dx.doi.org/10.1038/ng.3281)
#### Published in *Nature Genetics*, May 2015

*****

#### Overview
Typhoid, or typhoid fever is a global epidemic caused by the bacterium *Salmonella enterica* serovar Typhi, otherwise known as *S Typhi*, and the emergence of a multi-drug-resistant (MDR) strain is a cause of global concern.  The H58 lineage has been spreading over Africa and Asia over the past 30 years, replacing antibiotic-sensitive strains, and making treatment of typhoid more difficult.  
The goal of this paper was to conduct a phylogenetic analysis of the MDR typhoid lineage of the *Salmonella enterica* serovar Typhi H58, and to determine its origin and pattern of spread.  The authors of this paper were able to clearly 

#### Data
Here, data was subgenomic and whole-genome sequencing of 1,832 different S Typhi samples collected between 1992 to 2013.  A total of five samples were sequenced on PacBio's SMRT platform, allowing for better resolution of plasmid integration into the bacterial chromosome.  The raw sequence data are available in the European Nucleotide Archive (ENA), under accession identifier [ERP001718](http://www.ebi.ac.uk/ena/data/view/ERP001718).  Most of these files are outputs from Illumina HiSeq machines, with outputs downloadable as FASTQ formats.  

#### Analysis
Once all of the bacterial genomes had been sequenced, the authors aligned the Illumina reads to the [*Salmonella enterica* subsp. enterica serovar Typhi reference genome CT18](http://www.genome.jp/kegg-bin/show_organism?org=sty), and performed *de novo* assembly with their PacBio reads.  SNP analysis was performed after alignment using consensus base quality.  

The phylogenetic analysis was conducted using the SNP analysis to build a maximum-likelihood phylogenetic tree.  The authors used a "generalized time-reversible model and a Gamma distribution to model site-specific rate variation".  The authors validated this maximum likelihood by bootstrapping their model 100x.  The authors also modeled separate trees using just the H58 isolates. Finally, the authors used the geographical regions from which the samples were obtained to inform construction of the maximum-likelihood tree, and a markov model was fitted to this tree (see fig 4).  

Lastly, the authors conducted a Bayesian analysis of the evolutionary dynamics of the H58 lineage, and tested for robustness by repeating 100 times with randomly permuted tip dates.  

#### Critique
I do not have a problem with most of the methodologies employed in this paper, however there are some small items I wish the authors had taken under consideration.  Firstly, I think it would have been interesting to see the authors map their Illumina reads to a *de novo* assembled genome from their PacBio 3rd gen sequences. Overall, I wish the authors had taken greater advantage of the PacBio sequences at their disposal.  It makes sense in a phylogenetic/phylogeographics/phylotemporal analysis, where volume counts, that you would focus on your Illumina HiSeq reads. However, the authors do not indicate what extra information they gained from these long reads, which I find odd. PacBio's SMRT long reads are particularly good at analyzing mobile elements, specifically antibiotic resistant plasmids.  
