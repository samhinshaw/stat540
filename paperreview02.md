### [FloReMi: Flow Density Survival Regression Using Minimal Feature Redundancy](http://www.dx.doi.org/10.1002/cyto.a.22734)

#### Published in *Cytometry Part A*, August 2015

*****

#### Goals & Findings
This paper was published in response to the [FlowCAP IV challenge](http://flowcap.flowsite.org/), initiated to better predict clinical outcome of patients from blood draw samples.  Specifically, the goal of the FlowCAP IV challenge was to predict the time until progression to AIDS for HIV patients, given peripheral blood mononuclear cells (PBMC) analyzed in different conditions by flow cytometry. The FlowCAP consortium is run in part by the BC Cancer Agency and University of British Columbia.  
This paper, the winner of the challenge, found that their method of minimal feature redundancy worked well during cross validation when combined with predictive models such as the Cox Proportional-Hazards model, the Additive Hazards model, and the Random Survival Forests model.  However, upon testing testing novel data, the authors found that even with their method, only the Random Survival Forest model did significantly better than random, mostly due to its resilience to overfitting.  

#### Data

The authors were provided with the FlowCAP IV challenge dataset, high-dimensional (multicolor) flow cytometry dataset with sixteen different markers:  
- FSC-A, FSC-H, SSC-A - three markers describing cells' size and shape  
- IFN&gamma;, TNF&alpha;, CD4, CD27, CD107-A, CD154, CD3, CCR7, IL2, CD8, CD57, CD45RO, V-Amine/CD14 - various immune markers  
The samples were split into two groups: stimulated with HIV antigen and unstimulated. 
The authors were provided first with a training data set with which to develop their algorithm/pipeline, and later provided with a test dataset on which to test the efficacy of their algorithm.  

#### Analysis

The FloReMi algorithm performs in four main steps.  

1. Preprocessing.  
The authors needed to automate a standard flow cytometry workflow because of the high dimensionality of the data, and in order for their research to be reproducible. First is quality control, inspecting uniformity of "events", or cells detected over time to remove clumps of cells or blank spots (~5.30% removed).  Next, the authors removed margin events cells that saturated the detector when fluorescing, or events below detection limits.  These such events are not reliable measurements (~2.30% removed).  Next the authors removed doublets (cell pairs) by computing the front scatter height to area ratio (FSC-A/FSC-H) to remove further unreliable readings (~4.45% removed).  Next, the authors performed two standard flow cytometry steps, first with compensation.  Compensating for crossover between excitation spectra of multiple fluorochromes is crucial in multicolor flow cytometry, particularly when measuring 13 features.  Transformation is simply a data transformation done for easier analysis.  Finally the authors gated on live T-cells (the CD14<sup>lo</sup>/CD3<sup>hi</sup> population) with flowDensity, an automated gating program developed at the BC Cancer Research Agency. CD14 is a monocyte/macrophage marker and CD3 is a T-cell marker.  

2. Feature Extraction  
Next, the authors moved to unsupervised learning with feature extraction.  This step is composed of three main parts.  First, flowDensity was used again to determine splits for 10 of the 16 features in the dataset (FSA-A, SSC-A, CD4, CD27, CD107-A, CD154, CCR7, CD8, CD57, and CD45RO).  FlowDensity determines best split based on density distribution, and splits between peaks, as this figure illustrates.  
![Flow Density](flowdensity.png)
In this step, they excluded FSC-H, CD3, and CD14, because they were already taken into account in preprocessing.  Additionally, IFN&gamma;, TNF&alpha;, and IL-2 were removed to reduce computation time.  Though this sacrifices a great deal of information, these intracellular stains often don't have clear peaks, which makes automatic gating very difficult.  Next, the flowType dynamic programming algorithm was used to determine subsets based on these splits, identifying many cell populations that would not be manually identifiable.  In total, roughly 60,000 subsets were identified.  
Finally, "features" of each subset were extracted--mean fluorescence intensity for each of the 13 immune markers (including IFN&gamma;, TNF&alpha;, and IL-2), as well as the percentage of cells.  This left the authors with 2.5 million features per patient (3^10 * 14) x 3 (three groups - stimulated, unstimulated, and difference).  
![Defined Subsets](definedsubsets.png)

3. Feature Selection. 
This step is a supervised machine learning method, and is key to the "minimal feature redundancy" highlighted in this paper's title.  The idea was to select features with high correlation to survival time, but because of the nature of survival analysis, a simple metric such as pearson correlation could not be used due to censored values (patients that had not yet progressed to AIDS).  Therefore, the authors proceeded with the Cox proportional-hazards model.  First, each feature of the 2.5 million were fed into the hazards model--this returns a p-value and concordance index.  Next, the features must be picked.  However, features highly correlated with each other have a negative impact on training and classification (the multicollinearity problem), so redundancy must be minimized.  To combat this, the authors started adding the highest correlated features first, comparing pearson correlation *between features* as they went, and used a cutoff of corr < 0.2 to rule out correlated features.  This cutoff may seem stringent, but with 2.5 million features per patient, the authors could afford to be stringent.  Also remember that these features are the MFI of specific subtypes of cells, described by combinations of markers


4. Survival Time Prediction.  Finally, with the features selected, the Cox proportional-hazards model is used to fit a logsitic regression model to our selected features, predicting survival time of patients.  


















1. Preprocessing.  An automated approach to a standard flow cytometry workflow composed of six parts.  
	+ Quality Control - Inspect uniformity of data over time (~5.30% removed)  
	+ Remove margin events - Remove min/max, and oversaturated events (~2.30% removed)  
	+ Remove doublets - Compute FSC-A/FSC-H ratio (~4.45% removed)  
	+ Compensation - Traditional flow preprocessing  
	+ Transformation - Traditional flow preprocessing  
	+ Select live T-cells - Automatic gating with flowDensity on CD14<sup>lo</sup>/CD3<sup>hi</sup>
2. Feature Extraction.  This is an unsupervised machine learning process composed of three parts.  
	+ Determine Splits  
	+ Take hi/lo intensities from automatic gating, and use them to find subsets  
	+ Extract features of each subset as defined by flowType  
3. Feature Selection. This step is a supervised machine learning method, and is key to the "minimal feature redundancy".  
4. Survival Time Prediction.  Finally, with the features selected, the Cox proportional-hazards model is used to fit a logsitic regression model to our selected features, predicting survival time of patients.  

2. Feature Extraction
- Determine splits 
	+ flowDensity for automatic gating in one dimension.  
	+ Wonderfully optimized; no clustering!  
- Define subsets from thresholds determined by flowDensity:  
![Defined Subsets](./definedsubsets.png)
- Compute all 14 features for each subtype of each sample for both stimulated, unstimulated, and diff between stim & unstim. 
![Equation](./equation.png)

<center> **~~ Interjection ~~**  </center>

What is the Cox proportional-hazards Model?  
- Survival time is described as a probability distribution  
- Hazard Ratio: ratio between chance for event in one group vs other group  
- “proportional-hazards” means you can have multiple groups  
	+ Age, treatment, risk factors, etcetera  
- Cox PH Regression will fit & tell you what groups matter  
- “Censored values” allow for events not yet detected to be fit in regression  
- Susceptible to highly correlated values  

<center>![Cox PH Model](./coxPH.png) </center>  

3. Feature Selection
- Reduce number of features to allow for regression  
- Can’t use Pearson correlation because of censored values  
- Computer p-value of Cox proportional-hazard  
	+ Feed 2.5 million features into the hazards model  
- Sort on p-value  
- Select only uncorrelated  
	+ Pick features iteratively, discard if corr > 0.2  

4. Survival Time Prediction  
- Compute concordance  
	+ 0.5 = Random  
	+ 1.0 = Perfect  
	+ 0.0 = Predicted perfectly... just opposite  
- TNFα related features not present  
- Add uncorrelated features until concordance does NOT improve  
![Concordance Index](./coxphmodel.png)  
- Random Survival Forest  
	+ Survival trees where each split “maximizes survival distance between daughter nodes”  
	+ Takes censored values into account  
	+ Same 13 features as Cox PH  
	+ Returns MORTALITY, not survival time  
	+ Scaled mortality survival time: 0-1  
- Regularization for Semiparametric Additive Hazards Regression  
	+ Performs its own feature selection similar to LASSO or elastic net  
	+ Authors took 100 “best” features & trained with 5  

Step Done: Results  
![Table 1](./table1.png)  
- Overfitting  
	+ Random Forests inherently resilient to overfitting  
	+ Susceptible to highly correlated data, distorts randomness of trees  
	+ BUT authors attempted to control for correlation  
![Overfitting](./overfitting.png)
![Data Spread](./dataspread.png)

Findings
- Findings are limited
- Partially due to dataset
- Authors admit “mostly negative markers”
- Could be redone with new markers
![Obvious Markers](./obviousmarkers.png)
- Emphasizes possible importance of CD4+/CD8+ (double positive) T cells


#### Critique

- Exclusion of IFN-γ,IL-2, and TNFα in Feature Extraction (3 days)  
	+ IL-2 especially, T cell growth factor, and importance to T<sub>regs</sub>  
- Still used in feature selection  
![Other Papers](./otherpapers.png)  
- Cox proportional-hazards Model has drawbacks  
- Assumption that hazard ratio doesn’t change over time  
	+ Can be untrue when hazards diverge (i.e. different treatment courses)
- Even with Random Survival Forest, Cox PH Model was used for feature selection
- “Scaling” mortality to survival time seems a bit sketchy
	+ But then again, results speak for themselves
>*"Regression forests are for nonlinear multiple regression. They allow the analyst to view the importance of the predictor variables."*
>*"Survival forests are a model-free approach to survival analysis. They allow the analyst to view the importance of the covariates as the experiment evolves in time"*


*****

#### References 

