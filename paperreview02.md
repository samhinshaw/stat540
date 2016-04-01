### [FloReMi: Flow Density Survival Regression Using Minimal Feature Redundancy](http://www.dx.doi.org/10.1002/cyto.a.22734)

#### Published in *Cytometry Part A*, August 2015

*****

#### Goals & Findings
This paper was published in response to the [FlowCAP IV challenge](http://flowcap.flowsite.org/), initiated to better predict clinical outcome of patients from blood draw samples.  Specifically, the goal of the FlowCAP IV challenge was to predict the time until progression to AIDS for HIV patients, given peripheral blood mononuclear cells (PBMC) analyzed in different conditions by flow cytometry. 
This paper, the winner of the challenge, found that their method of minimal feature redundancy worked well during cross validation when combined with predictive models such as the Cox Proportional-Hazards model, the Additive Hazards model, and the Random Survival Forests model.  However, upon testing testing novel data, the authors found that even with their method, only the Random Survival Forest model did significantly better than random, mostly due to its resilience to overfitting.  

#### Data

The authors were provided with the FlowCAP IV challenge dataset, high-dimensional (multicolor) flow cytometry dataset with sixteen different markers:  
- FSC-A, FSC-H, SSC-A - three markers describing cells' size and shape  
- IFN&gamma;, TNF&alpha;, CD4, CD27, CD107-A, CD154, CD3, CCR7, IL2, CD8, CD57, CD45RO, V-Amine/CD14 - various immune markers  

#### Analysis

The FloReMi algorithm performs in four main steps.  

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
4. Survival Time Prediction.  Finally, with the features selected, the Cox Proportional Hazards model is used to fit a logsitic regression model to our selected features, predicting survival time of patients.  

2. Feature Extraction
- Determine splits 
	+ flowDensity for automatic gating in one dimension.  
	+ Wonderfully optimized; no clustering!  
- Define subsets from thresholds determined by flowDensity:  
![Defined Subsets](./definedsubsets.png)
- Compute all 14 features for each subtype of each sample for both stimulated, unstimulated, and diff between stim & unstim. 
![Equation](./equation.png)

<center> **~~ Interjection ~~**  </center>

What is the Cox Proportional Hazards Model?  
- Survival time is described as a probability distribution  
- Hazard Ratio: ratio between chance for event in one group vs other group  
- “Proportional Hazards” means you can have multiple groups  
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
- Cox Proportional Hazards Model has drawbacks  
- Assumption that hazard ratio doesn’t change over time  
	+ Can be untrue when hazards diverge (i.e. different treatment courses)
- Even with Random Survival Forest, Cox PH Model was used for feature selection
- “Scaling” mortality to survival time seems a bit sketchy
	+ But then again, results speak for themselves
>*"Regression forests are for nonlinear multiple regression. They allow the analyst to view the importance of the predictor variables."*
>*"Survival forests are a model-free approach to survival analysis. They allow the analyst to view the importance of the covariates as the experiment evolves in time"*


*****

#### References 

