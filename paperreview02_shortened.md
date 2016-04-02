## [FloReMi: Flow Density Survival Regression Using Minimal Feature Redundancy](http://www.dx.doi.org/10.1002/cyto.a.22734)

### Published in *Cytometry Part A*, August 2015

*****

### Goals & Findings
This paper was published in response to the [FlowCAP IV challenge](http://flowcap.flowsite.org/), initiated to better predict clinical outcome of patients from blood draw samples.  Specifically, the goal of the FlowCAP IV challenge was to predict the time until progression to AIDS for HIV patients, given peripheral blood mononuclear cells (PBMC) analyzed in different conditions by flow cytometry. The FlowCAP consortium is run in part by the BC Cancer Agency and University of British Columbia.  
This paper, the winner of the challenge, found that their method of minimal feature redundancy worked well during cross validation when combined with predictive models such as the Cox Proportional-Hazards model, the Additive Hazards model, and the Random Survival Forests model.  However, upon testing testing novel data, the authors found that even with their method, only the Random Survival Forest model did significantly better than random, mostly due to its resilience to overfitting.  

### Data

The authors were provided with the FlowCAP IV challenge dataset, high-dimensional (multicolor) flow cytometry dataset with sixteen different markers:  

- FSC-A, FSC-H, SSC-A - three markers describing cells' size and shape  
- IFNγ, TNFα, CD4, CD27, CD107-A, CD154, CD3, CCR7, IL2, CD8, CD57, CD45RO, V-Amine/CD14 - various immune markers  
The samples were split into two groups: stimulated with HIV antigen and unstimulated. 
The authors were provided first with a training data set with which to develop their algorithm/pipeline, and later provided with a test dataset on which to test the efficacy of their algorithm.  

### Analysis

The FloReMi algorithm performs in four main steps.  

#### 1. Preprocessing.  

The authors needed to automate a standard flow cytometry workflow because of the high dimensionality of the data, and in order for their research to be reproducible.  As this is a standard flow cytometry workflow, I will not touch on it further.  All I will mention is that the authors used flowDensity, an automated gating program developed at the BC Cancer Research Agency.  

#### 2. Feature Extraction  
Next, the authors moved to unsupervised learning with feature extraction.  This step is composed of three main parts.  

- Determine Splits  
- Determine Subsets  
- Extract Features  
![Defined Subsets](definedsubsets.png)

#### 3. Feature Selection.  
This step is a supervised machine learning method, and is key to the "minimal feature redundancy" highlighted in this paper's title.  The idea was to select features with high correlation to survival time, but because of the nature of survival analysis, a simple metric such as pearson correlation could not be used due to censored values (patients that had not yet progressed to AIDS).  Therefore, the authors proceeded with the Cox proportional-hazards model.  

#### 4. Survival Time Prediction.  
Finally, with the features selected, the authors evaluated three different regression techniques to predict survival time (time to AIDS) of patients.  The authors performed leave-one-out cross validation on their training sets, creating a final model with the entire training dataset.  Finally, the model's performance was evaluated with p-values and concordance index.  The algorithm was tested against:

- Cox proportional hazards  
- Random survival forests  
- Additive hazards

### Critique
Interestingly, the model predicts that some of the best features for predicting survival time are "negative" features. Looking at **Table 3**, you can see that the vast majority of the subset identifiers are negative populations, with 2/13 of the subsets being entirely composed of negative markers.  This is worrisome, as it leads me to believe that there *are* other markers not studied here that better identify these populations.  Just off the top of my head, the exclusion of IFN-γ, IL-2, and TNFα in the feature extraction step could contribute to this problem.  

Continuing with that thought, the exclusion of IFN-γ, IL-2, and TNFα in Feature Extraction is particularly worrysome due to their functional importance.  IL-2 especially, which is a T cell growth factor, and of importance to T<sub>regs</sub>, may play a role in HIV. In order to reduce the computation time necessary when including these cytokines, the authors could explore novel algorithms for automatic gating aside from flowDensity, which may have better luck with distributions that are not cleanly bimodal.  Though these cytokines were still used in feature selection, none of them appear in the top 13 list used in the Cox PH & Random Survival Forests models.  From an experimental design viewpoint, it would be interesting to repeat this study with new markers/stains.  

Another weakness of this algorithm is its reliance on the Cox proportional-hazards model, which has its drawbacks.  Particularly, the Cox PH model assumes that the hazard ratio for any given patient does not change over time.  I believe this particular assumption is often untrue, as hazards may often diverge due to environmental/lifestyle factors, or simply due to change in treatment--which happens often with HIV patients.  

Even in the Random Survival Forests approach, the Cox PH model was used for feature selection.  There is still a possible weakness in the feature selection step here, as random forests are particularly susceptible to highly correlated data which distorts the randomness of trees.  Even though the authors attempted to control for highly correlated features, they only used a pearson correlation to this end.  I believe they could have looked at the identifiers of the subsets as well, and restricted their overlap.  For example, the first two features used were the percentage of cells in the unstimulated group with CD4/CD27 double negative population.  Though there were other differences between these groups, I find it hard to believe that they could be so different as to warrant being the top two features.  Perhaps the authors could have 

Finally, “scaling” mortality (0-1 scale) to survival time seems a bit sketchy.  That being said, though it is less than ideal, the results speak for themselves.  However, it would be interesting if this method could be adapted to allow for greater performance, as random forests seem to have the best performance overall.  


*****

