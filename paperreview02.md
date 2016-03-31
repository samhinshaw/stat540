### [FloReMi: Flow Density Survival Regression Using Minimal Feature Redundancy](http://www.dx.doi.org/10.1002/cyto.a.22734)

#### Published in *Cytometry Part A*, August 2015

*****

#### Overview
The goal of the FlowCAP IV challenge was to predict the time until progression to AIDS for HIV patients, a task which was manually studied by Ganesan, et al.<sup>[1](#references)</sup>

#### Data

#### Analysis
The FloReMi algorithm performs in four main steps.  
1. Preprocessing.  I will not touch on this further as this step is part of a standard flow cytometry workflow, and is not computationally-based.  
2. Feature Extraction.  This step is unsupervised machine learning, and sets up feature selection.  
3. Feature Selection. This step is a supervised machine learning method, and is key to the "minimal feature redundancy".  
4. Survival Time Prediction.  Finally, with the features selected, the Cox Proportional Hazards model is used to fit a logsitic regression model to our selected features, predicting survival time of patients.  

2. Feature Extraction
- Determine splits -flowDensity for automatic gating in one dimension.  
	+ Wonderfully optimized; no clustering!  
- Define subsets from thresholds determined by flowDensity:  
![Defined Subsets](./definedsubsets.png)
- Compute all 14 features for each subtype of each sample for both stimulated, unstimulated, and diff between stim & unstim. 
![Equation](./equation.png)

<center> ~~ Interjection ~~  </center>

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
- Cox PH model
- Random survival forest
- Additive hazards model

#### Critique


*****

#### References 

1. Ganesan A, Chattopadhyay PK, Brodie TM, Qin J, Gu W, Mascola JR, Michael NL,Follmann DA, Roederer M. Immunologic and virologic events in early HIV infectionpredict subsequent rate of progression. J Infect Dis 2010;201:272–284.