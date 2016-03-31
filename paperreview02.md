### [FloReMi: Flow Density Survival Regression Using Minimal Feature Redundancy](http://www.dx.doi.org/10.1002/cyto.a.22734)

#### Published in *Cytometry Part A*, August 2015

*****

#### Overview
The goal of the FlowCAP IV challenge was to predict the time until progression to AIDS for HIV patients, a task which was manually studied by Ganesan, et al.<sup>1</sup>

#### Data

#### Analysis
The FloReMi algorithm performs in four main steps.  
1. Preprocessing.  I will not touch on this further as this step is part of a standard flow cytometry workflow, and is not computationally-based.  
2. Feature Extraction.  This step is unsupervised machine learning, and sets up feature selection.  
3. Feature Selection. This step is a supervised machine learning method, and is key to the "minimal feature redundancy".  
4. Survival Time Prediction.  Finally, with the features selected, the Cox Proportional Hazards model is used to fit a logsitic regression model to our selected features, predicting survival time of patients.  

#### Critique

*****
1. Ganesan A, Chattopadhyay PK, Brodie TM, Qin J, Gu W, Mascola JR, Michael NL,Follmann DA, Roederer M. Immunologic and virologic events in early HIV infectionpredict subsequent rate of progression. J Infect Dis 2010;201:272–284.