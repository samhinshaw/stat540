# Empirical Bayes, Linear Models, and Limma

Look at model residuals to see if you've underfit.  If residuals are too low, you might have overfit?  
- F test (nested models)  
- Likelihood Ratio Test (LRT) (nested models)  
- AIC (do not need nested model)  
- BIC (do not need nested model)  

MODEL SELECTION <- important

{base} function `lm()` can do "multivariate multiple regression"
Don't do this!!  This is very computationally inefficient.  

What do we need to do large scale inference?  
- Moderated t- and F-statistics (limma) 
- empirical Bayes (eBayes)  
Limma also has methods for moderating multiple testing  

BUT what's so special about these moderated statistics?
In too many experiments, n is low, but dimensionality is FAR higher.  
Volcano plot!
- y-axis is log2 fold change  
- x-axis is -log10 of p-value  
- visually represents differential expression AND significance  
	+ so, may have highly differentially expressed genes that is not significant  
- in particular, some genes will be very significant (large -log10 p value, or small p-value), but are not biologically relevant.  these genes will have low fold change, but high significance.  
- This is why we use a 2FC(log2) threshold on the data.  

We are underestimating variance, which puts greater importance on these genes that have falseley low variance, but also very little change in expression. 
SO, we can "moderate" the variance by adding a constant, or even weight variance (weighting the constant (So) more or the true variance (Sg) more)  
This moderated t-test statistic is derived from the empirical bayes formula  
- estimated covariance matrix?  
- contrast matrix!  
use contrast matrix to go from mean-centering parameterization to your treatment-effect parameterization  

Bayesian statistics!  
- we want to incorporate prior knowledge  
- maximum likelihood estimation versus maximum a *posteriori* (from Bayesian statistics)  
- <Insert Bayes' Theorem Here>  
- We need to come up with a prior distribution! Obviously we don't know this.  BUT, we can use some formula in the slides, but those require their own parameters, which we'll call **hyperparameters**.  

SO, what's our moderated variance?
It's a combination of our prior evidence and our observed evidence (true variance)
If d0 is 0, then you're ignoring your prior
If d0 is closer to infinity, you're completely ignoring your data and only going with the prior. (making DE calls based on fold-change alone?)

if d0 is very high, you're "overwriting" the variance for every gene

originally have per-gene variance based on residuals
limma gives you variance across all of your genes/samples

SAM - significance analysis of microarrays
Check out linear models book by Robert Gentleman

R code!

```r
library(limma)
```

We need our design matrix and contrast matrix.  The contrast matrix allows us to determine which samples to compare.  FORTUNATELY, we can easily generate this
- can also just be an identity matrix.

Next week is multiple test correction
- in particular analysis of image analysis, fMRI
