# BayesSenMC
R package for accounting for misclassification with Bayesian methods
Author: James (Jinhui) Yang, Haitao Chu, Lifeng Lin

This package is based on the methods proposed in 
1) Chu et al.'s paper (2006): https://www.ncbi.nlm.nih.gov/pubmed/16843678, which introduces a Bayesian approach to deal with misclassification of exposure in a case-control study
2) Chu et al.'s paper (2010): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3035476/, which introduces a Generalized Linear Mixed Effects Model to provide parameter estimates of Se and Sp based on meta-reviews of diagnostics accuracy.

BayesSenMC provides users with 6 different models (as outlined in Chu's paper) that estimates the adjusted odds ratio in a case-control study based on the different priors of sensitivity and specificity as well as historic data of related studies. 

The models with zero and constant misclassification (constant Se and Sp) are able to be used given the 2x2 table of a case-control study. However, other models require parameters from a NLMIXED procedure. Similarly, the graphing function can be used with stanfit objects that are returned from any of the model functions.
