# BayesSenMC
R package for accounting for misclassification with Bayesian methods
Author: James (Jinhui) Yang

This package is based on the methods proposed in Chu et al.'s paper (2013): https://www.ncbi.nlm.nih.gov/pubmed/16843678
BayesSenMC provides users with 6 different models (as outlined in Chu's paper) that estimates the adjusted odds ratio in a case-control study based on the different priors of sensitivity and specificity as well as historic data of related studies. 

The models with zero and constant misclassification are able to be used given the 2x2 table of a case-control study. However, other models require parameters from a NLMIXED procedure, which has not been implemented yet (TODO). Similarly, the graphing function can be used with input of stanfit objects that are returned from any of the model functions.
