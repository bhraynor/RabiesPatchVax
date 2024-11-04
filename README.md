# RabiesPatchVax
## Overview
This repo contains code for a metapopulation model for staggered, pulse vaccination 
campaigns. This code accompanies the manuscript Bellotti et al, 
"Challenging a paradigm: staggered versus single-pulse mass dog vaccination 
strategy for rabies elimination."

## Table of contents
**src folder** contains source code to run analyses for the manuscript  

- Patch_HierarchicalModel.Rmd is the Bayesian Hierarchical Model to estimate death 
rates of rabid dogs across different patches (microreds) with rendered output (html) 
files.

- Simulations_ToyModel.Rmd includes an analysis of the simplified toy model with 
rendered output (html) files.

-Historical_simulation_fitting.Rmd includes analysis to fit simulated data 
(included in manuscript supplement) to Arequipa case data with rendered output (html).

-Simulations_StochasticVaccinationScenarios.Rmd includes codes to run stochastic 
model simulations and output for the main analysis and sensitvity analysis (html files).

**data_minimal folder** contains the deidentified data sets needed to run the codes

**R folder** contains R functions called in the analysis codes

## Questions
- contact Brinkley Bellotti, brinkley.bellotti@gmail.com