# Tuning Random Forests for Causal Inference Under Cluster-Level Unmeasured Confounding

Youmi Suk and Hyunseung Kang

## Overview

This paper focuses on a machine learning (ML) method based on random forests known as Causal Forests and presents five simple modifications for tuning Causal Forests so that they are robust to cluster-level unmeasured confounding. Our simulation study finds that adjusting the default tuning procedure wit the propensity score from fixed effects logistic regression or using variables that are centered to their cluster means with a corresponding propensity score produces estimates that are more robust to cluster-level unmeasured confounding. Also, when these parametric propensity score models are mis-specified, our modified ML methods remain robust to bias from cluster-level unmeasured confounders compared to existing parametric approaches based on propensity score weighting. We conclude by demonstrating our proposals in a real data study concerning the effect of taking an eighth-grade algebra course on math achievement scores from the Early Childhood Longitudinal Study (ECLS). 

Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using the kindergarten cohort of ECLS (ECLS-K) data. 

## Simulation Study

* `DGP_clusterOVB.R`  

   This `R` file includes data generating codes for multilevel observational data with cluster-level unmeasured confounding.
 
```R
twolevel.normalU
twolevel.uniformU
twolevel.binary 
twolevel.experror 
```

* `Simulation_robustCF.R`
 
   This `R` file includes simulation codes with our proposed modifcations for Causal Forests. In addition to different DGPs, one can change other simulation parameters: `smpl.size`, `o.val`, and `m.val`.  For more details on simulation condtions, see our paper, https://psyarxiv.com/36w72/.


## ECLS-K Data Study

* `ECLSK_Algebra_complete.csv`

  This is our complete data. The original ECLSK 1998-99 dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. For more information, see [Tourangeau et al. (2009)](https://nces.ed.gov/pubs2009/2009003.pdf).

* `ECLSK_Algebra_robustCF.R` 
 
   This `R` file can be used to replicate our data analysis.
