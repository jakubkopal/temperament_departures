# Code supplement to Early-childhood temperament deviations mark psychiatric risk into early adulthood
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F862615-informational
)]([https://doi.org/10.1101/862615](https://doi.org/10.1101/2022.04.23.489093))

This repository contains the code used to process and analyse the data presented in the "Early-childhood temperament deviations mark psychiatric risk into early adulthood" paper. 

## Abstract

<div style="margin-left: 40px;" align="justify">
Early-childhood temperament is associated with mental health outcomes decades later. Temperament reflects early-emerging individual differences in emotional and behavioral tendencies. These differences are relatively stable across development and shaped by both genetic and environmental influences. However, the consequences of departures from expected developmental trajectories remain largely unexplored. Using data from more than 50,000 children in the Norwegian Mother, Father and Child Cohort Study, we modeled longitudinal temperament trajectories at 1.5, 3, and 5 years of age and quantified deviations from expected development. Multivariate pattern analysis revealed latent dimensions linking these deviations to clinical diagnoses, with ADHD as the most prominent outcome. Time-to-event analysis showed that these dimensions were associated with a higher hazard of ADHD diagnosis across childhood and adolescence. Finally, genetic analyses identified loci jointly associated with temperament trajectories and ADHD, revealing age-dependent genetic effects. Together, these findings show that deviations from temperament trajectories in early childhood capture transdiagnostic vulnerability across development. Early temperament monitoring may thus serve as an indicator of later mental health risk.
</div>


<c>![Figure 1](https://github.com/jakubkopal/temperament_deviations/figures/fig1.png)</c>


## Resources and Scripts
Findings reported in the article are based on the analysis scripts in the `Scripts` folder.

1.   `Scripts/FEMA-Long_GWAS.m` - an example MATLAB script leveraging FEMA-Long to run longitudinal GWAS on EAS temperament traits.
2.   `Scripts/PLS_EAS_Diagnosis.ipynb` - performs canonical Partial Least Squares to link temperament deviations with psychiatric diagnoses. Includes permutation testing of the number of significant PLS dimensions and bootstrap testing of PLS loadings.
3.   `Scripts/Survival_PLS-score.ipynb` - performs time-to-event analysis using a latent temperament profile to predict age of ADHD diagnosis.
