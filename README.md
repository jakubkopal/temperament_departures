# Departures from predicted temperament trajectories are associated with later psychiatric diagnoses
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F862615-informational
)]([https://doi.org/10.1101/862615](https://doi.org/10.1101/2022.04.23.489093))

This repository contains the code used to process and analyse the data presented in the "Departures from predicted temperament trajectories are associated with later psychiatric diagnoses" paper. 

## Abstract

<div style="margin-left: 40px;" align="justify">
Early-childhood temperament is associated with later mental health, yet it continues to develop across the first years of life, so single measurements may miss information in its trajectory. Whether departures from typical developmental trajectories carry psychiatric risk beyond static trait levels remains unknown. Using data from more than 50,000 children in the Norwegian Mother, Father and Child Cohort Study, we modeled longitudinal temperament at 1.5, 3, and 5 years of age with the FEMA-Long mixed-effects framework. We then quantified child's departure from the predicted trajectories. Multivariate analysis revealed two transdiagnostic dimensions linking the departures to psychiatric diagnoses across childhood and adolescence. Higher scores on the first dimension were associated with increased hazard of subsequent ADHD diagnosis (hazard ratio = 1.54), and higher scores on the second with an increased hazard of Asperger syndrome (hazard ratio = 1.64). To examine the genetic basis of these associations, we performed longitudinal GWAS of temperament and conjunctional FDR analysis to detect loci shared with ADHD. Examining how each SNP's effect changed across early childhood, we found that some strengthened and others attenuated with age. These findings show that departures from predicted temperament trajectories capture transdiagnostic psychiatric risk with a shared genetic basis. Longitudinal, trajectory-based monitoring could help identify children at elevated psychiatric risk.
</div>


<c>![Figure 1](https://github.com/jakubkopal/temperament_departures/blob/main/figures/fig5.png)</c>


## Resources and Scripts
Findings reported in the article are based on the analysis scripts in the `Scripts` folder.

1.   `Scripts/FEMA-Long_GWAS.m` - an example MATLAB script leveraging FEMA-Long to run longitudinal GWAS on EAS temperament traits.
2.   `Scripts/PLS_EAS_Diagnosis.ipynb` - performs canonical Partial Least Squares to link temperament deviations with psychiatric diagnoses. Includes permutation testing of the number of significant PLS dimensions and bootstrap testing of PLS loadings.
3.   `Scripts/Survival_PLS-score.ipynb` - performs time-to-event analysis using a latent temperament profile to predict age of ADHD diagnosis.
