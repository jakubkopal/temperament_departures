# Early-childhood temperament trajectories map onto transdiagnostic psychiatric risk
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F862615-informational
)]([https://doi.org/10.1101/862615](https://doi.org/10.1101/2022.04.23.489093))

This repository contains the code used to process and analyse the data presented in the "Departures from predicted temperament trajectories are associated with later psychiatric diagnoses" paper. 

## Abstract

<div style="margin-left: 40px;" align="justify">
Early-childhood temperament is associated with later mental health. Temperament continues to develop throughout the first years of life, and a single assessment cannot capture its trajectory. Whether departures from an individual's developmental trajectory carry psychiatric risk remains unknown. Using data from more than 50,000 children in the Norwegian Mother, Father and Child Cohort Study, we modeled longitudinal temperament at 1.5, 3, and 5 years of age with the FEMA-Long mixed-effects framework. We then quantified each child's departure from their predicted trajectory. Multivariate analysis revealed two transdiagnostic dimensions linking trajectory departures to psychiatric diagnoses across childhood and adolescence. Higher scores on the first dimension were associated with an increased hazard of subsequent ADHD diagnosis (hazard ratio = 1.54), and higher scores on the second with an increased hazard of Asperger syndrome (hazard ratio = 1.64). To examine the genetic basis of these associations, we performed longitudinal GWAS of temperament and conjunctional FDR analysis to detect loci shared between temperament and the associated disorders. SNP effects changed across early childhood, with some strengthening and others attenuating with age. These findings show that departures from predicted temperament trajectories reflect transdiagnostic psychiatric risk with a shared genetic basis. Longitudinal, trajectory-based monitoring could help identify children at elevated psychiatric risk.
</div>


<c>![Figure 1](https://github.com/jakubkopal/temperament_departures/blob/main/figures/fig5.png)</c>


## Resources and Scripts
Findings reported in the article are based on the analysis scripts in the `Scripts` folder.

1.   `Scripts/FEMA-Long_GWAS.m` - an example MATLAB script leveraging FEMA-Long to run longitudinal GWAS on EAS temperament traits.
2.   `Scripts/PLS_EAS_Diagnosis.ipynb` - performs canonical Partial Least Squares to link departures from temperament trajectories with psychiatric diagnoses. Includes permutation testing of the number of significant PLS dimensions and bootstrap testing of PLS loadings.
3.   `Scripts/Survival_PLS-score.ipynb` - performs time-to-event analysis using a latent temperament profiles to estimate the risk of ADHD diagnosis.
