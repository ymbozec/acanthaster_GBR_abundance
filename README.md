# Estimating Crown-of-Thorns Starfish (CoTS) population size across Australia’s Great Barrier Reef (GBR)

This repository contains scripts and data to estimate the population size of Acanthaster (Crown-of-Thorns Starfish) on the GBR.
Code and data are needed to reproduce the results shown in:

**Popovic I, Bergeron LA, Bozec Y-M, Waldvogel A-M, Howitt SM, Damjanovic K, Patel F, Cabrera MG, Wörheide G, Uthicke S, Riginos C** (2024) High germline mutation rates, but not extreme population outbreaks, influence genetic diversity in a keystone coral predator. *PLOS Genetics*.

# Instructions

1) In the folder **Calibration_CoTS_density_from_mantatow/** build the calibration model of the density of noncryptic CoTS from manta tow counts based on empirical data (*MTC_vs_SSC.csv*) using the R script *COTS_MTC_calibration.R*. The resulting model (R object of class ‘lm’) is stored in *Calibration_model.RData*.
  
2) In the folder **Bootstrap_generation/**, generate 500 bootstrap samples of reef-level CoTS manta-tow counts from the annual monitoring surveys (*LTM_COTS_MTC.mat*) using the Matlab script *CREATE_BOOTSTRAP.m*. The resulting bootstrap samples (Matlab array) are stored in *BOOTSTRAP_MTC.mat*.
   
3) Bootstrap confidence intervals of CoTS population size are then computed in the R script *COTS_BOOTSTRAP_CI_R*. The script also produces the graphics shown in Figure 5.

# 
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)
