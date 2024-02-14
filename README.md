# Estimating Crown-of-Thorns Starfish (CoTS) population size across Australia’s Great Barrier Reef (GBR)

This repository contains scripts and data to estimate census population size of *Acanthaster* cf. *solaris* (Crown-of-Thorns Starfish) on the GBR.
Code and data are needed to reproduce the results shown in:

[**Popovic I, Bergeron LA, Bozec Y-M, Waldvogel A-M, Howitt SM, Damjanovic K, Patel F, Cabrera MG, Wörheide G, Uthicke S, Riginos C** (2024) High germline mutation rates, but not extreme population outbreaks, influence genetic diversity in a keystone coral predator. *PLOS Genetics* 20(2): e1011129.](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011129)

## Instructions

1) In the folder **Calibration_CoTS_density_from_mantatow/** build the calibration model of the density of noncryptic CoTS from manta tow counts based on empirical data (*MTC_vs_SSC.csv*) using the R script *COTS_MTC_calibration.R*. The resulting model (R object of class ‘lm’) is stored in *Calibration_model.RData*.
  
2) In the folder **Bootstrap_generation/**, generate 500 bootstrap samples of reef-level CoTS manta-tow counts from the annual monitoring surveys (*LTM_COTS_MTC.mat*) using the Matlab script *CREATE_BOOTSTRAP.m*. The resulting bootstrap samples (Matlab array) are stored in *BOOTSTRAP_MTC.mat*.
   
3) Bootstrap confidence intervals of CoTS population size are then computed in the R script *COTS_BOOTSTRAP_CI_R*. The script also produces the graphics shown in Figure 5.

## Data sources
- CoTS calibration data (*MTC_vs_SSC.csv*) were extracted from:
  
    **Moran PJ, De'ath, G** (1992) Suitability of the manta tow technique for estimating relative and absolute abundances of crown-of-thorns starfish (Acanthaster planci L.) and corals. Marine and Freshwater Research, 43(2):357-379.
    
- CoTS survey data (*LTMP_COTS_MTC.csv*) were obtained from the Australian Institute of Marine Science's Long Term Monitoring Program (https://eatlas.org.au/gbr/ltmp-data)

- 3D geomorphic areas of individual reefs (*GEOM_3D_AREAS.csv*) were derived from:

    **Roelfsema et al.** (2021) How much shallow coral habitat is there on the Great Barrier Reef? Remote Sensing, 13(21):4343
    
    **Castro-Sanguino et al.** (2023) Control efforts of crown‐of‐thorns starfish outbreaks to limit future coral decline across the Great Barrier Reef. Ecosphere, 14:e4580.

- Spatial coordinates of the centroids of individual reefs and GBR management area designation from the ecosystem model ReefMod-GBR (*GBR_REEF_POLYGONS_2022.csv*): https://github.com/ymbozec/REEFMOD.6.8_GBR
  
    **Bozec et al.** (2022) Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92:e01494.

- Reef polygons from the GBR Reef Features shapefile available at: https://geoportal.gbrmpa.gov.au/
  
     **Great Barrier Reef Marine Park Authority** (1998) Great Barrier Reef Features (Version 1.4).
  
#
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)
