## Compute and plot (Figure 5) the bootstrap confidence intervals of CoTS population
## size annually from 1991 to 2022.
##
## First, the 500 bootstrapped annual samples of reef-level manta tow counts (compiled
## in 'BOOTSTRAP_MTC.mat') are converted into density estimates using the calibration
## model ('Calibration_model.RData'). Then, mean density is calculated for each annual
## bootstrap sample and turned into population size after multiplication with the total
## 3D area of CoTS habitat estimated from geomorphic maps (Roelfsema et al. 2021).
## Confidence intervals of total population size are determined by taking the 2.5th 
## and 97.5th percentiles of each annual distribution between 1991 and 2022.

## -----------------------------------------------------------------------------
## Yves-Marie Bozec, The University of Queensland, Australia (y.bozec@uq.edu.au)
## March 2023
## -----------------------------------------------------------------------------

library(ggplot2)
library(plotrix)
library(ggthemes)
library(patchwork); 
library(R.matlab)
library(comparator) # for the harmonic mean function
library(sf) # for mapping

rm(list=ls())
graphics.off()

## -----------------------------------------------------------------------------
# Load the calibration model to convert COTS per tow into abundance per unit area
# (M1 of class 'lm')
load('Calibration_CoTS_density_from_mantatow/Calibration_model.RData')

# Load the reef-level manta tow counts (AIMS Long Term Monitoring Program) (2649 x 5)
# This gives the mean number of CoTS per tow for each monitored reef between 1991 and 2022
# Australian Institute of Marine Science (AIMS). (2015). AIMS Long-term Monitoring Program:
# Crown-of-thorns starfish and benthos Manta Tow Data (Great Barrier Reef). https://doi.org/10.25845/5c09b0abf315a.
LTMP_MTC = read.table('LTMP_COTS_MTC.csv', header=T, sep=",")
# Dataframe with as many rows (2649) as reef-level observations per year
# SampledYear: year of monitoring (from original tow date)
# ReefID: Reef IDs beginning 10- to 24- are generally within the Marine Park, 
#      those beginning 26- and 27- are from the Moreton Bay area and those beginning 99- are from Torres Straight
# MngmtRegion: GBRMPA designation of management areas (Far Northern, Cairns/Cooktown, Townsville/Whitsunday , Mackay/Capricorn 
# N: number of conducted tows
# MEAN_COUNT: mean count per tow (number of counted CoTS divided by the number of tows)
                                                          
# Load the bootstrap samples of reef-level COTS per tow (2649 x 500)
BOO_MTC = readMat('Bootstrap_generation/BOOTSTRAP_MTC.mat')$BOOTSTRAP.COUNTS

# ------------------------------------------------------------------------------------------
# Estimate CoTS density from manta tow counts
# ------------------------------------------------------------------------------------------
# Use model M1 to predict, from reef-level CoTs per tow (over N tows), the abundance of COTS 
# in an equivalent 200m x 12m path surveyed by SCUBA swim search
# For the observed counts (MEAN_COUNT), just use deterministic predictions from model M1:
#   PRED = predict(M1, type='response', newdata=MEAN_COUNT, se.fit = TRUE)
# This gives a determistic estimate of absolute density of noncryptic COTS for a 200m x 12m area
#
# For the bootstrap samples (BOO_MTC), generate stochastic predictions as obtained with N tows:
#   PRED = predict(M1, type='response', newdata=BOO_MTC, se.fit = TRUE)
#   rnorm(1,PRED$fit, sqrt(PRED$se.fit^2+PRED$residual.scale^2)/sqrt(N))
# This gives a random estimate of absolute density of noncryptic COTS for a 200m x 12m area
#
# Predictions must be back-transformed (cubic), and scaled to a 1 km2 area
# WARNING: takes ~15 min to process!!

# First process the observed counts (using the deterministic model)
D_obs = matrix(NaN,length(LTMP_MTC$MEAN_COUNT),1) # observed (predicted) CoTS density per km2

for (r in (1:length(LTMP_MTC$MEAN_COUNT)))
{
  PRED = predict(M1, type='response', newdata=data.frame(X=LTMP_MTC$MEAN_COUNT[r]^(1/3)), se.fit = TRUE)
  D_obs[r,1] = (1e6/(200*12))*PRED$fit^3
}

# Then process the boostrap samples
D_boo = matrix(NaN, nrow(BOO_MTC),ncol(BOO_MTC)) # bootstrapped (predicted) CoTS density per km2

for (c in (1:ncol(BOO_MTC)))
{
  print(paste('processing bootstrap distri #',c,sep=''))
  for (r in (1:nrow(BOO_MTC)))
  {
    PRED = predict(M1, type='response', newdata=data.frame(X=BOO_MTC[r,c]^(1/3)), se.fit = TRUE)
    PRED_ran = rnorm(1,PRED$fit, sqrt(PRED$se.fit^2+PRED$residual.scale^2)/sqrt(LTMP_MTC$N[r]))
    D_boo[r,c] = (1e6/(200*12))*PRED_ran^3
  }
}
rm(PRED,PRED_ran, r, c)

# Replace by 0 the negatives in the stochastic boo densities (can happen with the generated noise)
for (i in (1:500))
{
  IDneg = which(D_boo[,i]<0)
  D_boo[IDneg,i]=0
}

# Collate with column of sampled years
D_boo = cbind(LTMP_MTC$SampledYear,D_boo)
# As many rows as reef-level observations per year
# Column #1 gives sampled year, #2 to #501 the stochastic densities (CoTS per km2) from bootstrap-generated MTC

# ------------------------------------------------------------------------------------------
# Bootstrap CI of reef-level density per year
# ------------------------------------------------------------------------------------------
# Estimate CI of total population over time.
# This takes all the bootstrap samples, calculate the mean density for each year
# then multiplies by total reef area of the GBR. This gives one distribution of 500 estimates
# of total population size from which a CI is calculated each year (percentile CI).

# List al the surveyed years
ALL_YEARS = unique(LTMP_MTC$SampledYear)

# Load estimates of 3D geomorphic habitat (Roelfsema et al 2021, Castro-Sanguino et al 2023)
# per management area and shelf position.
# $all_geom: cumulative 3D area of all geomorphic habitat classes (10 classes)
# $coral_geom: cumulative 3D area of coral-suitable geomorphic classes (4 classes)
GEOM_3D_AREAS = read.table('Spatial_data/GEOM_3D_AREAS.csv', header=T, sep=",")

# Calculate total geomorphic 3D area across the 4 management regions
# We only consider the coral-suitable geomorphic habitats (coral_geom)
TOT_GEOM = sum(GEOM_3D_AREAS$coral_geom) # for 3,104 reefs (including inshore)

## Function to calculate CI limits of population size in the entire system
myCIpop <- function(ALL_YEARS, DATA_BOO, TOT_AREA)
{
  # ALL_YEARS is the list of unique years of monitoring
  # DATA_BOO the bootstrapped densities with the first column giving the year of sampling
  # TOT_AREA is the total area (km2) of reef habitat
  CI_POP = matrix(NA,ncol=2,nrow=32)
  MEAN = matrix(NA,ncol=1,nrow=32) 
  MEDIAN = matrix(NA,ncol=1,nrow=32)
  
  for (id_year in (1:32))
  {
    I = which(DATA_BOO[,1]==ALL_YEARS[id_year]) #identify the rows corresponding to the current year
    # calculate mean reef-level densities (colMeans) for each bootstrap sample
    Mean_Dens = colMeans(as.matrix(DATA_BOO[I,2:501])) # 500 estimates of mean reef-level density (GBR wide)
    # then multiply by total area to get 500 estimates of total abundance
    ALL_BOOT = Mean_Dens*TOT_AREA

    CI_POP[id_year,] = quantile(ALL_BOOT,probs=c(0.025, 0.975)) # percentile CI
    MEAN[id_year,] = mean(ALL_BOOT)
    MEDIAN[id_year,] = quantile(ALL_BOOT,probs=0.5)
  }

  X = as.data.frame(cbind(ALL_YEARS,CI_POP,MEAN,MEDIAN))
  names(X)=c("YEARS", "CI0025", "CI0975","MEAN","MEDIAN") 
  return(X)
}

# Compute CI of total population
CI_POP_GEOM = myCIpop(ALL_YEARS, D_boo, TOT_GEOM)
# First column gives the 2.5th percentile of the 500 annual estimates of population size
# Second column the 97.5th percentile
# Third column the mean population size (super close to the observed mean)
# Fourth column is the median


# ------------------------------------------------------------------------------------------
# SUMMARY STATS ACROSS ALL YEARS
# ------------------------------------------------------------------------------------------
# Calculate the harmonic means of the CI limits across all years of monitoring
# These are expressed in millions starfish

HARMONIC_MEAN_GEOM = data.frame(CI0025 = round(hmean(CI_POP_GEOM[,2])/1e6,1),
                                CI0975 = round(hmean(CI_POP_GEOM[,3])/1e6,1),
                                MEAN = round(hmean(CI_POP_GEOM[,4])/1e6,1),
                                MEDIAN = round(hmean(CI_POP_GEOM[,5])/1e6,1))


# ------------------------------------------------------------------------------------------
# Plot confidence limits (FIGURE 5B)
# ------------------------------------------------------------------------------------------

P2 <-ggplot(data=CI_POP_GEOM, aes(x=YEARS)) +
  geom_point(aes(x=YEARS, y=MEDIAN/1e6),size=4) +
  geom_errorbar(aes(ymin=CI0025/1e6, ymax=CI0975/1e6), width=0) +
  geom_hline(yintercept = HARMONIC_MEAN_GEOM$CI0025, linetype = 2, linewidth = 0.8, col ='DarkCyan', alpha = 0.8) +
  geom_hline(yintercept = HARMONIC_MEAN_GEOM$CI0975, linetype = 2, linewidth = 0.8, col ='Chocolate', alpha = 0.8) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 0.25, colour='black'), axis.text = element_text(size=14, colour='black'), axis.title = element_text(size=22)) +
  theme(panel.grid = element_line(colour = 'white'), panel.background = element_rect(linewidth = 0.5, colour='black')) +
  labs(x="",y="Population size (million starfish)") +
  theme(axis.title.y = element_text(margin=margin(r=20)), plot.margin = margin(t = 20,  r = 20, b = 10, l = 30)) + # Top, Right, Bottom, Left margins
  ggtitle('B') + theme(plot.title = element_text(size = 20, face = 'bold'))

# Mean density and 95% confidence interval for every year of monitoring
# blue and dashed line: harmonic mean of the 2.5th percentile of mean annual density
# red dashed line: harmonic mean of the and 97.5th percentile

# ------------------------------------------------------------------------------------------
# Map of the GBR with monitored reefs (FIGURE 5A)
# ------------------------------------------------------------------------------------------

# Load the coordinates (centroids) of all GBR individual reefs (Bozec et al. 2022)
GBR_REEFS = read.table("Spatial_data/GBR_REEF_POLYGONS_2022.csv", header=T, sep=",")

monitored_reefs = unique(LTMP_MTC[,2:3])
REEF_POSITIONS = subset(merge(monitored_reefs, GBR_REEFS, all.x=T), select=c(LAT,LON))

# Map of entire QLD with reefs
# Great Barrier Reef Marine Park Authority 1998, Great Barrier Reef Features (Version 1.4).
# Retrieved from https://geoportal.gbrmpa.gov.au/.
myMap_shp = read_sf('Spatial_data/Great_Barrier_Reef_Features.shp')
head(myMap_shp)

# Merge the map with Management Areas defined in GBR_REEFS (=AREA_DESCR)
myReef_shp = subset(merge(myMap_shp, GBR_REEFS, all.x=T), select=AREA_DESCR)
select_reefs = which(myReef_shp$AREA_DESCR=='out'|myReef_shp$AREA_DESCR=='Far Northern'|
                     myReef_shp$AREA_DESCR=='Cairns/Cooktown'|myReef_shp$AREA_DESCR=='Townsville/Whitsunday'|
                     myReef_shp$AREA_DESCR=='Mackay/Capricorn')

# Locate main QLD towns
Towns_X = c(145.755504, 145.46191, 145.2502, 146.815703,  145.4621, 149.181966, 150.509423, 151.250934)
Towns_Y = c(-16.92745,  -16.48261, -15.4776, -19.262102, -14.6687, -21.144880, -23.381031, -23.843128)
Towns_XY = data.frame(LON = Towns_X,LAT = Towns_Y)
Towns_XY$Lab = c('Cairns', 'Port Douglas', 'Cooktown', 'Townsville', 'Lizard', 'Mackay', 'Rockhampton', 'Gladstone')

bkgd_col = 'Beige'

MyMap = ggplot(myMap_shp) +
  geom_sf(data=myMap_shp, fill = bkgd_col, linewidth=0.1) +
  geom_sf(data=myReef_shp[select_reefs,], color ='SteelBlue', fill='SteelBlue', linewidth=0.1) +
  geom_point(data = REEF_POSITIONS, x=REEF_POSITIONS$LON, y=REEF_POSITIONS$LAT, fill = 'orange', color='black', size=1.5,shape=21) +
  geom_point(data = Towns_XY[-c(2,3,5,8),], x = Towns_XY$LON[-c(2,3,5,8)], y = Towns_XY$LAT[-c(2,3,5,8)], fill = 'white', size=1, shape=24) +
  annotate(geom = 'text', x = Towns_XY$LON[-c(2,3,5,8)]-0.2, y = Towns_XY$LAT[-c(2,3,5,8)],
           label = Towns_XY$Lab[-c(2,3,5,8)], fontface = 'italic', color = 'black', size =4, hjust = 1) +
  annotate(geom = 'text', x = 151.7, y = -11.5, label = "Australia's Great Barrier Reef", fontface = 'italic', color='SteelBlue',
           size = 6.5, hjust = 1) +
  coord_sf(xlim = c(142.1,153.5), ylim = c(-25.3,-10.3), expand = FALSE) +
  scale_y_continuous(name='',breaks = seq(-28, -10, by = 4)) + 
  scale_x_continuous(name='',breaks = seq(140, 160, by = 4)) + 
  theme(panel.grid.major = element_line(linewidth = 0.01)) + 
  theme_bw() +
  theme(panel.background = element_rect(fill = 'LightBlue'))+
  theme(axis.ticks = element_line(linewidth = 0.25, colour='black'),axis.text = element_text(size=14, colour='black')) +
  theme(panel.grid = element_line(linewidth = 0, linetype = 0, colour='red'), panel.background = element_rect(linewidth = 0.5, colour='black')) +
  geom_rect(aes(xmin=142.5, xmax=146.8, ymin=-24.9, ymax=-23.55), fill='white', linewidth = 0.1, colour='black') +
  geom_point(x=143, y=-24.5, fill = 'orange', color='black', size=2,shape=21) +
  geom_point(x=143, y=-24, fill = 'SteelBlue', color='SteelBlue', size=2,shape=22) +
  annotate(geom = 'text', x = 143.3, y = -24.5, label = "CoTS monitoring", size =5, hjust = 0, vjust = 0.5) +
  annotate(geom = 'text', x = 143.3, y = -24, label = "individual reefs", size =5, hjust = 0, vjust = 0.5) +
  ggtitle('A') + theme(plot.title = element_text(size = 20, face = 'bold'))
  

# ------------------------------------------------------------------------------------------
# Combine map and plot of population sizes (FIGURE 5)
# ------------------------------------------------------------------------------------------
ggsave('FIGURE_5.png', plot=MyMap + P2, width = 15, height = 9)
