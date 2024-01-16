## Calibration model of the density of noncryptic CoTS from manta tow counts.
## This re-visits the regression model of Moran & De'ath (1992) AUST MAR FRESHW
## to provide an operational tool enabling estimation of the density of Acanthaster
## starfish (CoTS) from manta tow counts. This essentially re-builds their 
## regression model of SCUBA swim counts (SSC) vs manta tow counts (MTC)
## and adds uncertainty around predictions (consistent with the model and data)
## based on the number of tows performed around the perimeter of a reef. 

## -----------------------------------------------------------------------------
## Yves-Marie Bozec, The University of Queensland, Australia (y.bozec@uq.edu.au)
## March 2023
## -----------------------------------------------------------------------------

library(ggplot2)
library(ggthemes)

rm(list=ls())
graphics.off()

## Number of CoTS recorded during manta tows and SCUBA swim searches
## (extracted from Moran & De'ath 1992, table 7)(extracted from Moran & De'ath 1992, table 7)
DATA = read.table("MTC_vs_SSC.csv", header=TRUE, sep="\t")

# Remove the totals (last row)
DATA = DATA[-31,]

# Rename columns fro convenience
names(DATA)[4]="MTC" # manta tow counts (total number of COTS detected on manta tows)
names(DATA)[5]="SSC" # swim search counts (total number of COTS detected on the equivalent area)

# Calculate mean counts per tow
DATA$MTC_per_tow = DATA$MTC/DATA$Number.of.tows
DATA$SSC_per_tow = DATA$SSC/DATA$Number.of.tows

## Reproduce Fig 5a & 5b of Moran & De'ath (1992) to check we've got the right data
plot(DATA$MTC, DATA$SSC,pch=19,cex=1.5,xlim=c(0,620),ylim=c(0,1400), 
     xlab='Manta Tow Count (MTC)',ylab='Swim Search Count')

plot(DATA$MTC_per_tow^(1/3), DATA$SSC_per_tow^(1/3),pch=19,cex=1.5,xlim=c(0,5),ylim=c(0,7), 
     xlab='Cube Root MTC',ylab='Cube Root SSC')

## Test the relationship between transformed preditor/response
## Linear model using number of tows as weights
DATA$X = DATA$MTC_per_tow^(1/3)
DATA$Y = DATA$SSC_per_tow^(1/3)
M1 = lm(Y~X, data=DATA,weights = Number.of.tows)

## Check summary stats are ~ the same as the published ones (ie, regression model is the same)
R2 = summary(M1)$r.squared # 0.913 vs 0.911 in Moran & De'ath
m = summary(M1)$coefficients[2] # 1.2008 vs 1.2008 in Moran & De'ath
p = summary(M1)$coefficients[1] # 0.8072 vs 0.8071 in Moran & De'ath

## -----------------------------------------------------------------------------
## Use the model for predicting SSC per tow (ie, for an equivalent area ~ 200 m x 12 m)
## from MTC expressed in a per tow basis (ie, as averaged across a reef)
## -----------------------------------------------------------------------------
# Use type='response' (not 'terms') to get predictions comparable to observations
# se.fit is the standard error for the mean prediction at given values of the predictor variables.
NEW = data.frame(X=seq(0,5,by=0.01))
NEW$pred = predict(M1, type='response',newdata=NEW, se.fit = TRUE)$fit
NEW$se.fit = predict(M1, type='response',newdata=NEW, se.fit = TRUE)$se.fit
NEW$res.scale =  predict(M1, type='response',newdata=NEW, se.fit = TRUE)$residual.scale

# PLOT THE CALIBRATION MODEL
# Relationship between manta tow counts and SCUBA swim counts from Moran and De’ath (1992a). 
# The observed data (black dots) are expressed as a per-tow basis (200×12 m, with 2 to 15 tows/transects supporting each observation).
# The regression model is used to generate deterministic predictions of SCUBA swim counts (line) or stochastic predictions 
# that are function of the number of tows conducted in a given area (e.g. the perimeter of a reef), illustrated here with 5 (red dots)
# and 25 tows (green dots). Increasing the number of tows decreases the dispersion around the deterministic model.
plot1 <-  ggplot(DATA, aes(x=MTC_per_tow^(1/3), y=SSC_per_tow^(1/3))) + theme_bw() +
  # Predictions if all samples were obtained with N=5 tows (this generates some variability around the deterministic prediction):
  geom_point( data=NEW, aes(x=X, y=rnorm(length(pred),pred, sqrt(se.fit^2+res.scale^2)/sqrt(5))), size=2.5,shape=19,col='firebrick',alpha=0.5)+
  # Predictions if all samples were obtained with N=25 tows (greater span around the deterministic prediction):
  geom_point( data=NEW, aes(x=X, y=rnorm(length(pred),pred, sqrt(se.fit^2+res.scale^2)/sqrt(25))),size=2.5,shape=19,col='turquoise4',alpha=0.5)+
  # Deterministc prediction:
  geom_line( data=NEW, aes(x=X, y=pred), linewidth = 1, color = 'gray20') +
  geom_point(size=5, fill='gray20',color='white',shape=21) +
  scale_y_continuous(name = "SCUBA swim count (cube root)", limits = c(0, 7)) + 
  scale_x_continuous(name = "Manta tow count (cube root)", limits = c(0, 5)) + 
  theme(axis.ticks = element_line(linewidth = 0.25),axis.text = element_text(size=12), axis.title = element_text(size=16))

ggsave('CALIBRATION.png', plot1, width = 8, height = 6)

## -----------------------------------------------------------------------------
## Export the calibration model for routine use
## -----------------------------------------------------------------------------
save(M1, file = "Calibration_model.RData")

# NOTES ON HOW TO USE THE CALIBRATION MODEL:
# 1) use model M1 to predict from average CoTs per tow (over N tows) the abundance of COTS as in SCUBA search swims over 200m x 12m
# 2) for uncertainty in predictions, just generate at random (normal distribution) with mean=.$fit and sd=sqrt(.$se.fit^2+.$residual.scale^2)/sqrt(N)
#   So, for a given reef-wide value of MTC (per tow) obtained with N tows:
#   a) PRED = predict(M1, type='response', newdata=MTC, se.fit = TRUE)
#   b) SSC = rnorm(1,PRED$fit, sqrt(PRED$se.fit^2+PRED$residual.scale^2)/sqrt(N))
# In doing so, the predictions accounts for the number of tows that would be used, ie the lower n the greater the variability in predictions
# In the plot above, the mean number of tows per reef was n=5, and with this value the random generation gives
# similar levels of variability around the fit.
# Output is a random estimate of absolute density of noncryptic COTS for a 200m x 12m equivalent area
# 3) don't forget to back-transform the prediction (cubic)
