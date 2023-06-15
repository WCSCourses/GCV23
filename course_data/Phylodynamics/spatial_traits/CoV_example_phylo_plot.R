# basic script to plot BEAST mcc trees with traits
# S J Lycett
# 15 June 2023

# this is using my custom code- but I will release it as a GPLv3 package at some point soon

# only need to do this once
#install.packages("ape")
#install.packages("maps")
#install.packages("mapdata")
#install.packages("mapproj")
#install.packages("conicfit")

# base trees
library(ape)

# maps
library(maps)
library(mapdata)
library(mapproj)

# needed to fit HPDs to ellipses
library(conicfit)

# SJL specific - set wd to where the R scripts and data are
#setwd("~/Documents/Conferences/2023-06-12_Wellcome/Phylodynamics/SARS-CoV-2_example_files/spatial_traits/")

# Lycett_phylo utility code
source("getEl.R")
source("calcDecimalDate.R")
source("read_MCC_tree.R")
source("get_BEAST_cols.R")
source("custom_map_movie.R")

# define tree file name
treeName <- "cov_net_sim_mper2_120genomes_TN93G4_strict_skygrid_traits_2_mcc.tre"

# read in tree using custom function (could use ggtree also)
tr <- read_latlon_mcc_tr(treeName)

# manual correct of decimal dates due to dates format (fixing warning message)
dateTxt <- apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="\\|", fromEnd=TRUE)
decDates<- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE))
tr$decDates <- decDates
tr$youngestTip <- max(decDates)
tr$nodeDists <- NULL
tr$nodeTimes <- NULL
tr <- nodeTimes(tr, youngestTip=tr$youngestTip)

# add traits to the tree object
tr <- addDiscreteTraits(tr)
tr <- addLatLonHPD(tr)
tr <- fit_HPDs_to_standard(tr,ltol=0.005)
tr <- addFullTraitSet(tr, propIndex=which(tr$propNames=="Country"))

# save tree object for later use etc
save(tr, file=paste0(gsub(".tre","",treeName),"_with_hpds.tr.Rdata",sep=""))

# plot tree on map (all time)
xlim=c(-15,20)
ylim=c(40,65)
op <- par(mar=c(0,0,0,0))
plot_mcc_tree_with_hpds(tr, propIndex=1, 
                        xlim=xlim, ylim=ylim, legpos="topleft")
par(op)

# plot a series of maps at specific times
timeTxt    <- c("2022-01-01","2022-01-15",
                "2022-02-01","2022-02-15",
                "2022-03-01","2022-03-15",
                "2022-04-01","2022-04-15",
                "2022-05-01","2022-05-15")
timepoints <- apply(as.matrix(timeTxt), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE)
tpts       <- apply(as.matrix(timepoints), 1, pts_and_hpds_at_time, tr=tr,
                    xlim=xlim, ylim=ylim)
for (j in 1:length(tpts)) {
  op <- par(mar=c(0,0,0,0))
  plot_at_time(tpt=tpts[[j]], timePt=timepoints[j], tr=tr, xlim=xlim, ylim=ylim)
  par(op)
}



