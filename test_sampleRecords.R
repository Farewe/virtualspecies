#################################
## Test sampleRecords()
## 
## TODO:  - set bias.strength argument in call to sampleRecords()
##        - test presence-only in call to sampleRecords()
##        - figure out why (in some previous script) the weights in my spat_nrec_hec 
##          datasets sum to much higher than 1.
## 
## author: Willson Gaul
## created 17 Nov 2017
## last modified: 8 Dec 2017
#################################

# working directory should already be set by being in the virtualspecies project

library(wgutil)
library(raster)
library(rgdal)
library(dismo)
library(ade4)
library(rworldmap)
library(virtualspecies)

load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/annual_precip_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/ETOPO1/elevation_hec_ETOPO1.RData")
load("~/Documents/Data_Analysis/UCD/sampling_distribution/spat_nrec_hec.RData")
## ----------------------- load coastline for masking -------------------------
# Load Ireland coastline
ir <- readOGR(dsn='../mapping/data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))
rm(ir)
## end data loading ----------------------------------------------------------

# turn spatialGridDataFrames into RasterStacks as required by virtualspecies
# colnames(krg_elev_predict@data)[1] <- "elevation"
elev_rst <- mask(elev_hec, ir_TM75)
mean_rr_rst <- mask(raster(krg_mean_rr_predict), ir_TM75)
var_stack <- stack(elev_rst, mean_rr_rst)

# generate environmental suitability
sp <- generateSpFromPCA(var_stack, 
                        means = c(2, 1), 
                        sds = c(2, 2))

plotResponse(sp)

# convert sp suitability to p/a
set.seed(200)
sp_pa <- convertToPA(sp, beta = 0.6, alpha = -0.1)
alt_sp_pa <- convertToPA(sp, beta = NULL, alpha = NULL)

# distribution bias 
# (not implementing right now)

# sample
# currently testing presence-absence.
# TODO: test presence-only
# To sample with spatial bias, use e.g. raster of moth sampling weights
# TODO set bias.strength
set.seed(351)
samp_1 <- sampleRecords(sp_pa, n = 100, 
                        type = "presence-absence", 
                        detection.probability = 0.7, 
                        bias = "manual", 
                        weights = spat_nrec_hec$`insect - moth`$rast)

samp_2 <- sampleRecords(sp_pa, n = 1000, 
                        type = "presence-absence", 
                        detection.probability = 0.7, 
                        bias = "manual", 
                        weights = spat_nrec_hec$bird$rast)

# set.seed(seed = attr(old_samp_3, "seed"), 
#          kind = attr(old_samp_3, "RNGkind"))
.Random.seed <- attr(old_samp_3, "seed")
samp_3 <- sampleRecords(sp_pa, n = 100, 
                        type = "presence only", 
                        detection.probability = 0.7, 
                        bias = "manual", 
                        weights = spat_nrec_hec$bird$rast)


################### checking --------------------------------------------------
library(tidyverse)
reduced_records_1 <- samp_1$sample.records %>%
  group_by(x, y) %>%
  summarise(sum(Observed))

reduced_records_2 <- samp_2$sample.records %>%
  group_by(x, y) %>%
  summarise(n_recs = sum(Observed))

reduced_records_3 <- samp_3$sample.records %>%
  group_by(x, y) %>%
  summarise(n_recs = sum(Observed))

ggplot(data = reduced_records_3, aes(x = x, y = y, size = n_recs)) + 
  geom_point()


##
old_samp_3 <- samp_3
new_samp_3 <- samp_3
reproduce_attempt_1 <- samp_3

old_samp_3$plots
new_samp_3$plots
samp_3$plots
reproduce_attempt_1$plots

identical(old_samp_3, new_samp_3)
identical(reproduce_attempt_1, samp_3)
