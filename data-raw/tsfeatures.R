library(tidyverse)
library(tsfeatures)

f1 <- read_csv("~/timeseriesData/comp-engine-export/comp-engine-export-metadata.20181003.csv")
f2 <- read_csv("~/timeseriesData/comp-engine-export/comp-engine-export-datapoints.20181003.csv")

birdsongs <- filter(f1, category=="Birdsong")
airtemp <- filter(f1, category=="Air Temperature")
music <- filter(f1, category=="Music")

d_birdsongs <- filter(f2, timeseries_id %in% birdsongs$timeseries_id)$datapoints
d_airtemp <- filter(f2, timeseries_id %in% airtemp$timeseries_id)$datapoints
d_music <- filter(f2, timeseries_id %in% music$timeseries_id)$datapoints

birdsongs <- lapply(strsplit(d_birdsongs, ","), as.numeric)
airtemp <- lapply(strsplit(d_airtemp, ","), as.numeric)
music <- lapply(strsplit(d_music, ","), as.numeric)

feat_birdsongs <- tsfeatures(birdsongs)
feat_airtemp <- tsfeatures(airtemp)
feat_music <- tsfeatures(music)

fb <- add_column(feat_birdsongs, type="birdsongs")
fa <- add_column(feat_airtemp, type="air temperature")
fm <- add_column(feat_music, type="music")

tsfeatureData <- rbind(fb, fa, fm) %>%
  select(-frequency, -nperiods, -seasonal_period)

use_data(tsfeatureData)
