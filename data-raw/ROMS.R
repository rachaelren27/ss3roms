# Provided by Kristin Marshall via email on 2021-12-07
# Future versions of the data will have more columns and more rows
# because a Dynamic Factor Analysis will condense the variables into
# a single variable and the time series will extend up to the current year.
# All downstream code will be able to accommodate new data sets.
ROMSdata <- utils::read.csv(file.path("data-raw", "Data_ROMS.predictors.csv"))
ROMSdata[is.na(ROMSdata)] <- -999
usethis::use_data(ROMSdata, overwrite = TRUE)
