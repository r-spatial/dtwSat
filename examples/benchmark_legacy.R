remotes::install_github("vwmaus/dtwSat")

library(dtwSat)
alpha = -0.1
beta = 50
from = "2009-09-01"
to = "2017-08-31"
by = "12 month"

# S4 objects for original implementation 
tw_patt = readRDS(system.file("lucc_MT/patterns/patt.rds", package = "dtwSat"))
tw_ts = twdtwTimeSeries(MOD13Q1.ts) 

# Table from csv for legacy version 
mn_patt <- lapply(dir(system.file("lucc_MT/patterns", package = "dtwSat"), pattern = ".csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE)
mn_ts <- read.csv(system.file("reduce_time/ts_MODIS13Q1.csv", package = "dtwSat"), stringsAsFactors = FALSE)

# Benchtmark
rbenchmark::benchmark(
  t1_s4_legacy  = twdtwClassify(x = tw_ts, y = tw_patt, weight.fun = logisticWeight(alpha, beta), from = from, to = to, by = by, legacy = TRUE),
  t2_s4_fast    = twdtwClassify(x = tw_ts, y = tw_patt, from = from, to = to, by = by, alpha = alpha, beta = beta, legacy = FALSE),
  t3_s4_fast_tw = twdtwClassify(x = tw_ts, y = tw_patt, from = from, to = to, by = by, alpha = alpha, beta = beta, legacy = FALSE),
  t4_s3_fast    = twdtwClassify(x = mn_ts, y = mn_patt, from = from, to = to, by = by, alpha = alpha, beta = beta, time.window = FALSE),
  t5_s3_fast_tw = twdtwClassify(x = mn_ts, y = mn_patt, from = from, to = to, by = by, alpha = alpha, beta = beta, time.window = TRUE)
)

# TODO: increase TS density... 
library(profvis)
profvis({
  replicate(500, twdtwApply(x = mn_ts, y = mn_patt[1], from = from, to = to, by = by, alpha = alpha, beta = beta, time.window = TRUE))
})

twdtwClassify(x = mn_ts, y = mn_patt[1], from = from, to = to, by = by, alpha = alpha, beta = beta, time.window = FALSE)

plotClassification(twdtwClassify(x = tw_ts, y = tw_patt, from = from, to = to, by = by, alpha = alpha, beta = beta, legacy = FALSE))
