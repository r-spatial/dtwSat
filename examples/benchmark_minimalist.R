library(dtwSat)
log_fun = logisticWeight(-0.1, 50)
from = "2009-09-01"
to = "2017-09-01"
by = "6 month"

# S4 objects for original implementation 
tw_patt = twdtwTimeSeries(MOD13Q1.patterns.list)
tw_ts = twdtwTimeSeries(MOD13Q1.ts) 

# Table from csv for minimalist version 
mn_patt <- lapply(dir(system.file("reduce_time/patterns", package = "dtwSat"), full.names = TRUE), read.csv, stringsAsFactors = FALSE)
mn_ts <- read.csv(system.file("reduce_time/ts_MODIS13Q1.csv", package = "dtwSat"), stringsAsFactors = FALSE)

# Benchtmark 
rbenchmark::benchmark(
  original = twdtwClassify(twdtwApply(x = tw_ts, y = tw_patt, weight.fun = log_fun), from = from, to = to, by = by)[[1]],
  minimalist = .twdtw_reduce_time(x = mn_ts, y = mn_patt, weight.fun = log_fun, from = from, to = to, by = by)  
)

