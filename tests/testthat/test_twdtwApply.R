# test_twdtwApply.R

#require(testthat)
#require(dtwSat)

# get the pattern's data
pname.list <- list()
pdata.list <- list()
pindex.list <- list()
for(i in 1:length(patterns.list)){
  pname.list[[i]] <- names(patterns.list)[[i]]
  pdata.list[[i]] <- coredata(patterns.list[[i]])
  pindex.list[[i]] <- index(patterns.list[[i]])
}

patterns.list.Pct <- list()
patterns.list.Plt <- list()
for(i in 1:length(patterns.list)){
  pind <- as.POSIXct(pindex.list[[i]])
  patterns.list.Pct[[i]] <- zoo(x = pdata.list[[i]], order.by = pind)
  pind <- as.POSIXlt(pindex.list[[i]])
  patterns.list.Plt[[i]] <- zoo(x = pdata.list[[i]], order.by = pind)
}

patterns_ts = twdtwTimeSeries(patterns.list)
patterns_ts.Pct = twdtwTimeSeries(patterns.list.Pct)
patterns_ts.Plt = twdtwTimeSeries(patterns.list.Plt)

expect_true(is(index(patterns_ts[[1]]), "Date"))
expect_true(is(index(patterns_ts.Pct[[1]]), "Date"))
expect_true(is(index(patterns_ts.Plt[[1]]), "Date"))
