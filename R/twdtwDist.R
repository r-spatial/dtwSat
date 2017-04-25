
.twdtwDist = function (x, y, ...) {
  d = do.call("rbind", lapply(x, function (xx) {
      sapply(y, function (yy) {
          res = twdtwApply(twdtwTimeSeries(xx), twdtwTimeSeries(yy), n=1, ...)
          res[[1]]$distance[1]
      })
  }))
  as.matrix(d)
}

# library(dtwSat)
# library(dtw)
# library(dtwclust)
# log_fun = logisticWeight(-0.1, 100)
# phi = proxy::dist(x=MOD13Q1.patterns.list, y=MOD13Q1.patterns.list[1:2], method="TWDTW", weight.fun=log_fun)
# phi
# a=twdtwApply(twdtwTimeSeries(MOD13Q1.patterns.list), twdtwTimeSeries(MOD13Q1.patterns.list), n=1, weight.fun=log_fun)
# a[[3]]
# 
# aux = c(MOD13Q1.patterns.list, MOD13Q1.patterns.list, MOD13Q1.patterns.list)
# a1 = tsclust(series = aux, distance = "DTW", trace = TRUE, type = "hierarchical", k = 3)
# log_fun = logisticWeight(-0.1, 100)
# a2 = tsclust(series = aux, distance = "TWDTW", trace = TRUE, type = "hierarchical", k = 3, weight.fun=log_fun, span = 300)
# a1
# a2

