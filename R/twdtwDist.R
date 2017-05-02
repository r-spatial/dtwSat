
.twdtwDist = function (x, y, ...) {
  d = do.call("rbind", lapply(x, function (xx) {
      sapply(y, function (yy) {
          res = twdtwApply(twdtwTimeSeries(xx), twdtwTimeSeries(yy), n=1, ...)
          res[[1]]$distance[1]
      })
  }))
  as.matrix(d)
}

