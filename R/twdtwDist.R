

.twdtwDist = function (object, ...) {
  d = lapply(seq_along(object), function (i) { 
    x = subset(object,  i)
    y = subset(object, -i)
    res = twdtwApply(x, y, n=1)
    res[]
  })
  d
}

