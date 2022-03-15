
twdtw.as.data.frame <- function(x){
  as.data.frame(date = index(x), x[,,drop = FALSE])
}

twdtw.as.data.frame.list <- function(x){
  lapply(x, twdtw.as.data.frame)
}

twdtw.as.zoo <- function(x){
  zoo(x[, !names(x) %in% "date", drop = FALSE], order.by = x$date)
}

twdtw.as.zoo.list <- function(x){
 lapply(x, twdtw.as.zoo)
}

data.frame.as.twdtw <- function(x){
  
}

list.as.twdtw <- function(x){
  zoo(x[, !names(x) %in% "date", drop = FALSE], order.by = x$date)
}
