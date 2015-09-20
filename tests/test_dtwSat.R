# Rscript for open boundary DTW time series analysis
# Author: Victor Maus
# Dec. 2014

library(dtwSat)

# Class names
names(query.list)

# Perform twdtw
alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)

# Print twdtw object
print(alig)

# Plot twdtw object
plot(alig, main="Soybean", show.dist = TRUE)

# Wavelet time series smoothing
sy = timeSeriesSmoothing(template$evi, frequency=16, method=c("wavelet",1))
plot(template$evi, xlab="Time", ylab="EVI")
lines(sy, col="red")

# Perform twdtw to a query list 
alig = mtwdtw(query.list, template, weight = "logistic", alpha = 0.1, beta = 50)
alig

# Plot alignment for all classese
par(mfrow=c(3,1))
query.names = names(query.list)
for(p in query.names){
  alig = twdtw(query.list[[p]], template, weight = "logistic", alpha = 0.1, beta = 50)
  plot(alig, ylab=p, show.dist = TRUE)
}
