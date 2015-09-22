# Rscript for open boundary DTW time series analysis
# Author: Victor Maus
# Dec. 2014

library(dtwSat)

# Class names
names(query.list)

# Perform twdtw
alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4, keep=TRUE)
is(alig, "dtwSat")

# # 
# # # Print twdtw object
alig

# 
# # Plot twdtw object
plot(alig, normalize=TRUE, show.dist = TRUE)

#
# # Wavelet time series smoothing
sy = timeSeriesSmoothing(template$evi, frequency=16, method=c("wavelet",1))
df.y = melt(data.frame(Time=index(template), Raw=template$evi), id="Time")
df.sy = melt(data.frame(Time=index(sy), Wavelet=sy), id="Time")
df = rbind(df.y, df.sy)
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() + 
  ylab("EVI")


# 
# # Perform twdtw to a query list 
# alig = mtwdtw(query.list, template, weight = "logistic", alpha = 0.1, beta = 50)
# alig
# 
# # Plot alignment for all classese
# par(mfrow=c(3,1))
# query.names = names(query.list)
# for(p in query.names){
#   alig = twdtw(query.list[[p]], template, weight = "logistic", alpha = 0.1, beta = 50)
#   plot(alig, ylab=p, show.dist = TRUE)
# }
