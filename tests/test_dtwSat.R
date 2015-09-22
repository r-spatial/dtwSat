# Rscript for open boundary DTW time series analysis
# Author: Victor Maus
# Dec. 2014

library(dtwSat)

# Class names
names(query.list)

# Perform twdtw
alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4, keep=TRUE)
is(alig, "dtwSat")

# # 
# # # Print twdtw object
alig

# 
# # Plot twdtw object
gp = plot(alig, normalize=TRUE, show.dist = TRUE)
gp
# ggsave("alig_soy.png", plot=gp, width = 8.9, height=5.9, units="in", family="Helvetica")

#
# # Wavelet time series smoothing
sy = timeSeriesSmoothing(template$evi, frequency=16, method=c("wavelet",1))
df.y = melt(data.frame(Time=index(template), Raw=template$evi), id="Time")
df.sy = melt(data.frame(Time=index(sy), Wavelet=sy), id="Time")
df = rbind(df.y, df.sy)
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() + 
  ylab("EVI")



# Perform twdtw to a query list 
malig = mtwdtw(query.list, template, weight = "logistic", alpha = 0.1, beta = 50)
class(malig)
dim(malig)
malig

# Plot alignment for all classese
gp.list = lapply(query.list, function(query){
  alig = twdtw(query, template, weight = "logistic", alpha = 0.1, beta = 50, alignments = 4, keep = TRUE)
  plot(alig, normalize = TRUE, show.dist = TRUE)  
})
grid.arrange(arrangeGrob(gp.list[[1]] + ggtitle(names(query.list)[1]) + theme(axis.title.x=element_blank(), legend.position="none"),
                         gp.list[[2]] + ggtitle(names(query.list)[2]) + theme(axis.title.x=element_blank(), legend.position="none"),
                         gp.list[[3]] + ggtitle(names(query.list)[3]) + theme(legend.position="none"),
                         nrow=3))


