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
gp1 = plot(alig, type="alignment", dimension="evi", alignment=1, shift=0.5)
gp2 = plot(alig, type="alignment", dimension="evi", alignment=2, shift=0.5)
grid.arrange(arrangeGrob(gp1 + ggtitle("Alignment 1") + theme(axis.title.x=element_blank(), legend.position="none"),
                         gp2 + ggtitle("Alignment 2") + theme(axis.title.x=element_blank(), legend.position="none"),
                         nrow=2))

# Plot alignment for all classese
gp.list = lapply(query.list, function(query){
  alig = twdtw(query, template, weight = "logistic", alpha = 0.1, beta = 50, alignments = 4, keep = TRUE)
  plot(alig, normalize = TRUE, show.dist = TRUE)  
})
grid.arrange(arrangeGrob(gp.list[[1]] + ggtitle(names(query.list)[1]) + theme(axis.title.x=element_blank()),
                         gp.list[[2]] + ggtitle(names(query.list)[2]) + theme(axis.title.x=element_blank()),
                         gp.list[[3]] + ggtitle(names(query.list)[3]) ,
                         nrow=3))

#
# # Wavelet filter
sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, boundary = "periodic")

# Plot raw EVI and filtered EVI
df.y = melt(data.frame(Time=index(template), Raw=template$evi), id="Time")
df.sy = melt(data.frame(Time=index(sy), Wavelet=sy$evi), id="Time")
df = rbind(df.y, df.sy)
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() + 
  ylab("EVI")

# Plot all filter bands
df = melt(data.frame(Time=index(sy), sy), id="Time")
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() 


# Perform twdtw to a query list 
malig = mtwdtw(query.list, template, weight = "logistic", alpha = 0.1, beta = 50)
class(malig)
dim(malig)
malig




