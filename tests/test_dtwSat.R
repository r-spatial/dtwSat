
library(dtwSat)
library(ggplot2)
library(reshape2)

# Class names
names(query.list)

# Perform twdtw
alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
              alpha = 0.1, beta = 50, alignments=4, keep=TRUE)
is(alig, "dtwSat")

# # 
# # # Print twdtw object
alig

# 
# # Plot twdtw object
# Plot path
gp1 = plot(alig, type="path", normalize=TRUE, show.dist = TRUE)
gp1
# Plot alignment
gp2 = plot(alig, type="alignment", dimension="evi", alignment=1, shift=0.5)
gp2
# ggsave("soya_paths.png", plot=gp1, width = 8.9, height=5.9, units="in", 
#         family="Helvetica")
# ggsave("soya_aligs.png", plot=gp2, width = 8.9, height=5.9, units="in", 
#         family="Helvetica")


#
# # Wavelet filter
sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
                      boundary = "periodic")
## Plot raw EVI and filtered EVI
df = data.frame(Time=index(template), value=template$evi, variable="Raw")
df = rbind( df, data.frame(Time=index(sy), value=sy$evi, variable="Wavelet filter") )
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() +
  ylab("EVI")

# ## Plot all filter bands
df = melt(data.frame(Time=index(sy), sy), id="Time")
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() 


# Perform twdtw to a query list 
malig = mtwdtw(query.list, template, weight = "logistic", 
               alpha = 0.1, beta = 50)
class(malig)
dim(malig)
malig



# # ## Other plot examples 
# # library(grid)
# # library(gridExtra)
# # gp1 = plotAlignment(alig, dimension="evi", alignment=1, shift=0.5)
# # gp2 = plotAlignment(alig, dimension="evi", alignment=2, shift=0.5)
# # grid.arrange(arrangeGrob(gp1 + ggtitle("Alignment 1") + 
# #                                         theme(axis.title.x=element_blank(),
# #                                               legend.position="none"),
# #                          gp2 + ggtitle("Alignment 2") + 
# #                                         theme(axis.title.x=element_blank(), 
# #                                               legend.position="none"),
# #             nrow=2))
# # 
# # ## Crate a list of dtwSat objects and plot them 
# # gp.list = lapply(query.list, function(query){
# #      alig = dtwSat(query, template, weight = "logistic", 
# #                     alpha = 0.1, beta = 50, alignments = 4, keep = TRUE)
# #      plotPath(alig, normalize = TRUE, show.dist = TRUE)  
# # })
# # grid.arrange(arrangeGrob(gp.list[[1]] + ggtitle(names(query.list)[1]) + 
# #                                                 theme(axis.title.x=element_blank(), 
# #                                                       legend.position="none"),
# #                          gp.list[[2]] + ggtitle(names(query.list)[2]) + 
# #                                                 theme(axis.title.x=element_blank(),
# #                                                       legend.position="none"),
# #                          gp.list[[3]] + ggtitle(names(query.list)[3]) + 
# #                                                       theme(legend.position="none"),
# #              nrow=3))
# 
