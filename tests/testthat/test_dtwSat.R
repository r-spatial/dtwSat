###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### dtwSat TESTS


# Show query names
names(query.list)


# Perform twdtw alignment
alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
              alpha = 0.1, beta = 50, alignments=4, keep=TRUE)


# Test if object is dtwSat
is(alig, "dtwSat")


# Show dtwSat object
show(alig)
print(alig)


# Plot cost matrix paths
gp1 = plot(x=alig, type="path", normalize=TRUE, show.dist=TRUE)
gp1


# Plot alignment
gp2 = plot(alig, type="alignment", attributeibute="evi", alignment=1, shift=0.5)
gp2


# Wavelet filter
sy = waveletSmoothing(x=template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")


# Plot raw EVI and filtered EVI

df = data.frame(Time=index(template), value=template$evi, variable="Raw")
df = rbind( df, data.frame(Time=index(sy), value=sy$evi, variable="Wavelet filter") )
gp = ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() + 
  theme(legend.position="bottom") +
  ylab("EVI")
gp
# ggsave("filter.png", plot=gp, width = 8.9, height=5.9/1.5, units="in",
#         family="Helvetica")

# Plot all filtered bands
df = melt(data.frame(Time=index(sy), sy), id="Time")
ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
  geom_line() 


# Perform twdtw to query list 
malig = mtwdtw(query.list, template, weight = "logistic", 
               alpha = 0.1, beta = 100)
class(malig)
getAlignments(malig)
      
# Classify interval
best_class = classfyIntervals(x=malig, from=as.Date("2009-09-01"), 
                              to=as.Date("2013-09-01"), by = "6 month",
                              normalized=TRUE, overlap=.7, threshold=Inf)
best_class


malig = mtwdtw(query.list, template, weight = "logistic", 
               alpha = 0.1, beta = 100)
 
gp = plotClassify(x=malig, from=as.Date("2009-09-01"),  
              to=as.Date("2013-09-01"), by = "6 month",
              normalized=TRUE, overlap=.7) 
gp
# ggsave("classify.png", plot=gp, width = 8.9, height=5.9/1.5, units="in",
#         family="Helvetica", type = "cairo-png")

