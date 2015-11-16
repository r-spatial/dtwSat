# ###############################################################
# #                                                             #
# #   (c) Victor Maus <vwmaus1@gmail.com>                       #
# #       Institute for Geoinformatics (IFGI)                   #
# #       University of Muenster (WWU), Germany                 #
# #                                                             #
# #       Earth System Science Center (CCST)                    #
# #       National Institute for Space Research (INPE), Brazil  #
# #                                                             #
# #                                                             #
# #   R Package dtwSat - 2015-09-01                             #
# #                                                             #
# ###############################################################
# 
# 
# ###############################################################
# #### dtwSat TESTS
# 
# 
# # Show pattern names
names(patterns.list)

# Perform twdtw alignment
weight.fun = function(x) 0.1*x
alig = twdtw(patterns=patterns.list["Soybean"], timeseries=template, step.matrix = symmetric1, 
             dist.method = "Euclidean", weight.fun = weight.fun, keep=TRUE)

is(alig, "dtwSat")
show(alig)
print(alig)


# Plot cost matrix paths
gp1 = plot(x=alig, type="path", show.dist=TRUE)
grid.arrange(gp1)


# Plot alignment
gp2 = plot(alig, type="alignment", attr="evi")
gp2


# Wavelet filter
sy = waveletSmoothing(timeseries = template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")


# Plot raw EVI and filtered EVI
gp = autoplot(sy, facets = NULL)
gp

# Plot all filtered bands
evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
gp = autoplot(evi, facets = NULL)
gp

# Normalize queries length
new.patterns.list = normalizePatterns(patterns = patterns.list, patterns.length = 23)
data.frame(Old.Length=sapply(patterns.list, nrow), New.length=sapply(new.patterns.list, nrow))

# Perform twdtw to patterns list 
weight.fun = logisticWeight(alpha=-0.1, beta=100)
malig = twdtw(patterns=new.patterns.list, timeseries=template,
              weight.fun = weight.fun, keep=TRUE)

# Classify interval
best_class = classifyIntervals(x=malig, from=as.Date("2009-09-01"), 
                              to=as.Date("2013-09-01"), by = "6 month",
                              overlap=.3, threshold=Inf)
best_class

gp = plotGroup(x=malig, from=as.Date("2009-09-01"),  
              to=as.Date("2013-09-01"), by = "6 month",
              overlap=.3)
gp
# # ggsave("classify.png", plot=gp, width = 8.9, height=5.9/1.5, units="in",
# #         family="Helvetica", type = "cairo-png")
# 
# 
# 
# # Perform twdtw to patterns list 
# malig = mtwdtw(patterns.list, timeseries = template.list[[2]], weight.fun = "logistic", 
#                alpha = 0.1, beta = 100, normalize.patterns=TRUE, patterns.length = 23)
# 
# # Classify interval
# best_class = classifyIntervals(x=malig, from=as.Date("2007-09-01"), 
#                                to=as.Date("2013-09-01"), by = "6 month",
#                                overlap=.3, threshold=Inf)
# best_class
# 
# gp = plotClassify(x=malig, from=as.Date("2007-09-01"),  
#                   to=as.Date("2013-09-01"), by = "6 month",
#                   overlap=.3)
# gp
# 
# 
# # Plot cost matrix 
# alig = twdtw(patterns.list[["Soybean"]], timeseries=template, weight.fun = "logistic", 
#              alpha = 0.1, beta = 100, keep=TRUE)
#  
# 
# gp = plotCostMatrix(x=alig, matrix.name="timeWeight")
# gp
# 
# gp = plotCostMatrix(x=alig, matrix.name="localMatrix")
# gp
#  
# gp = plotCostMatrix(x=alig, matrix.name="costMatrix")
# gp
# 
# 
# # Test bands order 
# patterns.list2 = lapply(patterns.list, function(qq) qq[,c("evi"),drop=FALSE])
# malig = mtwdtw(patterns = patterns.list2, timeseries = template.list[[2]], weight.fun = "logistic", 
#                alpha = 0.1, beta = 100, normalize.query=TRUE, query.length = 23, keep=TRUE)
# # Classify interval
# best_class = classifyIntervals(x=malig, from=as.Date("2007-09-01"), 
#                                to=as.Date("2013-09-01"), by = "6 month",
#                                overlap=.3, threshold=Inf)
# best_class
# 
# 
# 
# 
