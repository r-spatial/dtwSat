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


# Show pattern names
names(patterns.list)

# 
plotPatterns(patterns.list)

# Perform twdtw alignment
zero_weight_fun = function(x) 0 
matches = twdtw(x=example_ts, patterns=patterns.list["Soybean"], step.matrix = symmetric1, 
             dist.method = "Euclidean", weight.fun = zero_weight_fun, keep=TRUE)

is(matches, "twdtw")
show(matches)
print(matches)
summary(matches)

# Plot cost matrix paths
gp1 = plot(x=matches, type="paths")
gp1


# Plot alignment
gp2 = plot(matches, type="alignments", attr="evi")
gp2


# Wavelet filter
sy = waveletSmoothing(x = example_ts, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")


# Plot raw EVI and filtered EVI
gp = plotTimeSeries(sy)
gp

# Plot all filtered bands
evi = merge(Raw=zoo(example_ts$evi), Wavelet=zoo(sy$evi))
evi = na.approx(evi)
gp = plotTimeSeries(evi)
gp

# Perform twdtw to patterns list 
log_fun = logisticWeight(alpha=-0.1, beta=100)
matches = twdtw(x=example_ts, patterns=patterns.list, 
              weight.fun = log_fun, keep=TRUE)

# Classify interval
best_class = classifyIntervals(x=matches, from=as.Date("2009-09-01"), 
                              to=as.Date("2013-09-01"), by = "6 month",
                              overlap=.3, threshold=Inf)
getAlignments(best_class)

gp = plotClassification(x=best_class)
gp

# Perform twdtw to patterns list 
matches = twdtw(x=example_ts.list[[2]], patterns = patterns.list, 
             weight.fun = log_fun, keep=TRUE)

# Classify interval
best_class = classifyIntervals(x=matches, from=as.Date("2007-09-01"), 
                               to=as.Date("2013-09-01"), by = "6 month",
                               overlap=.4, threshold=Inf)
getAlignments(best_class)

gp = plotClassification(x=best_class)
gp


# Apply twdtw to a list of time series 
matches = lapply(example_ts.list, twdtw, patterns = patterns.list, 
             weight.fun = log_fun, keep=TRUE)

gp.list = lapply(matches, plotClassification, from=as.Date("2007-09-01"),  
               to=as.Date("2013-09-01"), by = "6 month", overlap=.4)

grid.arrange(grobs=gp.list, ncol=1)


# Plot cost matrix 
matches = twdtw(example_ts.list[[1]], patterns=patterns.list[1], 
             weight.fun = log_fun, keep=TRUE)
 
gp = plotCostMatrix(x=matches, matrix.name="timeWeight")
gp

gp = plotCostMatrix(x=matches, matrix.name="localMatrix")
gp
 
gp = plotCostMatrix(x=matches, matrix.name="costMatrix")
gp

 
# Test one band  
patterns.list2 = lapply(patterns.list, function(qq) qq[,c("evi"), drop=FALSE])
matches = twdtw( x = example_ts.list[[2]], patterns = patterns.list2, 
              weight.fun = log_fun, keep=TRUE)
# Classify interval
best_class = classifyIntervals(x=matches, from=as.Date("2007-09-01"), 
                               to=as.Date("2013-09-01"), by = "6 month",
                               overlap=.4, threshold=Inf)
getAlignments(best_class)

gp = plotClassification(x=matches, from=as.Date("2007-09-01"),  
               to=as.Date("2013-09-01"), by = "6 month",
               overlap=.4)
gp




