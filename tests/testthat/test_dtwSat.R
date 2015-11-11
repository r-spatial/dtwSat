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
alig = twdtw(query.list[["Soybean"]], timeseries=template, dist.method = "Euclidean",
             weight = "logistic", alpha = 0.1, beta = 100, span=180, keep=TRUE)

is(alig, "dtwSat")
show(alig)
print(alig)


# Plot cost matrix paths
gp1 = plot(x=alig, type="path", show.dist=TRUE)
gp1


# Plot alignment
gp2 = plot(alig, type="alignment", n=1, attr="evi", shift=0.5)
gp2


# Wavelet filter
sy = waveletSmoothing(timeseries = template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")


# Plot raw EVI and filtered EVI
gp = autoplot(sy, facets = NULL) + xlab("Time")
gp

# Plot all filtered bands
evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
gp = autoplot(evi, facets = NULL) + xlab("Time")
gp

# Normalize queries length
new.query.list = normalizeQuery(query = query.list, query.length = 23)
data.frame(Old.Length=unlist(lapply(query.list, nrow)), New.length=unlist(lapply(new.query.list, nrow)))

# Perform twdtw to query list 
malig = mtwdtw(query.list, timeseries=template, weight = "logistic", 
               alpha = 0.1, beta = 100, normalize=TRUE, 
               query.length = 23, keep=TRUE, fast=TRUE)
# Classify interval
best_class = classifyIntervals(x=malig, from=as.Date("2009-09-01"), 
                              to=as.Date("2013-09-01"), by = "6 month",
                              overlap=.3, threshold=Inf)
best_class

gp = plotClassify(x=malig, from=as.Date("2009-09-01"),  
              to=as.Date("2013-09-01"), by = "6 month",
              overlap=.3)
gp
# ggsave("classify.png", plot=gp, width = 8.9, height=5.9/1.5, units="in",
#         family="Helvetica", type = "cairo-png")



# Perform twdtw to query list 
malig = mtwdtw(query.list, template.list[[2]], weight = "logistic", 
               alpha = 0.1, beta = 100, normalize=TRUE, query.length = 23)

# Classify interval
best_class = classifyIntervals(x=malig, from=as.Date("2007-09-01"), 
                               to=as.Date("2013-09-01"), by = "6 month",
                               overlap=.3, threshold=Inf)
best_class

gp = plotClassify(x=malig, from=as.Date("2007-09-01"),  
                  to=as.Date("2013-09-01"), by = "6 month",
                  overlap=.3)
gp


# Plot cost matrix 
alig = twdtw(query.list[["Soybean"]], timeseries=template, weight = "logistic", 
             alpha = 0.1, beta = 100, keep=TRUE)
 

gp = plotCostMatrix(x=alig, matrix.name="timeWeight")
gp

gp = plotCostMatrix(x=alig, matrix.name="localMatrix")
gp
 
gp = plotCostMatrix(x=alig, matrix.name="costMatrix")
gp


# Test bands order 
query.list2 = lapply(query.list, function(qq) qq[,c("evi"),drop=FALSE])
malig = mtwdtw(query = query.list2, timeseries = template.list[[2]], weight = "logistic", 
               alpha = 0.1, beta = 100, normalize=TRUE, query.length = 23, keep=TRUE)
# Classify interval
best_class = classifyIntervals(x=malig, from=as.Date("2007-09-01"), 
                               to=as.Date("2013-09-01"), by = "6 month",
                               overlap=.3, threshold=Inf)
best_class




