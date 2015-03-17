# Rscript for open boundary DTW time series analysis
# Author: Victor Maus
# Dec. 2014



# Set workspace and processing parameters
setwd("~/workspace/dtwSat/tests")
# LIBRARYPATH = "../../data/Tlibrary/embrapa/Tlibrary-"
LIBRARYPATH = "../../data/Tlibrary/embrapa/Tlibrary++"
# LIBRARYPATH = "../../data/Tlibrary/embrapa/amazon_simulation_all"
URLSERVER = "http://www.dpi.inpe.br/mds/mds"
COVERAGES = "MOD13Q1"
DATASETS = "evi"
FROM="2000-01-01"
TO="2014-08-31"
# COORDINATES = data.frame(longitude=-52.415303449,latitude=-12.3755423768)
# COORDINATES = data.frame(longitude=-52.4549374155,latitude=-12.3384079189)
COORDINATES = data.frame(longitude=-52.40429,latitude=-12.35496)
# COORDINATES = data.frame(longitude=-56.725483,latitude=-12.217708)
# COORDINATES = data.frame(longitude=-56.7065251011,latitude=-12.2115047254)
# COORDINATES = data.frame(longitude=-56.8204205399,latitude=-12.2681083867)
# COORDINATES = data.frame(longitude=-56.8058059977,latitude=-12.2547399261)
# COORDINATES = data.frame(longitude=-56.5851917,latitude=-12.259192635)

#######################################################################################
#
### Get patterns/signatures/queries whatever you want to call it
TemporalPaternsFiles = dir(LIBRARYPATH,pattern=paste(DATASETS,".csv",sep=""))
TemporalPatterns.list = lapply(paste(LIBRARYPATH,TemporalPaternsFiles,sep="/"),
                               read.table, as.is=TRUE, header=FALSE, sep=",")
classes = unlist(strsplit(TemporalPaternsFiles, split="\\.csv"))
# classes = unlist(unique(lapply(strsplit(classes, split="\\."), function(name) name[1])))
classes = unlist(lapply(strsplit(classes, split="\\."), function(name) name[1]))
names(TemporalPatterns.list) = classes


# Step 1. Download time series
datasets = DATASETS
if(COVERAGES=="MOD13Q1")
  datasets = c(DATASETS,"day2")
timeseries = getTimeSeries(wtssClient(URLSERVER), coverages=COVERAGES, datasets=datasets,
                           longitude=COORDINATES[1], latitude=COORDINATES[2],
                           from=FROM, to=TO)

# Step 2. Pre-processing: Compute regular time series and smoothing
timeline = as.Date(timeseries[[1]]$datasets$timeline)
IDs = which(FROM<=timeline & timeline<=TO)
ty = timeline[IDs]
if(COVERAGES=="MOD13Q1")
  ty = as.Date(timeseries[[1]]$datasets$day2[IDs])
y = as.numeric(timeseries[[1]]$datasets[,DATASETS][IDs])
template = timeSeriesSmoothing(ty, y, timeline, method=c("wavelet",1))
gp = ggplot(data = data.frame(x=ty, y=y), aes( x = as.Date(x), y = y )) + 
  ylim(c(0,1)) + 
  geom_line(size = 0.8, colour="black", fill="black") + 
  geom_line(size = 1.0, data=data.frame(x=as.Date(index(template)), y=as.numeric(template)), color="blue") + 
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 year") ) +
  xlab("Time") + 
  ylab(toupper(DATASETS)) +
  theme(text=element_text(size = 10, family="Helvetica"), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, 
        family="Helvetica"), 
        axis.text.y = element_text(size = 10, family="Helvetica"),
        legend.position=c(0.8,.8)) + 
  guides(colour=guide_legend(title="Time series"))
print(gp)

# Step 3. Processing: Compute dtw and other metrics
results = lapply(seq_along(TemporalPatterns.list), function(j){
  # j = 1
  patternName=names(TemporalPatterns.list)[j]
  x = as.numeric(TemporalPatterns.list[[j]][,2])
  tx = as.Date(TemporalPatterns.list[[j]][,1], origin="1970-01-01")
  query = zoo(x, tx)
  out = timeSeriesAnalysis2(query, template, method="logistictimeweight", threshold=Inf, normalize=TRUE, 
                                           delay = NULL, theta=NULL, alpha=0.1, beta=100)
  return(out)
})
names(results) = names(TemporalPatterns.list)

# Step 4. Post-processing
finalClassification = timeSeriesClassifier(results, from=FROM, to=TO, 
                                           by = 12, overlapping = 0.15, threshold=1, 
                                           sortBydtw=TRUE, aggregateByClass=TRUE)
# results

View(data.frame(finalClassification))

