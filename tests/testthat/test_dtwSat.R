
# Set workspace
setwd("~/workspace/AllBlue/test")
# URLSERVER = "http://chronos.dpi.inpe.br:6543/mds" # Local access to Chronos webservice access
URLSERVER = "http://www.dpi.inpe.br/mds/mds" # External access to Chronos webservice access
DATASETS = "evi2"
LIBRARYPATH = "../data/Tlibrary"
coords = data.frame(longitude=-52.415303449,latitude=-12.3755423768)

### Get patterns/signatures
TemporalPaternsFiles = dir(LIBRARYPATH,pattern=paste(DATASETS,".csv",sep=""));
TemporalPatterns.list = lapply(paste(LIBRARYPATH,TemporalPaternsFiles,sep="/"), read.table, as.is=TRUE, header=FALSE, sep=",")
classes = unlist(strsplit(TemporalPaternsFiles, split="\\.csv"))
classes = unlist(unique(lapply(strsplit(classes, split="\\."), function(name) name[1])))
names(TemporalPatterns.list) = classes

### Get time series 
timeseries = getTimeSeries(wtssClient(URLSERVER), coverages="MOD13Q1", datasets=c(unlist(DATASETS), "day2"), latitude=coords[1,2], longitude=coords[1,1], from="2000-01-01", to="2014-12-31")
timeline = as.Date(timeseries[[1]]$datasets$timeline)
ty = as.Date(timeseries[[1]]$datasets$day2)
y = as.numeric(timeseries[[1]]$datasets[,DATASETS])

# Pre-proc: Compute regular time series and smoothing
template = timeSeriesSmoothing(ty, y, timeline, method=c("wavelet",1))
plot(ty, y, type="l")
lines(template, col="green")

results = lapply(seq_along(TemporalPatterns.list), function(j){
  # j = 1
  patternName=names(TemporalPatterns.list)[j]
  x = as.numeric(TemporalPatterns.list[[j]][,2])
  tx = as.Date(TemporalPatterns.list[[j]][,1], origin="1970-01-01")
  query = zoo(x, tx)
  return(timeSeriesAnalysis(query, template, theta=0.8, satStat=FALSE))
})
names(results) = names(TemporalPatterns.list)
results



