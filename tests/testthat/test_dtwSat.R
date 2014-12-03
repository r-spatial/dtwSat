# Rscript for open boundary DTW time series analysis
# Author: Victor Maus
# Dec. 2014



# Set workspace and processing parameters
setwd("~/workspace/AllBlue/test")
URLSERVER = "http://www.dpi.inpe.br/mds/mds"
DATASETS = "evi2"
LIBRARYPATH = "../data/Tlibrary"
coords = data.frame(longitude=-52.415303449,latitude=-12.3755423768)
# coords = data.frame(longitude=-52.4549374155,latitude=-12.3384079189)
# coords = data.frame(longitude=-52.40429,latitude=-12.35496)


#######################################################################################
#
### Get patterns/signatures/queries whatever you want to call it
TemporalPaternsFiles = dir(LIBRARYPATH,pattern=paste(DATASETS,".csv",sep=""));
TemporalPatterns.list = lapply(paste(LIBRARYPATH,TemporalPaternsFiles,sep="/"),
                               read.table, as.is=TRUE, header=FALSE, sep=",")
classes = unlist(strsplit(TemporalPaternsFiles, split="\\.csv"))
classes = unlist(unique(lapply(strsplit(classes, split="\\."), function(name) name[1])))
names(TemporalPatterns.list) = classes


# for(g in G){
    # Step 1. Get time series 
    timeseries = getTimeSeries(wtssClient(URLSERVER), coverages="MOD13Q1", 
                               datasets=c(unlist(DATASETS), "day2"), 
                               latitude=coords[1,2], longitude=coords[1,1], 
                               from="2000-01-01", to="2014-12-31")
    timeline = as.Date(timeseries[[1]]$datasets$timeline)
    ty = as.Date(timeseries[[1]]$datasets$day2)
    y = as.numeric(timeseries[[1]]$datasets[,DATASETS])
    
    # Step 2. Pre-processing: Compute regular time series and smoothing
    template = timeSeriesSmoothing(ty, y, timeline, method=c("wavelet",1))
    gp = ggplot(data = data.frame(x=ty, y=y), aes( x = as.Date(x), y = y )) + 
        ylim(c(0,1)) + 
        geom_line(size = .5, colour="black", fill="black") + 
        geom_line(size = .8, data=data.frame(x=as.Date(index(template)), y=as.numeric(template)), color="green") + 
        scale_x_date(labels = date_format("%Y"), breaks = date_breaks("year")) +
        xlab("Time") + 
        ylab("EVI2") +
        theme(text=element_text(size = 10, family="Helvetica"), plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 10, family="Helvetica"), axis.text.y = element_text(size = 10, family="Helvetica"))
    print(gp)
    
    # Step 3. Processing: Compute dtw and other metrics
    results = lapply(seq_along(TemporalPatterns.list), function(j){
      # j = 1
      patternName=names(TemporalPatterns.list)[j]
      x = as.numeric(TemporalPatterns.list[[j]][,2])
      tx = as.Date(TemporalPatterns.list[[j]][,1], origin="1970-01-01")
      query = zoo(x, tx)
      return(timeSeriesAnalysis(query, template, theta=0.8, span=0.5, satStat=FALSE))
    })
    names(results) = names(TemporalPatterns.list)
    
    # Step 4. Post-processing
    finalClassification = timeSeriesClassifier(results, from="2000-09-01", to="2014-08-31", 
                         by = 12, overlapping = 0.0, threshold=4.0)

# }

View(finalClassification)




