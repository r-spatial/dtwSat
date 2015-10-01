---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

dtwSat
=====

### Time-Weighted Dynamic Time Warping for remote sensing time series analysis
The dtwSat provides a Time-Weighted Dynamic Time Warping (TWDTW) algorithm to measure similarity between two temporal sequences. This adaptation of the classical Dynamic Time Warping (DTW) algorithm is flexible to compare events that have a strong time dependency, such as phenological stages of cropland systems and tropical forests. This package provides methods for visualization of minimum cost paths, time series alignment, and time intervals classification.

### Install

```r
devtools::install_github("vwmaus/dtwSat")
```


### Quick demo

This dome performs a dtwSat analysis and show the results.

```r
library(dtwSat)
```

```r
names(query.list)
```

```
## [1] "Soybean" "Cotton"  "Maize"
```

```r
alig = twdtw(query.list[["Soybean"]], template, 
             weight = "logistic", alpha = 0.1, beta = 50, alignments=4, keep=TRUE) 
print(alig)
```

```
## Time-Weighted DTW alignment object
## Number of alignments: 4 
##   query       from         to distance normalizedDistance
## 1     1 2011-10-04 2012-01-28 3.956483         0.03140066
## 2     1 2012-10-06 2013-02-15 4.008838         0.03181617
## 3     1 2009-09-13 2010-03-05 4.539202         0.03602541
## 4     1 2010-10-20 2011-03-18 5.528445         0.04387655
```

### Plot examples

Load libraries

```r
library(dtwSat)
library(ggplot2)
library(reshape2)
library(gridExtra)
```

Plot alignments

```r
gp1 = plot(alig, type="alignment", attribute="evi", alignment=1, shift=0.5) + 
          ggtitle("Alignment 1") +
		      theme(axis.title.x=element_blank())
```

```
## Error in eval(expr, envir, enclos): could not find function "ggtitle"
```

```r
gp2 = plot(alig, type="alignment", attribute="evi", alignment=2, shift=0.5) +
          ggtitle("Alignment 2") + 
          theme(legend.position="none")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggtitle"
```

```r
grid.arrange(gp1,gp2,nrow=2)
```

```
## Error in eval(expr, envir, enclos): could not find function "grid.arrange"
```


Plot path for all classese

```r
gp.list = lapply(query.list, function(query){
  				alig = twdtw(query, template, weight = "logistic", alpha = 0.1, 
  				             beta = 100, alignments = 4, keep = TRUE)
  				plot(alig, normalize = TRUE, show.dist = TRUE) + 
  				  theme(axis.title.x=element_blank(),
  				        legend.position="none")
})
```

```
## Error in FUN(X[[i]], ...): could not find function "theme"
```

```r
grid.arrange(gp.list[[1]] + ggtitle(names(query.list)[1]),
             gp.list[[2]] + ggtitle(names(query.list)[2]),
             gp.list[[3]] + ggtitle(names(query.list)[3]),
             nrow=3)
```

```
## Error in eval(expr, envir, enclos): could not find function "grid.arrange"
```

Plot classification

```r
malig = mtwdtw(query.list, template, weight = "logistic", 
        alpha = 0.1, beta = 100)
 
gp1 = plot(x=malig, type="classify", from=as.Date("2009-09-01"),  
     to=as.Date("2013-09-01"), by = "6 month",
     normalized=TRUE, overlap=.7) 

gp2 = plot(x=malig, type="classify", attribute = c("evi","ndvi"),
           from=as.Date("2009-09-01"), to=as.Date("2013-09-01"), 
           by = "6 month", normalized=TRUE, overlap=.7)

grid.arrange(gp1,gp2,nrow=2)
```

```
## Error in eval(expr, envir, enclos): could not find function "grid.arrange"
```



Plot wavelet smoothing

```r
sy = waveletSmoothing(x=template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")
df1 = data.frame(Time=index(template), value=template$evi, variable="Raw")
df1 = rbind(df1, data.frame(Time=index(sy), value=sy$evi, variable="Wavelet filter") )
gp1 = ggplot(df1, aes(x=Time, y=value, group=variable, colour=variable)) +
   		geom_line() + 
  		# theme(legend.position="bottom") +
  		ylab("EVI")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggplot"
```

```r
df2 = melt(data.frame(Time=index(sy), sy), id="Time")
```

```
## Error in eval(expr, envir, enclos): could not find function "melt"
```

```r
gp2 = ggplot(df2, aes(x=Time, y=value, group=variable, colour=variable)) +
   		geom_line() + 
  		ylab("Value") 
```

```
## Error in eval(expr, envir, enclos): could not find function "ggplot"
```

```r
grid.arrange(gp1,gp2,nrow=2)
```

```
## Error in eval(expr, envir, enclos): could not find function "grid.arrange"
```




