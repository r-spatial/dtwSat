---
output:
  md_document:
    variant: markdown_github
---

<!-- 
# Edit and run 
library(knitr)
knit(input="README.Rmd", output = "README.md")
-->

dtwSat
=====

### Time-Weighted Dynamic Time Warping for remote sensing time series analysis
dtwSat package provides a Time-Weighted Dynamic Time Warping (TWDTW) algorithm to measure similarity between two temporal sequences. This adaptation of the classical Dynamic Time Warping (DTW) algorithm is flexible to compare events that have a strong time dependency, such as phenological stages of cropland systems and tropical forests. This package provides methods for visualization of minimum cost paths, time series alignment, and time intervals classification.

### Install

```r
devtools::install_github("vwmaus/dtwSat")
```


### Quick demo

This dome performs a dtwSat analysis and show the results.

```r
library(dtwSat, quietly = TRUE)
names(query.list)
```

```
## [1] "Soybean" "Cotton"  "Maize"
```

```r
alig = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
             weight = "logistic", alpha = 0.1, beta = 50, span=180, keep=TRUE) 
is(alig, "dtwSat")
```

```
## [1] TRUE
```

```r
print(alig)
```

```
## Time-Weighted DTW alignment object
## Number of alignments: 4 
##   query       from         to distance
## 1     1 2011-10-27 2012-02-07 2.580247
## 2     1 2009-10-01 2010-03-05 2.795275
## 3     1 2010-11-03 2011-01-31 2.850574
## 4     1 2012-10-06 2013-02-15 3.056848
```

### Plot examples

Plot alignments

```r
library(dtwSat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
alig = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
             weight = "logistic", alpha = 0.1, beta = 50, span=180, keep=TRUE) 
gp1 = plot(alig, type="alignment", n=1, attr="evi", shift=0.5) + 
          ggtitle("Alignment 1") +
		      theme(axis.title.x=element_blank())
gp2 = plot(alig, type="alignment", n=2, attr="evi", shift=0.5) +
          ggtitle("Alignment 2") + 
          theme(legend.position="none")
grid.arrange(gp1,gp2,nrow=2)
```

![plot of chunk define-demo-plot-alignments](figure/define-demo-plot-alignments-1.png) 


Plot path for DTW and TWDTW

```r
library(dtwSat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
alig1 = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
              weight = NULL, span=180, keep = TRUE)
gp1 = plot(alig1, show.dist = TRUE) + 
  				  theme(axis.title.x=element_blank(), legend.position="none")
alig2 = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
              weight = "linear", theta = 1, span=180, keep = TRUE)
gp2 = plot(alig2, show.dist = TRUE) + 
  				  theme(axis.title.x=element_blank(), legend.position="none")
alig3 = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
              weight = "logistic", alpha = 0.1, beta = 100, span=180, keep = TRUE)
gp3 = plot(alig3, show.dist = TRUE) + 
  				  theme(axis.title.x=element_blank(), legend.position="none")

grid.arrange(gp1 + ggtitle("DTW"),
             gp2 + ggtitle("TWDTW - Linear weight"),
             gp3 + ggtitle("TWDTW - Nonlinear weight"),
             nrow=3)
```

![plot of chunk define-demo-twdtw-x-dtw](figure/define-demo-twdtw-x-dtw-1.png) 

Plot path for all classese

```r
library(dtwSat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
gp.list = lapply(seq_along(query.list), function(i){
  				alig = twdtw(query.list[[i]], waveletSmoothing(template), 
  				             weight = "logistic", alpha = 0.1, 
  				             beta = 100, span=180, keep = TRUE)
  				plot(alig, show.dist = TRUE) + 
  				                      ggtitle(names(query.list)[i]) + 
  				                      theme(axis.title.x=element_blank(), 
  				                      legend.position="none")
})
grid.arrange(grobs=gp.list, nrow=3)
```

![plot of chunk define-demo-plot-paths](figure/define-demo-plot-paths-1.png) 

Plot classification

```r
library(dtwSat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
malig = mtwdtw(query.list, waveletSmoothing(template), 
               normalize = TRUE, query.length = 23,
               weight = "logistic", alpha = 0.1, beta = 100)
 
gp1 = plot(x=malig, type="classify", from=as.Date("2009-09-01"),  
     to=as.Date("2013-09-01"), by = "6 month", overlap=.3) 

gp2 = plot(x=malig, type="classify", attr = c("evi","ndvi"),
           from=as.Date("2009-09-01"), to=as.Date("2013-09-01"), 
           by = "6 month", overlap=.3)

grid.arrange(gp1,gp2,nrow=2)
```

![plot of chunk define-demo-plot-classification](figure/define-demo-plot-classification-1.png) 



Plot wavelet smoothing

```r
library(dtwSat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
sy = waveletSmoothing(x=template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")

gp1 = autoplot(sy, facets = NULL) + xlab("Time")

evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))

gp2 = autoplot(evi, facets = NULL) + xlab("Time")

grid.arrange(gp1,gp2,nrow=2)
```

![plot of chunk define-demo-plot-smoothing](figure/define-demo-plot-smoothing-1.png) 




