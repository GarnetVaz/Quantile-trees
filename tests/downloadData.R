#!/usr/bin/env R
## Code downloads datasets from the web for testing purpose.
## Datasets include :
## Wine quality
## Communities and crime
## Cement quality

## Set download directory name
cur.dir <- getwd()
dir.create(file.path(cur.dir,'qtree_datasets'),showWarnings=False)
data.dir <- paste(cur.dir,'qtree_datasets',sep='')

if(!file.exists(data.dir)) {
    ## Download the wine dataset from http://www3.dsi.uminho.pt/pcortez/wine/
    wine.url <- url <- 'http://www3.dsi.uminho.pt/pcortez/wine/winequality.zip'
    download.file(url=url,destfile='./qtree_datasets/winequality.zip')
    unzip('qtree_datasets/winequality.zip')
    whitewine.loc <- './qtree_datasets/winequality/winequality-white.csv'
    redwine.loc <- './qtree_datasets/winequality/winequality-red.csv'
    whitewine.data <- read.csv(whitewine.loc,header=)

    ## Download the communities and crime dataset from UCI database.
    crime.url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.data'
    download.file(url=crime.url,destfile='./qtree_datasets/crime.data')
    crime.alldata <- read.table('./qtree_datasets/crime.data',sep=',',na.strings='?')
    ## Eliminate columns with NA
    na.cols <- (sapply(crime.alldata, function(x) sum(is.na(x))) > 0)
    na.cols[1:4] <- TRUE                    # Categorical variables

    crime.subset <- crime.alldata[,!na.cols]
    names(crime.subset) <- paste('X',seq(1L:dim(crime.subset)[2]),sep="")
    save.loc
    save(crime.subset,file='')
}
