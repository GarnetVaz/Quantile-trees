#!/usr/bin/env R
## Code downloads datasets from the web for testing purpose.
## Datasets include :
## Wine quality
## Communities and crime
## Concrete strength (TODO)

## Set download directory name
cur.dir <- getwd()
dir.create(file.path(cur.dir,'qtree_datasets'),showWarnings=FALSE)
data.dir <- paste(cur.dir,'/qtree_datasets',sep='')

if(!file.exists(data.dir)) {
    ## Download the wine dataset from http://www3.dsi.uminho.pt/pcortez/wine/
    wine.url <- 'http://www3.dsi.uminho.pt/pcortez/wine/winequality.zip'
    download.file(url=wine.url,destfile='./winequality.zip')
    unzip('./winequality.zip')
    whitewine.loc <- './winequality/winequality-white.csv'
    whitewine.save.loc <- './qtree_datasets/wine-quality-white.RData'
    redwine.loc <- './winequality/winequality-red.csv'
    redwine.save.loc <- './qtree_datasets/wine-quality-red.RData'
    whitewine.data <- read.csv(whitewine.loc,header=TRUE,sep=';')
    save(whitewine.data,file=whitewine.save.loc)
    redwine.data <- read.csv(redwine.loc,header=TRUE,sep=';')
    save(redwine.data,file=redwine.save.loc)

    ## Download the communities and crime dataset from UCI database.
    crime.url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.data'
    download.file(url=crime.url,destfile='./crime.data')
    crime.alldata <- read.table('./crime.data',sep=',',na.strings='?')
    ## Eliminate columns with NA
    na.cols <- (sapply(crime.alldata, function(x) sum(is.na(x))) > 0)
    na.cols[1:4] <- TRUE                    # Categorical variables
    crime.data <- crime.alldata[,!na.cols]
    names(crime.data) <- paste('X',seq(1L:dim(crime.data)[2]),sep="")
    crime.save.loc <- './qtree_datasets/crime.RData'
    save(crime.data,file=crime.save.loc)

    ## Clean up
    wine.files <- list.files('./winequality')
    for (i in wine.files) {
        fname <- paste('./winequality',i,sep="/")
        file.remove(fname)
    }
    file.remove('./winequality')
    file.remove('./winequality.zip')
    file.remove('./crime.data')
}
