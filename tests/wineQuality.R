#!/usr/bin/env R
library('qtree',quietly=TRUE)

## Download the wine dataset from http://www3.dsi.uminho.pt/pcortez/wine/
url <- 'http://www3.dsi.uminho.pt/pcortez/wine/winequality.zip'
download.file(url=url,destfile='winequality.zip')
unzip('winequality.zip')

white <- read.csv('./winequality/winequality-white.csv')
