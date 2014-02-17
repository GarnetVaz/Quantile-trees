#!/usr/bin/env R

library('qtree',quietly=TRUE)

## Download the wine dataset from http://www3.dsi.uminho.pt/pcortez/wine/
url <- 'http://www3.dsi.uminho.pt/pcortez/wine/winequality.zip'
download.file(url=url,destfile='winequality.zip')
unzip('winequality.zip')

white <- read.csv('./winequality/winequality-white.csv',sep=';',header=TRUE)
red <- read.csv('./winequality/winequality-red.csv',sep=';',header=TRUE)
form <- as.formula(quality ~ .)

mod.white.qtree <- qtree(form, data = white, mindev = 0.001, minsize = 2, mincut = 4, tau = 0.5)
mod.red.qtree <- qtree(form, data = red, mindev = 0.001, minsize = 2, mincut = 4, tau = 0.5)

mod.white.tree <- tree(form, data = white, mindev = 0.001, minsize = 4, mincut = 2)
mod.red.tree <- tree(form, data = red, mindev = 0.001, minsize = 4, mincut = 2)
