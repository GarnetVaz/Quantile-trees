#!/usr/bin/env R
## Time-stamp: <2014-02-19 14:01:19 garnet>
## Simple code to test working of qtree.
## For crossvalidation results check the other tests.
library('tree',quietly=TRUE)
library('qtree',quietly=TRUE)

## Download datasets if not exists.
source('./downloadData.R')
rm(list=ls())
load('./qtree_datasets/wine-quality-white.RData')
load('./qtree_datasets/wine-quality-red.RData')

form <- as.formula(quality ~ .)

## Note that the meaning of minsize and mincut is flipped from the 'tree' package.
mod.white.qtree <- qtree(form, data = whitewine.data,
                         mindev = 0.1, minsize = 2, mincut = 4, tau = 0.5)
mod.red.qtree <- qtree(form, data = redwine.data,
                       mindev = 0.1, minsize = 2, mincut = 4, tau = 0.5)

mod.white.tree <- tree(form, data = whitewine.data,
                       mindev = 0.1, minsize = 4, mincut = 2)
mod.red.tree <- tree(form, data = redwine.data,
                     mindev = 0.1, minsize = 4, mincut = 2)

lad <- function(y, yhat) {
    return(sum(abs(y-yhat))/length(yhat))
}
## To predict for in sample : mod.red.qtree$frame$yval[mod.red.qtree$where]

## Simple in sample testing to ensure correctness while modifying code.
insample.white.qtree <- lad(predict(mod.white.qtree,newdata=whitewine.data),whitewine.data$quality)
insample.red.qtree <- lad(predict(mod.red.qtree,newdata=redwine.data),redwine.data$quality)

cat('In sample white wine MAD is ', insample.white.qtree, '\n')
cat('In sample red wine MAD is ', insample.red.qtree, '\n')

## Output: Master branch
## In sample white wine MAD is  0.4897918
## In sample red wine MAD is  0.3692933
