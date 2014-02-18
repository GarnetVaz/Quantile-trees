#!/usr/bin/env R
## Time-stamp: <2014-02-17 17:26:27 garnet>
## Simple code to test working of qtree.
## For crossvalidation results check the other tests.

library('qtree',quietly=TRUE)

## Download datasets if not exists.
source('./downloadData.R')

load('./qtree_datasets/wine-quality-white.RData')
load('./qtree_datasets/wine-quality-red.RData')

form <- as.formula(quality ~ .)

## Note that the meaning of minsize and mincut is flipped from the 'tree' package.
mod.white.qtree <- qtree(form, data = whitewine.data, mindev = 0.001, minsize = 2, mincut = 4, tau = 0.5)
mod.red.qtree <- qtree(form, data = redwine.data, mindev = 0.001, minsize = 2, mincut = 4, tau = 0.5)

mod.white.tree <- tree(form, data = white, mindev = 0.001, minsize = 4, mincut = 2)
mod.red.tree <- tree(form, data = red, mindev = 0.001, minsize = 4, mincut = 2)
