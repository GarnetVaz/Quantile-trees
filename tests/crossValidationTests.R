#!/usr/bin/env R
## Run cross-validation on LAD vs OLS trees

rm(list=ls())
library('caret')
library('tree')
library('qtree')

## Load datasets
stopifnot(file.exists('./qtree_datasets'))

load('./qtree_datasets/wine-quality-white.RData')
load('./qtree_datasets/wine-quality-red.RData')
load('./qtree_datasets/crime.RData')

datasets <- list('whitewine' = whitewine.data, 'redwine' = redwine.data, 'crime' = crime.data)

computeLAD <- function(ytrue, yhat) {
    N <- length(ytrue)
    return(sum(abs(ytrue-yhat))/N)
}

computeMSD <- function(ytrue, yhat) {
    N <- length(ytrue)
    return(sum((ytrue-yhat)^2)/N)
}

crossValidate <- function(dat, formula, folds) {
    errorTreeLAD <- vector(mode="numeric", length=length(folds))
    errorTreeMSD <- vector(mode="numeric", length=length(folds))
    errorQtreeLAD <- vector(mode="numeric", length=length(folds))
    errorQtreeMSD <- vector(mode="numeric", length=length(folds))
    for(i in 1L:length(folds)) {
        train <- dat[-folds[[i]],]
        test <- dat[folds[[i]],]
        ## Build qtree model
        mod.qtree <- qtree(formula, data = train,
                           mindev = 0.01, mincut = 10, minsize = 5,tau=0.5)
        mod.predict.qtree <- predict(mod.qtree, newdata = test)
        errorQtreeLAD[i] <- computeLAD(test[,dim(test)[2]], mod.predict.qtree)
        errorQtreeMSD[i] <- computeMSD(test[,dim(test)[2]], mod.predict.qtree)
        ## Build tree model
        mod.tree <- tree(formula, data = train,
                         mindev = 0.01, mincut = 5, minsize = 10)
        mod.predict.tree <- predict(mod.tree, newdata = test)
        errorTreeLAD[i] <- computeLAD(test[,dim(test)[2]], mod.predict.tree)
        errorTreeMSD[i] <- computeMSD(test[,dim(test)[2]], mod.predict.tree)
    }
    return(list(treeLAD = errorTreeLAD,
                treeMSD = errorTreeMSD,
                qtreeLAD = errorQtreeLAD,
                qtreeMSD = errorQtreeMSD))
}

for(k in 1L:length(datasets)) {
    cat("Starting simulation for", names(datasets)[k],"\n")
    allData <- datasets[[k]]
    oldNames <- names(allData)
    oldNames[length(oldNames)] <- "y"
    names(allData) <- oldNames
    numRuns <- 10
    numFolds <- 10
    form <- as.formula(y ~.)
    result <- list()
    for(i in 1L:numRuns) {
        folds <- createFolds(1L:dim(allData)[1],numFolds)
        result[[i]] <- crossValidate(allData, form, folds)
    }

    cat("Mean of tree LAD is", mean(sapply(result, function(x) mean(x$treeLAD))),"\n")
    cat("Mean of tree MSD is", mean(sapply(result, function(x) mean(x$treeMSD))),"\n")
    cat("Mean of qtree LAD is", mean(sapply(result, function(x) mean(x$qtreeLAD))),"\n")
    cat("Mean of qtree MSD is", mean(sapply(result, function(x) mean(x$qtreeMSD))),"\n")

    cat("SD of tree LAD is", mean(sapply(result, function(x) sd(x$treeLAD))),"\n")
    cat("SD of tree MSD is", mean(sapply(result, function(x) sd(x$treeMSD))),"\n")
    cat("SD of qtree LAD is", mean(sapply(result, function(x) sd(x$qtreeLAD))),"\n")
    cat("SD of qtree MSD is", mean(sapply(result, function(x) sd(x$qtreeMSD))),"\n")
    cat("Done with simulation for", names(datasets)[k],"\n\n")
}
