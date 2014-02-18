#!/usr/bin/env R
## Small test case
set.seed(1)

library('tree')
library('qtree')
## Simple dataset
x1 <- c(2,4,3,7,8)
x1 <- x1 + rnorm(5)
x2 <- c(12,19,17,3,5)
x2 <- x2 + rnorm(5)
y <- c(-3, 2, 0, 1, 12)
y <- y + rnorm(5)


data.simple <- data.frame(x1,x2,y)

mod.data <- qtree(y ~ ., data = data.simple, tau = 0.5, minsize = 1, mincut = 2)
