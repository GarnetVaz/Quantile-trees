#!/usr/bin/env R
## Small test case
set.seed(1)

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

## noarma output:
## 1) root 5 9.1090 -0.6212
##   2) x1 < 3.17400735590602 2 0.4335 -1.0550 *
##   3) x1 > 3.17400735590602 3 7.1700  2.3900
##     6) x1 < 8.46239428697658 2 5.3680  7.7570 *
##     7) x1 > 8.46239428697658 1 0.0000 -1.2150 *
