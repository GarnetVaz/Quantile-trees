require('tree')                         # ensures that printing and prediction works.

qtree <- function (formula,
                   data,
                   mindev=0.01,
                   mincut=5,
                   minsize=10,
                   tau=0.5,
                   na.action = na.pass,
                   model = FALSE,
                   x = FALSE,y = TRUE,
                   wts = TRUE,...)
{

	if (is.data.frame(model)){
		m <- model
		model <- FALSE
	}
	else{
		m <- match.call(expand.dots = FALSE)
		m$mindev <- m$mincut <- m$minsize <- m$tau <- m$model <- m$x <- m$y <- m$wts <- NULL
		m[[1L]] <- as.name("model.frame.default")
		m <- eval.parent(m)
	}

	CALL = match.call()
	Terms <- attr(m,"terms")
	if (any(attr(Terms, "order") > 1)){
        	stop("qtree cannot handle interaction terms")
	}
	if(any(lapply(m,class)=="factor"))
		stop("qtree cannot handle categorical input yet!")
	# extract response variable
    Y <- model.extract(m, "response")
    if (is.matrix(Y) && ncol(Y) > 1L)
        stop("trees cannot handle multiple responses")
    w <- model.extract(m, "weights")
    if (!length(w))
        w <- rep(1, nrow(m))
    if (any(yna <- is.na(Y))) {
        Y[yna] <- 1
        w[yna] <- 0
    }
    ## if (2*mincut < minsize) {
    ##     stop("minsize should be atleast 2*mincut")
    ## }
    offset <- attr(Terms, "offset")
    if (!is.null(offset)) {
        offset <- m[[offset]]
        Y <- Y - offset
    }

    X <- as.matrix(m[,-1])
	xlevels0 <- attr(X, "column.levels")
    if (is.null(xlevels0)) {
        xlevels0 <- rep(list(NULL), ncol(X))
        names(xlevels0) <- dimnames(X)[[2L]]
    }

	nobs <- length(Y)
    if (nobs == 0L)
        stop("no observations from which to fit a model")
	# define tuning parameters
        mylist <- qtreeCPP(X,Y,mindev,mincut,minsize,tau)
    ##     mylist = .Call("qtreeCPP", X, Y, mindev, mincut,
    ## minsize, tau, PACKAGE = "qtree")

	ourtree = with(mylist, {
        splits = NULL
        splits = rbind(splits, replicate(2, replicate(length(val),
            "")))
        colnames(splits) = c("cutleft", "cutright")
        indadd = which(valguide == 1)
        splits[indadd, 1] = paste("<", val[indadd], sep = "")
        splits[indadd, 2] = paste(">", val[indadd], sep = "")
        splits = as.array(splits)
        yval = as.numeric(round(yval, 5))
        dev = as.numeric(round(dev, 4))
		# columns names can be different from "X1"... "Xn"
        varind = (var!="<leaf>")
        varind2 = as.integer(substr(var[varind] ,2,10000L))
        var[varind] = colnames(X)[varind2]
        var = factor(var, levels = c("<leaf>", names(xlevels0)))
        mydataframe = data.frame(var = var, n = n, dev = dev, yval = yval)
        mydataframe$splits = splits
        rownames(mydataframe) = nodeID
        varleaf = which(var == "<leaf>")
        mywhere = integer()
        for (i in c(1:length(varleaf))) {
            indices = leaflist[[i]] + 1
            mywhere[indices] = varleaf[i]
        }
        browser()
        names(mywhere) = c(1:length(Y))
        otree = list(frame = mydataframe, where = mywhere, terms = Terms,
			call = CALL)
        attr(otree$where, "names") <- row.names(m)
	    if (length(n) > 1L)
    	    class(otree) <- "tree"
    	else class(otree) <- c("singlenode", "tree")
    	attr(otree, "xlevels") <- xlevels0
	    if (is.logical(model) && model)
    	    otree$model <- m
	    if (x)
    	    otree$x <- X
	    if (y)
    	    otree$y <- Y
	    if (wts)
    	    otree$weights <- w
		otree
	})
	return(ourtree)
}
