#' Estimate treatment effect via the R-lasso, as proposed by Nie and Wager (2017)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validated

rlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se")) {
    
    lambda.choice = match.arg(lambda.choice)
    
    
    nobs = nrow(X)
    pobs = ncol(X)
    
    if (is.null(nfolds)) {
        nfolds = floor(max(3, min(10, sum(W==0)/5, sum(W==1)/5)))
    }
    
    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(nfolds), length = length(W)))
    
    y.fit = crossfit.cv.glmnet(X, Y, foldid = foldid, lambda.choice = lambda.choice, alpha = alpha)
    y.hat = crossfit.predict(X, y.fit, foldid)
    
    w.fit = crossfit.cv.glmnet(X, W, foldid = foldid, lambda.choice = lambda.choice, alpha = alpha)
    w.hat = crossfit.predict(X, w.fit, foldid)
    
    Y.tilde = Y - y.hat
    X.tilde = (W - w.hat) * X
    
    tau.fit = crossfit.cv.glmnet(X.tilde, Y.tilde, foldid = foldid, lambda.choice = lambda.choice, alpha = alpha)
    tau.hat = crossfit.predict(X, tau.fit, foldid)
    
    return(data.frame(tau.hat = tau.hat, y.hat = y.hat, w.hat = w.hat))
}

crossfit.predict = function(X, lasso.obj, foldid) {
    sapply(1:length(foldid), function(i) {
        sum(lasso.obj$cv.betas[foldid[i],] * c(1, X[i,1]))
    })
}
