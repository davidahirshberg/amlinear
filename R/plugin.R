#' Estimate plugin balancing weights
#'
#' @param X the input features
#' @param W the effect variable (real valued)
#' @param w.hat optional precomupted estimate of E[W | X]
#' @param alpha tuning paramter for glmnet
#' @param standardize whether X should be standardized for the lasso
#' 
#' @return balancing weights
#' 
#' @export balance_plugin
balance_plugin = function(X, W, w.hat = NULL, alpha = 1, standardize = TRUE) {
    
    if (is.null(w.hat)) {
        w.fit = crossfit.cv.glmnet(X, W, alpha = alpha)
        w.hat = crossfit.predict(w.fit)
    }
    
    var.fit = crossfit.cv.glmnet(X, (W - w.hat)^2, alpha = alpha, standardize = standardize)
    var.hat = crossfit.predict(var.fit)
    
    gamma.raw = (W - w.hat) / var.hat
    gamma.centered = gamma.raw - mean(gamma.raw)
    gamma.scaled = gamma.centered / sqrt(mean(gamma.centered^2) * mean((W - w.hat)^2))
    gamma.scaled
}