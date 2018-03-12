#' Estimate average partial effect via augmented balancing
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param balance.method how the balancing weights gamma are derived
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x, w) = E[Y | X = x, W = w]
#' @param alpha tuning paramter for glmnet
#' @param scale.X whether non-binary features should be noramlized
#' @param verbose whether the optimizer should print progress information
#'
#' @return ATE estimate, along with standard error estimate
#'
#' @export average_partial_effect
average_partial_effect = function(X, Y, W,
                                  balance.method = c("minimax", "plugin"),
                                  zeta=0.5,
                                  fit.method = c("elnet", "none"),
                                  alpha=1,
                                  scale.X = TRUE,
                                  verbose = FALSE) {
    
    balance.method = match.arg(balance.method)
    fit.method = match.arg(fit.method)
    
    if (fit.method == "none") {
        warning(paste("Setting fit.method to none is not recommended.",
                      "Resulting confidence intervals may not be valid."))
    }
    
    if (scale.X) {
        scl = apply(X, 2, sd, na.rm = TRUE)
        is.binary = apply(X, 2, function(xx) sum(xx == 0) + sum(xx == 1) == length(xx))
        scl[is.binary] = 1
        X.scl = scale(X, center = FALSE, scale = scl)
    } else {
        X.scl = X
    }
    
    # Compute regression adjustment
    if (fit.method = "none") {
        y.hat = rep(0, length(Y))
        tau.hat = rep(0, length(Y))
    } else if (fit.method = "elnet") {
        lasso.out = rlasso(X, Y, W, alpha)
    } else {
        stop("Unrecognized fitting method.")
    }
    
    # Compute balancing weights
    if (balance.method == "minimax") {
        gamma = balance_minimax(X, W, zeta)
    } else if (balance.method == "plugin") {
        gamma = balance_plugin(X, W, alpha)
    } else {
        stop("Unrecognized balancing method.")
    }
    
    # Compute point estimate and standard errors
    point.estimate = mean(tau.hat + gamma * (Y - y.hat))
    Vhat = mean((tau.hat - point.estimate)^2 + gamma^2 * (Y - y.hat)^2)
    standard.error.estimate = sqrt(Vhat / length(W))
    
    return(c(point.estimate=point.estimate,
             standard.error.estimate=standard.error.estimate))
}
