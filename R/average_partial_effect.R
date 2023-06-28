#' Estimate average partial effect via augmented balancing
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param balance.method how the balancing weights gamma are derived
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fitted.model optional pre-fitted regression adjustment
#' @param alpha tuning paramter for glmnet
#' @param standardize whether non-binary features should be noramlized
#' @param solver convex optimzer used by CVXR for minimax weights
#' @param verbose whether the optimizer should print progress information
#'
#' @return ATE estimate with standard error estimate. Also returns ``linear''
#'         point estimate of the form sum gamma_i Yi, as in Donoho (1994), for comparison.
#'
#' @export average_partial_effect
average_partial_effect = function(X, Y, W,
                                  balance.method = c("minimax", "plugin"),
                                  zeta=0.5,
                                  fitted.model = NULL,
                                  alpha=1,
                                  standardize = TRUE,
                                  solver = c("ECOS", "SCS"),
                                  verbose = TRUE) {
    solver = match.arg(solver) 
    balance.method = match.arg(balance.method)
    
    if (standardize) {
        scl = apply(X, 2, sd, na.rm = TRUE)
        is.binary = apply(X, 2, function(xx) sum(xx == 0) + sum(xx == 1) == length(xx))
        scl[is.binary] = 1
        X.scl = scale(X, center = FALSE, scale = scl)
    } else {
        X.scl = X
    }
    
    # Compute regression adjustment
    if (is.null(fitted.model)) {
        lasso.out = rlasso(X.scl, Y, W, alpha, standardize = FALSE)
    } else {
        lasso.out = fitted.model
    }
    tau.hat = lasso.out$tau.hat
    w.hat = lasso.out$w.hat
    
    # Compute balancing weights
    if (balance.method == "minimax") {
        gamma = balance_minimax(X.scl, W, zeta, solver = solver, verbose = verbose)
    } else if (balance.method == "plugin") {
        gamma = balance_plugin(X.scl, W, w.hat, alpha, standardize = FALSE)
    } else {
        stop("Unrecognized balancing method.")
    }
    
    # Compute point estimate and standard errors
    m.hat = lasso.out$y.hat + (W - lasso.out$w.hat) * lasso.out$tau.hat
    point.estimate = mean(tau.hat + gamma * (Y - m.hat))
    Vhat = mean((tau.hat - point.estimate)^2 + gamma^2 * (Y - m.hat)^2)
    standard.error.estimate = sqrt(Vhat / length(W))
    ret = c(point.estimate=point.estimate,
            standard.error.estimate=standard.error.estimate,
            linear.point.estimate=mean(gamma * Y))
    
    return(ret)
}
