#' Estimate minimax balancing weights
#'
#' @param X the input features
#' @param W the effect variable (real valued)
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param solver convex optimzer used by CVXR
#' @param verbose whether convex optimizer should print verbose output
#' 
#' @return balancing weights
#' 
#' @export balance_minimax
balance_minimax = function(X, W, zeta, solver = c("ECOS", "SCS"), verbose = TRUE) {
    solver = match.arg(solver)
    nobs = nrow(X)
    pobs = ncol(X)
    gg = CVXR::Variable(nobs + 2)
    objective = (1 - zeta) * sum(gg[1:nobs]^2) + zeta * sum(gg[nobs + 1:2]^2)
    contraints = list(
        sum(gg[1:nobs]) == 0,
        t(X) %*% gg[1:nobs] <= gg[nobs + 1],
        -t(X) %*% gg[1:nobs] <= gg[nobs + 1],
        sum(W * gg[1:nobs]) == 1,
        t(X) %*% (W * gg[1:nobs]) <= colMeans(X) + gg[nobs + 2],
        - t(X) %*% (W * gg[1:nobs]) <= - colMeans(X) + gg[nobs + 2]
    )
    cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
    cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
    result = cvx.output$getValue(gg)
    gamma = nobs * result[1:nobs]
    gamma
}