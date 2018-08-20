rm(list = ls())

library(glmnet)
library(CVXR)


balance_minimax = function(M,
                           balance.target,
                           zeta = 1/2,
                           solver = c("ECOS", "SCS"),
                           verbose = TRUE) {
    solver = match.arg(solver)
    nobs = nrow(M)
    pobs = ncol(M)
    gg = CVXR::Variable(nobs + 1)
    objective = (1 - zeta) * sum(gg[1:nobs]^2) + zeta * sum(gg[nobs + 1]^2)
    contraints = list(
        t(M) %*% gg[1:nobs] - balance.target <= gg[nobs + 1],
        -t(M) %*% gg[1:nobs] + balance.target <= gg[nobs + 1]
    )
    cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
    cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
    result = cvx.output$getValue(gg)
    gamma = nobs * result[1:nobs]
    gamma
}



#n = 400
#p = 800
#confound = TRUE

focal.idx = 1

nrep = 40
nn = c(200, 400, 800, 1600)
conf = c(FALSE, TRUE)

results = Reduce(rbind, lapply(nn, function(n) {
    p = 2 * n
    Reduce(rbind, lapply(conf, function(confound) {
        Reduce(rbind, lapply(1:nrep, function(rrr) {
            X = matrix(rnorm(n * p), n, p)
            
            if (confound) {
                X[,focal.idx] = sqrt(10 / 3) * rowMeans(X[, 11:20]) + sqrt(2 / 3) * rnorm(n)
            }
            
            theta = 20 /((1:p) + 5)^2
            theta[focal.idx] = -0.1
            prob = 1/(1 + exp(-X %*% theta))
            Y = rbinom(n, 1, prob)
            
            ape = theta[focal.idx] * mean(prob * (1 - prob))
            
            lasso = cv.glmnet(X, Y, family = "binomial",
                              keep = TRUE)
            
            intercept.hat = coef(lasso)[1]
            theta.hat = coef(lasso)[-1]
            theta.hat.focal = theta.hat[focal.idx]
            
            # extract out-of-fold (= prevalidated) predictions from the lasso
            prob.hat.lasso = lasso$fit.preval[,!is.na(colSums(lasso$fit.preval))][, lasso$lambda == lasso$lambda.1se]
            ape.hat.plugin = theta.hat.focal * mean(prob.hat.lasso * (1 - prob.hat.lasso))
            
            X.plus = cbind(1, X)
            psi = prob.hat.lasso * (1 - prob.hat.lasso)
            psi.prime = prob.hat.lasso * (1 - prob.hat.lasso) * (1 - 2 * prob.hat.lasso)
            balance.target = theta.hat.focal * t(X.plus) %*% psi.prime / n + 
                mean(psi) * as.numeric(0:p == focal.idx)
            M = as.numeric(psi) * X.plus
            
            gamma.hat = balance_minimax(M, balance.target, verbose = FALSE)
            
            ape.hat.debias = mean(gamma.hat * (Y - prob.hat.lasso))
            ape.hat = ape.hat.plugin + ape.hat.debias
            
            A.hat = sort(unique(c(focal.idx, which(theta.hat > 0))))
            logreg.wz = glm(Y ~ X[,A.hat], family = binomial)
            prob.hat.wz = predict(logreg.wz, type = "response")
            tilde.theta.focal = coef(logreg.wz)[-1][which(A.hat == focal.idx)]
            ape.hat.wz = as.numeric(tilde.theta.focal * mean(prob.hat.wz * (1 - prob.hat.wz)))
            
            out = c(confound=confound, n=n, p=p, oracle=ape, plugin=ape.hat.plugin, WZ=ape.hat.wz, DEB=ape.hat)
            print(out)
            out
        }))
    }))
}))

colMeans(results - mean(results[,"oracle"]))
sqrt(colMeans((results - mean(results[,"oracle"]))^2))

write.csv(results, "results.csv")
