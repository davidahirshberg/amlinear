set.seed(1)

test_that("cross-fitting removes optimism", {
    
    p = 403
    n = 404
    beta = 1/(1:p)
    sigma = 0.2
    
    X = matrix(rnorm(n*p), n, p)
    Y = X %*% beta + sigma * rnorm(n)
    
    X.test = matrix(rnorm(n*p), n, p)
    Y.test = X.test %*% beta + sigma * rnorm(n)
    
    lasso = glmnet::cv.glmnet(X, Y)
    Y.hat = predict(lasso, X)
    Y.test.hat = predict(lasso, X.test)
    
    lasso.train = mean((Y - Y.hat)^2)
    lasso.test = mean((Y.test - Y.test.hat)^2)
    
    classo = crossfit.cv.glmnet(X, Y)
    cY.hat = crossfit.predict(classo)
    cY.test.hat = predict(classo, X.test)
    
    classo.train = mean((Y - cY.hat)^2)
    classo.test = mean((Y.test - cY.test.hat)^2)
    
    expect_equal(classo.test / lasso.test, 1, tol = 0.15)
    expect_equal(classo.train / lasso.test, 1, tol = 0.15)
    
    expect_lt(lasso.train / lasso.test, 0.7)
    
})

test_that("rlasso is accurate", {
    
    p = 603
    n = 304
    
    beta.e = 1/(p:1)
    beta.m = 1/(4 + p:1)
    sigma.Y = 4
    sigma.W = 2
    beta.tau = rep(0, p+1)
    beta.tau[1:5] = c(1:5)/5
    
    X = matrix(rnorm(n*p), n, p)
    TAU = cbind(1, X) %*% beta.tau
    W = X %*% beta.e + sigma.W * rnorm(n)
    Y = X %*% beta.m + W * TAU + sigma.Y * rnorm(n)
    
    rr = rlasso(X, Y, W)
    rerror = mean((rr$tau.hat - TAU)^2)
    rerror.coef = sum((beta.tau - rr$tau.beta)^2)
    
    expect_lt(rerror.coef / sum(beta.tau^2), 0.5)
    
    #
    # For quasi-orthogonal designs, the joint lasso is a strong baseline
    #
    
    jj = cv.glmnet(cbind(W, X, as.numeric(W) * X), Y, penalty.factor = c(0, rep(1, 2 * p)))
    jj.beta = coef(jj)
    jj.beta.tau = jj.beta[c(2, 2 + p + 1:p)]
    jerror = mean((cbind(1, X) %*% (jj.beta.tau - beta.tau))^2)
    jerror.coef = sum((beta.tau - jj.beta.tau)^2)
    
    expect_lt(rerror.coef / jerror.coef, 2)
    expect_lt(rerror / jerror, 2)
    
})
