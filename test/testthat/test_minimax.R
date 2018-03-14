p = 203
n = 604

beta.e = 1/(p:1)
beta.m = 1/(4 + p:1)
sigma.Y = 2
sigma.W = 4
beta.tau = rep(0, p+1)
beta.tau[1:5] = c(5:1)/2

X = matrix(rnorm(n*p), n, p)
TAU = cbind(1, X) %*% beta.tau
W = X %*% beta.e + sigma.W * (0.3 + 1/ (1 + exp(X %*% beta.tau[-1]))) * rnorm(n)
Y = X %*% beta.m + W * TAU + sigma.Y * rnorm(n)

APE = average_partial_effect(X, Y, W, balance.method = "minimax", verbose = FALSE)

test_that("plugin estimator is accurate", {
    expect_equal(as.numeric(APE["point.estimate"]), mean(TAU), tolerance = 0.3)
    expect_lt(abs(as.numeric(APE["point.estimate"]) - mean(TAU)), 3 * as.numeric(APE["standard.error.estimate"]))
})


sY = Y + 1
APEsY = average_partial_effect(X, sY, W, balance.method = "minimax", verbose = FALSE)

test_that("plugin estimator is Y-traslantion invariant", {
    expect_equal(APEsY["point.estimate"], APE["point.estimate"], tolerance = 0.05)
    expect_equal(APEsY["standard.error.estimate"], APE["standard.error.estimate"], tolerance = 0.02)
})

sW = W + 1
APEsW = average_partial_effect(X, Y, sW, balance.method = "minimax", verbose = FALSE)

test_that("plugin estimator is W-traslantion invariant", {
    expect_equal(APEsW["point.estimate"], APE["point.estimate"], tolerance = 0.05)
    expect_equal(APEsW["standard.error.estimate"], APE["standard.error.estimate"], tolerance = 0.02)
})

aY = 2 * Y
aW = 2 * W
APEaWY = average_partial_effect(X, aY, aW, balance.method = "minimax", verbose = FALSE)

test_that("plugin estimator is scale invariant", {
    expect_equal(APEaWY["point.estimate"], APE["point.estimate"], tolerance = 0.08)
    expect_equal(APEaWY["standard.error.estimate"], APE["standard.error.estimate"], tolerance = 0.04)
})
