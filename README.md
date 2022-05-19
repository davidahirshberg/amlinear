# amlinear: Augmented Minimax Linear Estimation

This package implements the augmented minimax linear estimator (AML) for estimating
average partial effects, as proposed by Hirshberg and Wager (2017). We consider a
setup where we observe data `(X, Y, W)` generated according
to a conditionally linear model,
```
Y = f(X) + W tau(X) + noise
```
and want to estimate the average effect `theta = E[tau(X)]` (both `W` and `Y` are taken
to be real valued). In the case where `W` is binary, then `theta` corresponds to the average
treatment effect under unconfoundedness, and our AML estimator is closely related to approximate 
residual balancing as implemented in <a href="https://github.com/swager/balanceHD">balanceHD</a>.

To install this package in R, run the following commands:
```R
library(devtools) 
install_github("davidahirshberg/amlinear")
```
Example usage:

```R
library(amlinear)

n = 300; p = 600; nclust = 10
beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))
clust.alpha = rep(c(0.3, 0.7), nclust/2)
cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
cluster = sample.int(nclust, n, replace = TRUE)

X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
W = rbeta(n, clust.alpha[cluster], 1 - clust.alpha[cluster])
Y = X %*% beta + rnorm(n, 0, 1) + 2 * W * X[,2]

tau.hat = average_partial_effect(X, Y, W)
print(paste("true tau:", round(mean(2 * cluster.center[,2]), 2)))
print(paste("point estimate:", round(tau.hat[1], 2)))
print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))

```

#### References
David Hirshberg and Stefan Wager.
<b>Augmented Minimax Linear Estimation.</b>
<i>The Annals of Statistics</i>, 49(6), 2021.
[<a href="http://arxiv.org/pdf/1712.00038.pdf">arxiv</a>]
