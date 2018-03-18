rm(list = ls())

library(amlinear)
source("utils.R")

args=(commandArgs(TRUE))
setup = as.numeric(args[1])
n = as.numeric(args[2])
p = as.numeric(args[3])
sigma = as.numeric(args[4])
k = as.numeric(args[5])
NREP = as.numeric(args[6])

if (setup == 1) {
    
    get.params = function(X, k) {
        a = rowMeans(X[,1:k]) * sqrt(k)
        x = sign(a) * a^2
        tau = -0.2
        w.mean = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
        mu = x - tau * (w.mean - 0.5)
        alpha = w.mean
        beta = (1 - w.mean)
        w.var = alpha * beta / (alpha + beta)^2 / (1 + alpha * beta)
        w.fun = function() rbeta(nrow(X), alpha, beta)
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = -0.2)
    }
    
} else if (setup == 2) {
    
    k = k - 2
    get.params = function(X, k) {
        pr = apply(2 * X[,1:k,drop=FALSE], 1, prod)
        mu = sign(pr) * sqrt(abs(pr)) / 2
        w.mean = 0.1 * sign(mu) + mu
        w.var = w.mean^2
        w.fun = function() w.mean * (1 + rnorm(nrow(X)))
        tau = pmax(X[,1] + X[,2], 0)/ 2
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = 0.5 / sqrt(pi))
    }
    
} else if (setup == 3) {
    
    get.params = function(X, k) {
        tau = rowMeans(cos(pi * X[,1:k] / 3))
        prob = 0.2 + tau^2
        mu = 4 * rowMeans(X) + 2 * prob
        w.fun = function() rpois(nrow(X), lambda = prob)
        list(mu=mu, tau=tau, w.mean=prob, w.var=prob, w.fun=w.fun, sigma.mult = 1, ape = 0.58)
    }
    
} else if (setup == 4) {
    
    k = k + 1
    get.params = function(X, k) {
        tau = sin(2 * pi * X[,1])
        rowm = rowMeans(X[,1:k]* sqrt(k))
        ln.mu = 1 / (1 + exp(-sign(rowm) * rowm^2))
        ln.sigma = 1/3
        w.mean = exp(ln.mu + ln.sigma^2/2)
        w.var = (exp(ln.sigma^2) - 1) * w.mean^2
        w.fun = function() exp(ln.mu + ln.sigma * rnorm(nrow(X)))
        mu = pmax(0, 2 * rowm)
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = 0)
    }
    
} else {
    
    stop("bad setup")
    
}

results.list = lapply(1:NREP, function(iter) {
    X = matrix(rnorm(n*p), n, p)
    params = get.params(X, k)
    W = params$w.fun()
    Y = params$mu + W * params$tau + sigma * params$sigma.mult * rnorm(n)
    btau = mean(params$tau)
    gamma.oracle = (W - params$w.mean) / params$w.var
    cmp = ape.comparison(X, Y, W, gamma.oracle = gamma.oracle)
    all.cmp = c(ape=params$ape, cape=btau, cmp)
    print(all.cmp)
    all.cmp
})

results = Reduce(rbind, results.list)

fnm = paste("results/output", setup, n, p, sigma, k, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)

print(colMeans(results^2))
print(colMeans((results - results[1,1])^2))
