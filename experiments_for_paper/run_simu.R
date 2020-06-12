rm(list = ls())
source("utils2.R")
source("simulators.R")

args=(commandArgs(TRUE))
setup = as.numeric(args[1])
n = as.numeric(args[2])
p = as.numeric(args[3])
sigma = as.numeric(args[4])
k = as.numeric(args[5])
NREP = as.numeric(args[6])

results.list = lapply(1:NREP, function(iter) {
    X = matrix(rnorm(n*p), n, p)
    params = simulators[[setup]](X, k)
    W = params$w.fun()
    Y = params$mu + W * params$tau + sigma * params$sigma.mult * rnorm(n)
    btau = mean(params$tau)
    gamma.oracle = (W - params$w.mean) / params$w.var
    cmp = ape.comparison(X, Y, W, gamma.oracle = gamma.oracle)
    cmp$ape = params$ape
    cmp$cape = btau
    print(cmp)
    cmp
})

results = Reduce(rbind, results.list)

fnm = paste("results/output", setup, n, p, sigma, k, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)

print(colMeans(results^2))
print(colMeans((results - results[1,1])^2))
