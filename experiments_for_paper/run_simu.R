rm(list = ls())

# Trick to check if we're on sherlock
sherlock = nchar(Sys.getenv("GROUP_SCRATCH")) > 0

if (!sherlock) {
  setup = 1
  n = 200
  p = 3
  sigma = 1
  k = 2
  NREP = 2
} else {
  setwd("/home/users/vitorh/amlinear/experiments_for_paper/")
  setup = sample(c(1, 2, 3, 4), 1)
  n = sample(c(200, 400, 800, 1600), 1)
  p = sample(c(6, 12), 1)
  sigma = sample(c(1), 1)
  k = sample(c(3, 4), 1)
  NREP = 10
}

source("utils2.R")
source("simulators.R")

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

# Saves filename with random suffix bit a the end
uniqstr = paste0(c(Sys.getenv(c('SLURM_JOB_ID', 'SLURM_LOCALID', 'SLURM_JOB_NAME')), sample(1e8, 1)), sep="")
fnm = paste("results/output", setup, n, p, sigma, k, NREP, "full", uniqstr, ".csv", sep="-")
write.csv(results, file=fnm)

