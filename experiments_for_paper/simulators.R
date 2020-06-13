simulators = list(
    function(X, k, sigma) {
        a = rowMeans(X[,1:k]) * sqrt(k)
        x = sign(a) * a^2
        tau = -0.2
        w.mean = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
        mu = x - tau * (w.mean - 0.5)
        alpha = w.mean
        beta = (1 - w.mean)
        w.var = alpha * beta / (alpha + beta)^2 / (1 + alpha + beta)
        w.fun = function() rbeta(nrow(X), alpha, beta)
	y.fun = function(W) mu + W * tau + sigma * rnorm(n)
	oracle.fun = function(W) (W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = 0.58, cpsi = function(W) mean(tau), estimand='clin')
    },
    function(X, k, sigma) {
	k = k - 2
        pr = apply(2 * X[,1:k,drop=FALSE], 1, prod)
        mu = sign(pr) * sqrt(abs(pr)) / 2
        w.mean = 0.1 * sign(mu) + mu
        w.var = w.mean^2
        w.fun = function() w.mean * (1 + rnorm(nrow(X)))
        tau = pmax(X[,1] + X[,2], 0)/ 2
	y.fun = function(W) mu + W * tau + sigma * rnorm(n)
	oracle.fun = function(W) (W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = 0.58, cpsi = function(W) mean(tau), estimand='clin')
    },
    function(X, k, sigma) {
        tau = rowMeans(cos(pi * X[,1:k] / 3))
        w.mean = 0.2 + tau^2
        mu = 4 * rowMeans(X) + 2 * w.mean
        w.fun = function() rpois(nrow(X), lambda = w.mean)
	y.fun = function(W) mu + W * tau + sigma * rnorm(n)
	oracle.fun = function(W) (W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = 0.58, cpsi = function(W) mean(tau), estimand='clin')
    },
    function(X, k, sigma) {
	k = k + 1
        tau = sin(2 * pi * X[,1])
        rowm = rowMeans(X[,1:k]* sqrt(k))
        ln.mu = 1 / (1 + exp(-sign(rowm) * rowm^2))
        ln.sigma = 1/3
        w.mean = exp(ln.mu + ln.sigma^2/2)
        w.var = (exp(ln.sigma^2) - 1) * w.mean^2
        mu = pmax(0, 2 * rowm)
        w.fun = function() exp(ln.mu + ln.sigma * rnorm(nrow(X)))
	y.fun = function(W) mu + W * tau + sigma * rnorm(n)
	oracle.fun = function(W) (W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = 0, cpsi = function(W) mean(tau), estimand='clin')
    },
    function(X, k, sigma) {
        a = rowMeans(X[,1:k]) * sqrt(k)
        x = sign(a) * a^2
        tau = -0.2
        w.mean = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
        mu = x - tau * (w.mean - 0.5)
        alpha = w.mean
        beta = (1 - w.mean)
        w.var = alpha * beta / (alpha + beta)^2 / (1 + alpha + beta)
        w.fun = function() rbeta(nrow(X), alpha, beta)
	y.fun = function(W) mu + W^3 * tau + sigma * rnorm(n)
	oracle.fun = function(W) 0*(W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = NA, cpsi = function(W) mean(3*W^2*tau), estimand='deriv')
    },
    function(X, k, sigma) {
        a = rowMeans(X[,1:k]) * sqrt(k)
        x = sign(a) * a^2
        tau = -0.2
        w.mean = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
        mu = x - tau * (w.mean - 0.5)
        alpha = w.mean
        beta = (1 - w.mean)
        w.var = alpha * beta / (alpha + beta)^2 / (1 + alpha + beta)
        w.fun = function() rbeta(nrow(X), alpha, beta)
	y.fun = function(W) mu + (W-w.mean)^2 * tau + sigma * rnorm(n)
	oracle.fun = function(W) 0*(W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = NA, cpsi = function(W) mean(2*(W-w.mean)*tau), estimand='deriv')
    },
    function(X, k, sigma) {
	k = k - 2
        pr = apply(2 * X[,1:k,drop=FALSE], 1, prod)
        mu = sign(pr) * sqrt(abs(pr)) / 2
        w.mean = 0.1 * sign(mu) + mu
        w.var = w.mean^2
        w.fun = function() w.mean * (1 + rnorm(nrow(X)))
        tau = pmax(X[,1] + X[,2], 0)/ 2
	y.fun = function(W) mu + sin(20*W) * tau + sigma * rnorm(n)
	oracle.fun = function(W) 0*(W - w.mean) / w.var 
        list(w.fun=w.fun, y.fun = y.fun, oracle.fun = oracle.fun, psi = NA, cpsi = function(W) mean(20*cos(20*W)*tau), estimand='deriv')
    }
)

