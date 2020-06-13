library(EQL)
library(CVXR)
library(grf)

# generate the ball L^2_s(gamma_d), where gamma_d = N(0,1)
# up do given order of hermite basis
# X is n by d
# if derivative[j], provide basis for derivative of basis element
get.basis = function(X, order, s, derivative=c()) {
  H = lapply(1:ncol(X), function(j) {
    Hj = sapply(1:order, function(k) hermite(X[,j], k, prob = TRUE)) / sqrt(factorial(j))
    if(j %in% derivative) { cbind(1, Hj[,1:(order-1)] %*% diag(2:order)) } 
    else { Hj }
  })
  polys = lapply(0:order, function(r) {
    if (r == 0){
      val = if(length(derivative)==0) { 1 } else { 0 }
      return(rep(val, nrow(X)))
    }
    partitions = combn(r + ncol(X) -1, ncol(X) - 1,
                       function(vec) c(vec, r + ncol(X)) - c(0, vec) - 1)
    elems = sapply(1:ncol(partitions), function(iter) {
      part = partitions[,iter]
      if(any(part[derivative]==0)) { return(rep(0,nrow(X))) }
      idx = which(part > 0)
      elem = H[[idx[1]]][,part[idx[1]]]
      if (length(idx) > 1) {
        for (id in idx[-1]) {
          elem = elem * H[[id]][,part[id]]
        }
      }
      elem * (1 + r)^(-s/2)
    })
    elems
  })
  Reduce(cbind, polys)
}

# check derivative by finite differencing
check.derivative.fd = function() {
    x=rnorm(6)
    dim(x)=c(3,2)
    dx = x*0; dx[,1]=1;
    eps=1e-6

    db=get.basis(x,order=5,s=2, derivative=1)
    b=get.basis(x,order=5,s=2)
    bp=get.basis(x+eps*dx, order=5,s=2)
    rel.error = abs((bp-b)/eps - db) / (abs(db)+eps^2)
    max(rel.error)
}

# minimax weights for a function class F := { sum_j phi_j(z) beta_j : || beta[!exact] ||_p <= 1 }.
# fBasis (n x p) with fBasis[i,j] = phi_j(Z_i)
# hBasis of length p with hBasis[j] = sum_i h(Z_i, phi_j) 
# p, a scalar
# sigma, a scalar
# exact, a list of the 'unbounded directions' in F for which we should have exact balance
# returns the minimax balancing weights
minimax.weights = function(fBasis, hBasis, p, sigma=1, exact=c()) {
  is.exact = 1:ncol(fBasis) %in% exact
  n = nrow(fBasis)
  gamma = CVXR::Variable(n)
  objective = CVXR::p_norm( hBasis[!is.exact] - t(fBasis[,!is.exact]) %*% gamma,p)^2 + sigma^2 + CVXR::norm2(t(gamma))^2
  constraints = if(length(exact) == 0) { list() } else { list(hBasis[is.exact] == t(fBasis[,is.exact]) %*% gamma) }
  cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
  cvx.output = solve(cvx.problem, solver = "ECOS", verbose = FALSE)
  result = cvx.output$getValue(gamma)
  return(result)
}


#crossfit.rlasso = function(Y, X, W, fold, ff) {
#    tau.fit = rlasso(X[fold != ff,],  Y[fold != ff], W[fold != ff])
#    Y.fit = cv.glmnet(X[fold != ff,], Y[fold != ff])
#    W.fit = cv.glmnet(X[fold != ff,], W[fold != ff])
#    list(tau = X[fold == ff,] %*% tau.fit$tau.beta,
#	 m = predict(Y.fit, X[fold==ff,]) + (W[fold == ff]-predict(Y.fit, X[fold == ff, ]))*tau.hat)
#}

crossfit.zero = function(Y,X,W,fold) {
    list(tau = rep(0,length(Y)), m = rep(0, length(Y)))
}
crossfit.causal.forest = function(Y,X,W,fold) {
    stopifnot(all(fold %in% c(1,2)))
    tau.hat = rep(NA, length(Y))
    m.hat = rep(NA, length(Y))
    for(ff in 1:2) {
	forest.Y = regression_forest(X[fold != ff,],Y[fold != ff])
	forest.W = regression_forest(X[fold != ff,],W[fold != ff])
	forest.tau = causal_forest(X[fold != ff,], Y[fold != ff], W[fold != ff], 
				   Y.hat = predict(forest.Y)$predictions, W.hat = predict(forest.W)$predictions)
    
	Y.hat = predict(forest.Y, X[fold == ff,])[,1]
	W.hat = predict(forest.W, X[fold == ff,])[,1]
	tau.hat[fold == ff] = predict(forest.tau, X[fold == ff,])[,1]
	m.hat[fold == ff]   = Y.hat + (W[fold == ff] - W.hat)*tau.hat[fold == ff]
    }
    list(tau=tau.hat, m=m.hat)
}


amle.average.derivative = function(Y, X, W, sigma, 
				   p = 1, s = (ncol(X)+1)/2 + 1, order=4, 
				   crossfit.regression=crossfit.causal.forest) {
    fBasis = get.basis(cbind(W,X), order, s)
    hBasis = get.basis(cbind(W,X), order, s, derivative=1)

    fold = 1 + (sample(1:length(W)) > length(W)/2)
    outcome.model = crossfit.regression(Y, X, W, fold)
    gamma = rep(NA, length(Y))
    for(ff in 1:2) {
	fBasis.aug = cbind(outcome.model$m[fold == ff], fBasis[fold == ff, ])
	hBasis.aug = c(sum(outcome.model$tau[fold == ff]), colSums(hBasis[fold == ff,]))
	gamma[fold == ff] = minimax.weights(fBasis.aug, hBasis.aug, p, sigma, exact=c(1,2))
    }
    psi.hat = mean(outcome.model$tau + gamma * (Y-outcome.model$m))
    V.hat.sample =  mean(gamma^2 * (Y - outcome.model$m)^2)
    V.hat.pop = V.hat.sample + mean(psi.hat - outcome.model$tau)^2
    attr(psi.hat, 'se.sample') = sqrt(V.hat.sample/length(Y))
    attr(psi.hat, 'se.pop') = sqrt(V.hat.pop/length(Y))
    psi.hat
}

amle.cond.linear = function(Y, X, W, sigma, 
				   p = 1, s = (ncol(X)+1)/2 + 1, order=4, 
				   crossfit.regression=crossfit.causal.forest) {
    fBasisX = get.basis(X, order, s)
    fBasis = cbind(fBasisX,    W*fBasisX)
    hBasis = cbind(0*fBasisX, fBasisX)

    fold = 1 + (sample(1:length(W)) > length(W)/2)
    outcome.model = crossfit.regression(Y, X, W, fold)
    gamma = rep(NA, length(Y))
    for(ff in 1:2) {
	fBasis.aug = cbind(outcome.model$m[fold == ff], fBasis[fold == ff, ])
	hBasis.aug = c(sum(outcome.model$tau[fold == ff]), colSums(hBasis[fold == ff,]))
	gamma[fold == ff] = minimax.weights(fBasis.aug, hBasis.aug, p, sigma, exact=c(1,2))
    }
    psi.hat = mean(outcome.model$tau + gamma * (Y-outcome.model$m))
    V.hat.sample =  mean(gamma^2 * (Y - outcome.model$m)^2)
    V.hat.pop = V.hat.sample + mean(psi.hat - outcome.model$tau)^2
    attr(psi.hat, 'se.sample') = sqrt(V.hat.sample/length(Y))
    attr(psi.hat, 'se.pop') = sqrt(V.hat.pop/length(Y))
    psi.hat
}


comparison = function(X, Y, W, order=NULL, gamma.oracle = rep(1, length(W))) {
    forest.tau = causal_forest(X,Y,W)
    forest.estimate = average_partial_effect(forest.tau)
    psi.forest = forest.estimate[[1]]
    attr(psi.forest, 'se.sample') = forest.estimate[[2]]
    
    tau.hat = predict(forest.tau, X)[,1]
    m.hat = forest.tau$Y.hat + (W - forest.tau$W.hat)*tau.hat 
    psi.oracle = mean(tau.hat + gamma.oracle * (Y - m.hat))
    attr(psi.oracle, 'se.sample') = sqrt(mean(gamma.oracle^2 * (Y - m.hat)^2)/length(Y))

    estimates = list(psi.oracle, psi.forest,
		     amle.average.derivative(Y,X,W,1,1, crossfit.regression=crossfit.zero),
		     amle.average.derivative(Y,X,W,1,1, crossfit.regression=crossfit.causal.forest),
		     amle.cond.linear(Y,X,W,1,1, crossfit.regression=crossfit.zero),
		     amle.cond.linear(Y,X,W,1,1, crossfit.regression=crossfit.causal.forest))
    weights = c('oracle', 'grf', 'deriv', 'deriv', 'clin', 'clin')
    outcome = c('grf',    'grf', 'none',   'grf',  'none', 'grf')
    data.frame(weights = weights, outcome = outcome, estimate=sapply(estimates, c), se.sample = sapply(estimates, function(e)attr(e,'se.sample')))
}
