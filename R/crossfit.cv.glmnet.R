crossfit.cv.glmnet <-
function (x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse", 
                                                           "deviance", "class", "auc", "mae"), nfolds = 10, foldid,
                                                           lambda.choice = c("lambda.1se", "lambda.min"),
            grouped = TRUE, keep = FALSE, parallel = FALSE, ...) 
{
  lambda.choice = match.arg(lambda.choice)

  if (missing(type.measure)) 
    type.measure = "default"
  else type.measure = match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.glmnet")
  N = nrow(x)
  if (missing(weights)) 
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped", 
    "keep"), names(glmnet.call), F)
  if (any(which)) 
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  glmnet.object = glmnet(x, y, weights = weights, offset = offset, 
    lambda = lambda, ...)
  glmnet.object$call = glmnet.call
  is.offset = glmnet.object$offset
###Next line is commented out so each call generates its own lambda sequence
# lambda=glmnet.object$lambda
 if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  }
  else nz = sapply(predict(glmnet.object, type = "nonzero"), 
         length)
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
#  if (parallel && require(foreach)) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
    {
      which = foldid == i
      # if (is.matrix(y))
      if (length(dim(y))>1) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
             offset = offset_sub, weights = weights[!which], 
             ...)
    }
  }
  else {
    for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y)) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE], 
               y_sub, lambda = lambda, offset = offset_sub, 
               weights = weights[!which], ...)
    }
  }
  # Choose lambda with the exported cv.glmnet, on the same folds. Earlier
  # versions called the glmnet internals cvtype(), cv.elnet() and getmin(),
  # which no longer exist in current releases of glmnet.
  cv.args = list(x = x, y = y, weights = weights, offset = offset,
                 lambda = glmnet.object$lambda, foldid = foldid,
                 grouped = grouped, keep = keep, ...)
  if (type.measure != "default")
    cv.args$type.measure = type.measure
  obj = do.call(glmnet::cv.glmnet, cv.args)
  
  lambda.opt = as.numeric(obj[lambda.choice])
  cv.betas.list = lapply(outlist, function(elem) coef(elem, s = lambda.opt))
  cv.betas = Reduce(rbind, lapply(cv.betas.list, FUN=as.numeric))
  
  obj[[length(obj) + 1]] = cv.betas
  names(obj)[length(obj)] = "cv.betas"
  
  obj[[length(obj) + 1]] = foldid
  names(obj)[length(obj)] = "foldid"
  
  obj[[length(obj) + 1]] = x
  names(obj)[length(obj)] = "x"
  
  obj
}

