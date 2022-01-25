phyloglm_ML <- function (formula, data = list(), phy, btol = 100, start.beta = NULL, ultrametric = T){
  # this is used in PLogML, PLogML_LR and PLogML_LR_p
  # ultrametric=T is by default
  require("logistf") 
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")
  phy = reorder(phy, "pruningwise")
  original.edge.length = phy$edge.length
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  
  mf = model.frame(formula = formula, data = data)
  if (nrow(mf) != length(phy$tip.label)) stop("the number of rows in the data does not match the number of tips in the tree.")
  if (is.null(rownames(mf))) {
    warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
  } else {
    ordr = match(phy$tip.label, rownames(mf))
    if (sum(is.na(ordr)) > 0) 
      stop("the row names in the data do not match the tip labels in the tree.\n")
    mf = mf[ordr, , drop = F]
  }
  
  X = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  dk = ncol(X)
  dis = pruningwise.distFromRoot(phy)
  
  if (sum(!(y %in% c(0, 1)))) stop("The model by Ives and Garland requires a binary response (dependent variable).")
  if (var(y) == 0) stop("the response (dependent variable) is always 0 or always 1.")
  
  btouch = 0
  proposedBetaSD = 0.05
  
  externalEdge = (des <= n)
  D = max(dis[1:n]) - dis[1:n]
  # D = D - mean(D)
  if(ultrametric) {
    phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
  }
  
  times <- pruningwise.branching.times(phy)
  names(times) <- (n + 1):(n + phy$Nnode)
  Tmax <- max(times)
  intern = which(phy$edge[, 2] > n)
  lok = rep(-1, N)
  lok[intern] = des[intern] - n
  
  transf.branch.lengths <- function(B, lL) {
    if (dk > 1) 
      g = X %*% B
    else g = rep(1, n) * B
    mu = as.vector(1/(1 + exp(-g)))
    p = mean(mu)
    alpha = 1/exp(lL)
    edge.length = numeric(N)
    distFromRoot <- exp(-2 * alpha * times)
    tmp = .C("transbranchlengths_IvesGarland2010", as.integer(N), 
             as.integer(des), as.integer(anc - n), as.integer(lok), 
             as.double(distFromRoot), as.integer(externalEdge), 
             as.double(mu), as.double(p), as.double(alpha), as.double(D), 
             el = as.double(1:N), di = as.double(1:n))
    edge.length = tmp$el
    diag = tmp$di
    root.edge = min(distFromRoot)
    if (any(is.nan(edge.length))) 
      stop("edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.")
    return(list(edge.length, root.edge, diag))
  }
  three.point.compute <- function(trans, y, X) {
    ole = 4 + 2 * dk + dk * dk
    tmp = .C("threepoint", as.integer(N), as.integer(n), 
             as.integer(phy$Nnode), as.integer(1), as.integer(dk), 
             as.integer(ROOT), as.double(trans[[2]]), as.double(trans[[1]]), 
             as.integer(des), as.integer(anc), as.double(as.vector(y)), 
             as.double(as.vector(X)), result = double(ole))$result
    return(list(vec11 = tmp[2], y1 = tmp[3], yy = tmp[4], 
                X1 = tmp[5:(4 + dk)], XX = matrix(tmp[(5 + dk):(ole - 
                                                                  dk)], dk, dk), Xy = tmp[(ole - dk + 1):ole], 
                logd = tmp[1]))
  }
  plogregBSEfunct <- function(B, lL) {
    g = X %*% B
    mu = as.vector(1/(1 + exp(-g)))
    temp = transf.branch.lengths(B, lL)
    dia = temp[[3]]
    comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
                                 (1 - mu) * X/dia)
    infoM = comp$XX
    covBSE = solve(infoM)
    BSE = sqrt(diag(covBSE))
    return(list(BSE = BSE, covBSE = covBSE, info = infoM))
  }
  npllh <- function(par) {
    g = X %*% par[1:dk]
    if (any(abs(g) >= btol)) {
      btouch <<- 1
      return(1e+10)
    }
    mu = as.vector(1/(1 + exp(-g)))
    temp = transf.branch.lengths(par[1:dk], 1)
    dia = temp[[3]]
    comp = three.point.compute(temp[1:2], numeric(n), mu * 
                                 (1 - mu) * X/dia)
    infoM = comp$XX
    llk <- .C("logistreglikelihood", as.integer(N), as.integer(n), 
              as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
              as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
              as.double(as.vector(mu)), as.integer(dk), as.double(1), loglik = double(1))$loglik
    if (dk == 1) 
      pllik = llk + log(abs(infoM))/2
    else pllik = llk + log(det(infoM))/2
    -pllik
  }
  llh <- function(mu) {
    .C("logistreglikelihood", as.integer(N), as.integer(n), 
       as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
       as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
       as.double(as.vector(mu)), as.integer(dk), as.double(1), 
       loglik = double(1))$loglik
  }
  
  if (is.null(start.beta)) {
    if(dk>1){
      fit = logistf(y ~ X - 1)
      startB = fit$coefficients
    }else{
      startB = mean(y)
    }
  } else {
    if (length(start.beta) != dk) stop(paste("start.beta should be of length", dk))
    startB = as.vector(start.beta)
    if (any(abs(X %*% startB) >= btol)) stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
  }
  
  opt <- optim(par = c(startB), fn = npllh, method = "L-BFGS-B", control = list(factr = 1e+12))
  B = opt$par[1:dk]
  lL = 1
  convergeflag = opt$convergence
  
  if (btouch == 1) warning("the boundary of the linear predictor has been reached during the optimization procedure.\nYou can increase this bound by increasing 'btol'.")
  plogregBSE = plogregBSEfunct(B, lL)
  results <- list(coefficients = B, sd = plogregBSE$BSE, vcov = plogregBSE$covBSE, convergence = convergeflag)
  
  if (results$converge) warning("phyloglm failed to converge.\n")
  
  names(results$coefficients) = colnames(X)
  colnames(results$vcov) = colnames(X)
  rownames(results$vcov) = colnames(X)
  results$linear.predictors = as.vector(X %*% results$coefficients)
  names(results$linear.predictors) = names(y)
  
  if (max(abs(results$linear.predictors)) + 0.01 > btol) warning("the linear predictor reaches its bound for one (or more) tip.")
  results$fitted.values = as.vector(1/(1 + exp(-results$linear.predictors)))
  results$mean.tip.height = Tmax
  results$logLik = llh(results$fitted.values)
  results$penlogLik = results$logLik + log(det(as.matrix(plogregBSE$info)))/2
  results$aic = -2 * results$logLik + 2 * (dk + 1)
  
  names(results$fitted.values) = names(y)
  results$residuals = y - results$fitted.values
  results$y = y
  results$n = n
  results$d = dk
  results$X = X
  results
}
