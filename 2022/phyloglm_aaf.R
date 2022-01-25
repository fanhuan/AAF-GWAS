phyloglm_aaf <- function (formula, data = list(), phy, btol = 100, start.beta = NULL, full.matrix = TRUE, Firth = T, ultrametric = T) 
{
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
		return(list(vec11 = tmp[2], y1 = tmp[3], yy = tmp[4], X1 = tmp[5:(4 + dk)], XX = matrix(tmp[(5 + dk):(ole - dk)], dk, dk), Xy = tmp[(ole - dk + 1):ole], logd = tmp[1]))
	}
	
	
	plogregfunct <- function(startB) {
		convergeflag = 0
		cB = startB
		difB = 100
		counter = 0
		ttozero = 10^6
		optss <- list(reltol = .Machine$double.eps^0.5, maxit = 1e+05, parscale = 1)
		while ((difB > 10^-6)  & (counter < 20)) {
			counter = counter + 1
			oldB = cB
			olddifB = difB
			opt <- optim(par = cB, fn = plogregBfunct, method = "L-BFGS-B", control = list(factr = 1e+12))
			cB = as.vector(opt$par)
			ttozero = as.numeric(opt$value)
			if (ttozero > 10^-2) {
				Btemp = rnorm(dk, startB, proposedBetaSD * pmax(abs(startB), rep(0.1, dk)))
				opt <- optim(par = Btemp, fn = plogregBfunct, method = "L-BFGS-B", control = list(factr = 1e+12))
				Btemp = as.vector(opt$par)
				newttozero = as.numeric(opt$value)
				if (newttozero < ttozero) {
					cB = Btemp
					ttozero = newttozero
				}
			}
			difB = sum((cB - oldB) * (cB - oldB))
			if (counter >= 10) 
				cB = (cB + oldB)/2
		}
		if (counter >= 19) 
			if ((max(abs(oldB - cB)) > 0.1) | 
					(ttozero > 0.5)) 
				convergeflag = 1
		return(list(B = cB, convergeflag = convergeflag))
	}
	plogregBfunct <- function(B) {
		if (dk > 1) {
			g = X %*% B
		}else{
			g = rep(1, n) * B
		}
		if (any(abs(g) >= btol)) {
			btouch <<- 1
			return(1e+06)
		}
		mu = as.vector(1/(1 + exp(-g)))
		
		# set alpha = 1
		temp = transf.branch.lengths(B, 0)
		dia = temp[[3]]
		comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * (1 - mu) * X/dia)
		logdetC = comp$logd + 2 * sum(log(dia)) - sum(log(mu * (1 - mu)))
		if (logdetC < -100 * log(10)) return(1e+06)
		Z = comp$Xy
		if (Firth == F){
			FirthC <- 0
		} else {
			if (dk == 1) 
				FirthC = (1 - 2 * mu)/2
			else {
				Dx = 0.1
				infoM = comp$XX
				invInfoM = solve(infoM)
				FirthC = rep(NA, dk)
				for (i in 1:dk) {
					dB = B
					dB[i] = dB[i] + Dx
					g = X %*% dB
					if (any(abs(g) >= btol)) 
						return(1e+06)
					mu = as.vector(1/(1 + exp(-g)))
					ttemp = transf.branch.lengths(dB, 0)
					tdiag = ttemp[[3]]
					tcomp = three.point.compute(ttemp[1:2], (y - mu)/tdiag, mu * (1 - mu) * X/tdiag)
					dinfoMp = tcomp$XX
					dB = B
					dB[i] = dB[i] - Dx
					g = X %*% dB
					if (any(abs(g) >= btol)) return(1e+06)
					mu = as.vector(1/(1 + exp(-g)))
					ttemp = transf.branch.lengths(dB, 0)
					tdiag = ttemp[[3]]
					tcomp = three.point.compute(ttemp[1:2], (y - mu)/tdiag, mu * (1 - mu) * X/tdiag)
					dinfoMm = tcomp$XX
					DinfoM = (dinfoMp - dinfoMm)/Dx/2
					FirthC[i] = sum(diag(invInfoM %*% DinfoM))/2
				}
			}
		}
		tozero = Z + FirthC
		return(sum(tozero^2))
	}
	plogregBSEfunct <- function(B) {
		g = X %*% B
		mu = as.vector(1/(1 + exp(-g)))
		temp = transf.branch.lengths(B, 0)
		dia = temp[[3]]
		comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * (1 - mu) * X/dia)
		infoM = comp$XX
		covBSE = solve(infoM)
		BSE = sqrt(diag(covBSE))
		return(list(BSE = BSE, covBSE = covBSE, info = infoM))
	}
	
	if (is.null(start.beta)) {
		# fit = glm(y ~ X - 1, family = binomial)
		fit = logistf(y ~ X - 1)
		startB = fit$coefficients
	} else {
		if (length(start.beta) != dk) stop(paste("start.beta should be of length", dk))
		startB = as.vector(start.beta)
		if (any(abs(X %*% startB) >= btol)) stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
	}

	plogreg = plogregfunct(startB)
	B = plogreg$B
	convergeflag = plogreg$convergeflag
	if (btouch == 1) 
		warning("the boundary of the linear predictor has been reached during the optimization procedure.\nYou can increase this bound by increasing 'btol'.")
	plogregBSE = plogregBSEfunct(B)
	results <- list(coefficients = B, BSE = plogregBSE$BSE, vcov = plogregBSE$covBSE, zB = B/plogregBSE$BSE, convergence = convergeflag)
	if (results$converge) warning("phyloglm failed to converge.\n")
	names(results$coefficients) = colnames(X)
	colnames(results$vcov) = colnames(X)
	rownames(results$vcov) = colnames(X)
	results$linear.predictors = as.vector(X %*% results$coefficients)
	names(results$linear.predictors) = names(y)
	if (max(abs(results$linear.predictors)) + 0.01 > btol) 
		warning("the linear predictor reaches its bound for one (or more) tip.")
	results$fitted.values = as.vector(1/(1 + exp(-results$linear.predictors)))
	names(results$fitted.values) = names(y)
	results$residuals = y - results$fitted.values
	results$y = y
	results$n = n
	results$d = dk
	results$formula = formula
	results$X = X
	results
}



###############################################################
# library(phylolm)

# nspp <- 50
# tre = rtree(nspp)
# x = rTrait(n=1,phy=tre)
# X = cbind(rep(1,nspp),x)
# y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
# dat = data.frame(trait01 = y, predictor = x)
# fit <- phyloglm_aaf(trait01~predictor, phy=tre, data=dat)
# fit$zB

# summary(phyloglm(trait01~predictor, phy=tre, data=dat))

# system.time(fit <- phyloglm_aaf(trait01~predictor, phy=tre, data=dat))
# system.time(fit <- phyloglm(trait01~predictor, phy=tre, data=dat))
# system.time(fit <- phyloglm(trait01~predictor, phy=tre, data=dat, method="logistic_IG10"))


