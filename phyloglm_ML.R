phyloglm_ML <- function (formula, data = list(), phy, method = c("logistic_MPLE", 
    "logistic_IG10", "poisson_GEE"), btol = 50, log.alpha.bound = 4, 
    start.beta = NULL, start.alpha = NULL, boot = 0, full.matrix = TRUE) 
{
    if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) 
        stop("the tree has no branch lengths.")
    if (is.null(phy$tip.label)) 
        stop("the tree has no tip labels.")
    method = match.arg(method)
    phy = reorder(phy, "pruningwise")
    original.edge.length = phy$edge.length
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    mf = model.frame(formula = formula, data = data)
    if (nrow(mf) != length(phy$tip.label)) 
        stop("the number of rows in the data does not match the number of tips in the tree.")
    if (is.null(rownames(mf))) 
        warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
    else {
        ordr = match(phy$tip.label, rownames(mf))
        if (sum(is.na(ordr)) > 0) 
            stop("the row names in the data do not match the tip labels in the tree.\n")
        mf = mf[ordr, , drop = F]
    }
    X = model.matrix(attr(mf, "terms"), data = mf)
    y = model.response(mf)
    dk = ncol(X)
    dis = pruningwise.distFromRoot(phy)
# browser()
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (sum(!(y %in% c(0, 1)))) 
            stop("The model by Ives and Garland requires a binary response (dependent variable).")
        if (var(y) == 0) 
            stop("the response (dependent variable) is always 0 or always 1.")
        btouch = 0
        proposedBetaSD = 0.05
        D = max(dis[1:n]) - dis[1:n]
        D = D - mean(D)
        externalEdge = (des <= n)
        
        ##########################
        # I deleted this
        phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
        ##########################
        
        times <- pruningwise.branching.times(phy)
        names(times) <- (n + 1):(n + phy$Nnode)
        Tmax <- max(times)
        intern = which(phy$edge[, 2] > n)
        lok = rep(-1, N)
        lok[intern] = des[intern] - n
    }
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
    plogregfunct <- function(startB, startlL) {
        convergeflag = 0
        clL = startlL
        cB = startB
        diflL = 100
        difB = 100
        counter = 0
        ttozero = 10^6
        optss <- list(reltol = .Machine$double.eps^0.5, maxit = 1e+05, 
            parscale = 1)
        while (((diflL > 10^-6) | (difB > 10^-6) | (ttozero > 
            10^-1)) & (counter < 20)) {
            counter = counter + 1
            oldlL = clL
            oldB = cB
            olddiflL = diflL
            olddifB = difB
            opt <- optim(par = clL, fn = function(par) {
                plogreglLfunct(cB, par)
            }, method = "L-BFGS-B")
            clL = as.numeric(opt$par)
            diflL = (clL - oldlL)^2
            if (counter >= 10) 
                clL = (clL + oldlL)/2
            opt <- optim(par = cB, fn = function(par) {
                plogregBfunct(par, clL)
            }, method = "L-BFGS-B", control = list(factr = 1e+12))
            cB = as.vector(opt$par)
            ttozero = as.numeric(opt$value)
            if (ttozero > 10^-2) {
                Btemp = rnorm(dk, startB, proposedBetaSD * pmax(abs(startB), 
                  rep(0.1, dk)))
                opt <- optim(par = Btemp, fn = function(par) {
                  plogregBfunct(par, clL)
                }, method = "L-BFGS-B", control = list(factr = 1e+12))
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
            if ((max(abs(c(oldlL - clL, oldB - cB))) > 0.1) | 
                (ttozero > 0.5)) 
                convergeflag = 1
        return(list(B = cB, lL = clL, convergeflag = convergeflag))
    }
    plogregBfunct <- function(B, lL) {
        if (dk > 1) 
            g = X %*% B
        else g = rep(1, n) * B
        if (any(abs(g) >= btol)) {
            btouch <<- 1
            return(1e+06)
        }
        mu = as.vector(1/(1 + exp(-g)))
        temp = transf.branch.lengths(B, lL)
        dia = temp[[3]]
        comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
            (1 - mu) * X/dia)
        logdetC = comp$logd + 2 * sum(log(dia)) - sum(log(mu * 
            (1 - mu)))
        if (logdetC < -100 * log(10)) 
            return(1e+06)
        Z = comp$Xy
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
                ttemp = transf.branch.lengths(dB, lL)
                tdiag = ttemp[[3]]
                tcomp = three.point.compute(ttemp[1:2], (y - 
                  mu)/tdiag, mu * (1 - mu) * X/tdiag)
                dinfoMp = tcomp$XX
                dB = B
                dB[i] = dB[i] - Dx
                g = X %*% dB
                if (any(abs(g) >= btol)) 
                  return(1e+06)
                mu = as.vector(1/(1 + exp(-g)))
                ttemp = transf.branch.lengths(dB, lL)
                tdiag = ttemp[[3]]
                tcomp = three.point.compute(ttemp[1:2], (y - 
                  mu)/tdiag, mu * (1 - mu) * X/tdiag)
                dinfoMm = tcomp$XX
                DinfoM = (dinfoMp - dinfoMm)/Dx/2
                FirthC[i] = sum(diag(invInfoM %*% DinfoM))/2
            }
        }
        tozero = Z + FirthC
        return(sum(tozero^2))
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
    iterate_beta <- function(beta) {
        difbeta = 1
        maxint = 10000
        count = 0
        curbeta = beta
        while ((difbeta > 1e-10) && (count < maxint)) {
            mu = as.vector(exp(X %*% curbeta))
            temp = transf.branch.lengths_poisson_GEE(curbeta)
            dia = temp[[3]]
            if (sum(which(mu == 0)) > 0) 
                break
            comp = three.point.compute(temp[1:2], (y - mu)/dia, 
                mu * X/dia)
            invI = solve(comp$XX)
            newbeta = curbeta + invI %*% comp$Xy
            count = count + 1
            difbeta = sum(abs(newbeta - curbeta))
            curbeta = newbeta
        }
        mu = as.vector(exp(X %*% curbeta))
        r = (y - mu)/sqrt(mu)
        phi = sum(r^2)/(n - dk)
        covBSE = phi * invI
        BSE = sqrt(diag(covBSE))
        if (difbeta > 1e-10) 
            convergeflag = 1
        else convergeflag = 0
        return(list(beta = as.vector(curbeta), BSE = BSE, covBSE = covBSE, 
            phi = phi, convergeflag = convergeflag))
    }
    if (is.null(start.beta)) {
        if (method %in% c("logistic_MPLE", "logistic_IG10")) {
# browser()
			if (dk > 1) {
				fit = logistf(y ~ X - 1)
				startB = fit$coefficients
			}else{
				fit = mean(y)
				startB = fit
			}
			
			
            # if (any(abs(X %*% startB) >= btol)) {
                # warning("The estimated coefficients in the absence of phylogenetic signal lead\n  to some linear predictors beyond 'btol'. Increase btol?\n  Starting from beta=0 other than intercept.")
                # startB = numeric(dk)
                # iint = match("(Intercept)", colnames(X))
                # if (!is.na(iint)) 
                  # startB[iint] = log(sum(y == 1)/sum(y == 0))
                # if (any(abs(X %*% startB) >= btol)) 
                  # startB[iint] = 0
            # }
        }
    }
    else {
        if (length(start.beta) != dk) 
            stop(paste("start.beta shoudl be of length", dk))
        if (method %in% c("logistic_MPLE", "logistic_IG10")) {
            startB = as.vector(start.beta)
            if (any(abs(X %*% startB) >= btol)) 
                stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
        }
    }
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (method == "logistic_MPLE") {
# browser()
            opt <- optim(par = c(startB), fn = npllh, 
                method = "L-BFGS-B", control = list(factr = 1e+12))
            B = opt$par[1:dk]
            lL = 1
            convergeflag = opt$convergence
        }
        if (btouch == 1) 
            warning("the boundary of the linear predictor has been reached during the optimization procedure.\nYou can increase this bound by increasing 'btol'.")
        plogregBSE = plogregBSEfunct(B, lL)
        results <- list(coefficients = B, 
            sd = plogregBSE$BSE, vcov = plogregBSE$covBSE, convergence = convergeflag)
    }
    if (method == "poisson_GEE") {
        res = iterate_beta(as.vector(start.beta))
        results <- list(coefficients = res$beta, scale = res$phi, 
            sd = res$BSE, vcov = res$covBSE, convergence = res$convergeflag)
    }
    if (results$converge) 
        warning("phyloglm failed to converge.\n")
    names(results$coefficients) = colnames(X)
    colnames(results$vcov) = colnames(X)
    rownames(results$vcov) = colnames(X)
    results$linear.predictors = as.vector(X %*% results$coefficients)
    names(results$linear.predictors) = names(y)
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (max(abs(results$linear.predictors)) + 0.01 > btol) 
            warning("the linear predictor reaches its bound for one (or more) tip.")
        results$fitted.values = as.vector(1/(1 + exp(-results$linear.predictors)))
        results$mean.tip.height = Tmax
        results$logLik = llh(results$fitted.values)
        results$penlogLik = results$logLik + log(det(as.matrix(plogregBSE$info)))/2
        results$aic = -2 * results$logLik + 2 * (dk + 1)
    }
    names(results$fitted.values) = names(y)
    results$residuals = y - results$fitted.values
    results$y = y
    results$n = n
    results$d = dk
    results$formula = formula
    results$call = match.call()
    results$method = method
    results$X = X
    results$boot = boot
    if ((boot > 0) && (method %in% c("logistic_MPLE", "logistic_IG10"))) {
        options(warn = -1)
        bootobject <- rbinTrait(n = boot, phy = phy, beta = results$coefficients, 
            alpha = results$alpha, X = X, model = "LogReg")
        ncoeff = length(results$coefficients)
        bootmatrix <- matrix(NA, boot, ncoeff + 1)
        colnames(bootmatrix) <- c(names(results$coefficients), 
            "alpha")
        for (i in 1:boot) {
            y = bootobject[, i]
            if (method == "logistic_IG10") {
                bootfit <- try(plogregfunct(startB, startlL), 
                  silent = TRUE)
                if (!inherits(bootfit, "try-error")) {
                  bootmatrix[i, 1:ncoeff] <- bootfit$B
                  bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$lL)
                }
            }
            if (method == "logistic_MPLE") {
                bootfit <- try(optim(par = c(startB, startlL), 
                  fn = npllh, method = "L-BFGS-B", control = list(factr = 1e+12)), 
                  silent = TRUE)
                if (!inherits(bootfit, "try-error")) {
                  bootmatrix[i, 1:ncoeff] <- bootfit$par[1:dk]
                  bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$par[dk + 
                    1])
                }
            }
        }
        ind.na <- which(is.na(bootmatrix[, 1]))
        if (length(ind.na) > 0) {
            bootmatrix <- bootmatrix[-ind.na, ]
            numOnes <- range(apply(bootobject[, ind.na], 2, sum))
        }
        bootmean <- apply(bootmatrix, 2, mean)
        bootsd <- apply(bootmatrix, 2, sd)
        bootconfint95 <- apply(bootmatrix, 2, quantile, probs = c(0.025, 
            0.975))
        bootmeanAlog <- mean(log(bootmatrix[, ncoeff + 1]))
        bootsdAlog <- sd(log(bootmatrix[, ncoeff + 1]))
        results$bootmean = bootmean
        results$bootsd = bootsd
        results$bootconfint95 = bootconfint95
        results$bootmeanAlog = bootmeanAlog
        results$bootsdAlog = bootsdAlog
        results$bootnumFailed = length(ind.na)
        if (full.matrix) 
            results$bootstrap = bootmatrix
        options(warn = 0)
    }
    class(results) = "phyloglm"
    results
}