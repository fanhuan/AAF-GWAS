binaryPGLM <- function (formula, data = list(), phy, s2 = 0.1, B.init = NULL, 
    tol.pql = 10^-6, maxit.pql = 200) 
{
    if (!inherits(phy, "phylo")) 
        stop("Object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) 
        stop("The tree has no branch lengths.")
    if (is.null(phy$tip.label)) 
        stop("The tree has no tip labels.")
    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    mf <- model.frame(formula = formula, data = data)
    if (nrow(mf) != length(phy$tip.label)) 
        stop("Number of rows of the design matrix does not match with length of the tree.")
    if (is.null(rownames(mf))) {
        warning("No tip labels, order assumed to be the same as in the tree.\n")
        data.names = phy$tip.label
    }
    else data.names = rownames(mf)
    order <- match(data.names, phy$tip.label)
    if (sum(is.na(order)) > 0) {
        warning("Data names do not match with the tip labels.\n")
        rownames(mf) <- data.names
    }
    else {
        tmp <- mf
        rownames(mf) <- phy$tip.label
        mf[order, ] <- tmp[1:nrow(tmp), ]
    }
    X <- model.matrix(attr(mf, "terms"), data = mf)
    y <- model.response(mf)
    if (sum(!(y %in% c(0, 1)))) {
        stop("binaryPGLM requires a binary response (dependent variable).")
    }
    if (var(y) == 0) {
        stop("The response (dependent variable) is always 0 or always 1.")
    }
    p <- ncol(X)
    Vphy <- vcv(phy)
    Vphy <- Vphy/max(Vphy)
    Vphy/exp(determinant(Vphy)$modulus[1]/n)
    if (!is.null(B.init) & length(B.init) != p) {
        warning("B.init not correct length, so computed B.init using glm()")
    }
    if (is.null(B.init) | (!is.null(B.init) & length(B.init) != 
        p)) {
        B.init <- t(matrix(glm(formula = formula, data = data, 
            family = "binomial")$coefficients, ncol = p))
    }
    B <- B.init
    b <- matrix(0, nrow = n)
    beta <- rbind(B, b)
    mu <- exp(X %*% B)/(1 + exp(X %*% B))
    XX <- cbind(X, diag(1, nrow = n, ncol = n))
    C <- s2 * Vphy
    est.B <- B
    oldest.B <- matrix(10^6, nrow = length(est.B))
    iteration <- 0
    exitflag <- 0
    rcondflag <- 0
    while (((t(est.B - oldest.B) %*% (est.B - oldest.B)/length(B) > 
        tol.pql^2)) & (iteration <= maxit.pql)) {
        iteration <- iteration + 1
        oldest.B <- est.B
        est.B.m <- B
        oldest.B.m <- matrix(10^6, nrow = length(est.B))
        iteration.m <- 0
        while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m)/length(B) > 
            tol.pql^2) & (iteration.m <= maxit.pql)) {
            iteration.m <- iteration.m + 1
            oldest.B.m <- est.B.m
            invW <- diag(as.vector((mu * (1 - mu))^-1))
            V <- invW + C
            if (sum(is.infinite(V)) > 0 | rcond(V) < 10^-10) {
                rcondflag <- rcondflag + 1
                B <- 0 * B.init + 0.001
                b <- matrix(0, nrow = n)
                beta <- rbind(B, b)
                mu <- exp(X %*% B)/(1 + exp(X %*% B))
                oldest.B.m <- matrix(10^6, nrow = length(est.B))
                invW <- diag(as.vector((mu * (1 - mu))^-1))
                V <- invW + C
            }
            invV <- solve(V)
            Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
            denom <- t(X) %*% invV %*% X
            num <- t(X) %*% invV %*% Z
            B <- as.matrix(solve(denom, num))
            b <- C %*% invV %*% (Z - X %*% B)
            beta <- rbind(B, b)
            mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
            est.B.m <- B
        }
    }
    convergeflag <- "converged"
    if (iteration >= maxit.pql | rcondflag >= 3) {
        convergeflag <- "Did not converge; try increasing maxit.pql or starting with B.init values of .001"
    }
    converge.test.B <- (t(est.B - oldest.B) %*% (est.B - oldest.B))^0.5/length(est.B)
    invW <- diag(as.vector((mu * (1 - mu))^-1))
    V <- invW + C
    invV <- solve(V)
    Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
    denom <- t(X) %*% invV %*% X
    num <- t(X) %*% invV %*% Z
    B <- solve(denom, num)
    b <- C %*% invV %*% (Z - X %*% B)
    beta <- rbind(B, b)
    mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
    H <- Z - X %*% B
    B.cov <- solve(t(X) %*% invV %*% X)
    B.se <- as.matrix(diag(B.cov))^0.5
    B.zscore <- B/B.se
    B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
    results <- list(formula = formula, B = B, B.se = B.se, B.cov = B.cov, 
        B.zscore = B.zscore, B.pvalue = B.pvalue, s2 = s2, 
        mu = mu, b = b, B.init = B.init, X = X, H = H, VCV = Vphy, 
        V = V, convergeflag = convergeflag, iteration = iteration, converge.test.B = converge.test.B, 
        rcondflag = rcondflag)
    class(results) <- "binaryPGLMM"
    results
}
