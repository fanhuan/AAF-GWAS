#read in the tree and convert it into a vcv
library(ape)

#AAFtree <- read.tree(file = '/Users/iveslab/Documents/projects/GWAS/Yeast/yeast32_asb_k17.tre')

AAFtree <- read.tree(file = '/Users/arives/Dropbox/traveller\ dropbox/Traveller\ AAF/Traveller\ GWAS/yeast32_asb_k17.tre')
AAFtree_VCV <- vcv(AAFtree)

#sort the vcv according to strain names alphabetically
h <- order(colnames(AAFtree_VCV))
VCV <- AAFtree_VCV[h,h]
p <- nrow(VCV) #number of species

#read in the trait value (sporulation rate of yeast)
X <- c(0,0,0,0.613,0.224,0.404,0,0.674,0.06,0,0.199,0.766,0.702,0,0,0.924,0.868,0.885,0.852,0.986,0.899,0.54,0.737,0.221,0,0.482,0.402,0,0.99,0.979,0,0)
X <- as.matrix(X, nrow=p, ncol=1)

#read in the 01 patterns (Y!)
#Y <- read.fwf(file='/Users/iveslab/Documents/projects/GWAS/Yeast/yeast32_asb_k17_patternOnly.txt', widths=array(1,c(1,p)),header=F)

########################################################################
# Huan, in this file, and the species in alphabetical order?
########################################################################
Y <- read.fwf(file='/Users/arives/Dropbox/traveller\ dropbox/Traveller\ AAF/Traveller\ GWAS/yeast32_asb_k17_patternOnly.txt', widths=array(1,c(1,p)),header=F)
dim(Y) # check dimensions of Y

#n, number of kmer patterns.
n <- nrow(Y) #STOPPED HERE
# [1] 113955

# compute scores
C <- as.matrix(VCV)
C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier

iC <- solve(C) #inverse of C
ones <- array(1,c(p,1)) # array of ones for the intercept
output <- array(0,c(n,1)) # collects the results for the n patterns
threshold <- 3 # minimum number of zeros or ones for includion of pattern
for(i in 1:n){
  y <- as.matrix(Y[i,]) 
  sumy <- sum(y) #sum of the pattern
  
  if (sumy < threshold | sumy > (p-threshold)) {
  	output[i] <- 0 #all patterns with only 1 one or zero are assigned score 0.
  	}  else{
    xx <- cbind(ones,X)
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% t(y)
    
    b <- solve(XiCX,XiCY)
    
    h <- t(y) - (xx %*% b)
    
    MSE <- t(h) %*% iC %*% h/(p-2)
    
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2,2])^.5
    output[i] <- b[2]/bSE
  }
}
par(mfrow=c(1,2))
hist(abs(output),breaks=50, main=NULL, xlab='Abs(scores)', ylab='Count')
plot(abs(output), type="l")

# proportion of patterns with scores > 1.96
mean(abs(output) > 1.96)

# patterns with the highest and lowest scores
minY <- t(Y[which.min(output),])
maxY <- t(Y[which.max(output),])
cbind(X, minY, maxY)

# switch direction for minY
cbind(X, 1-minY, maxY)

output.VCV <- output

# save results
Pattern <- data.frame(score=output, pattern=Y)
name.order <- order(colnames(AAFtree_VCV))
names(Pattern)[2:(1+p)] <- colnames(AAFtree_VCV)[name.order]

# This file include the species names in the first row; change col.names=F if you don't want them
write.table(Pattern,file='/Users/arives/Dropbox/traveller\ dropbox/Traveller\ AAF/Traveller\ GWAS/yeastKmer_GWAScore_29Nov14.txt', sep="\t",col.names=T, row.names=F, quote=F)

#############################################################################################################
# I wanted to look at the scores in more detail

score.order <- order(abs(output), decreasing=T)
ordered.scores <- output[score.order]
ordered.patterns <- Y[score.order,]

cbind(rbind(0,X),t(cbind(ordered.scores[1:10], ordered.patterns[1:10,])))
cor(ordered.scores[1:10], ordered.patterns[1:10,])


#############################################################################################################
# I'm going to try this with the empirical covariance matrix; you can ignore this.
#############################################################################################################

# construct empirical covariance matrix
C <- as.matrix(cov(Y))

# compute scores
iC <- solve(C) #inverse of C
ones <- array(1,c(p,1)) # array of ones for the intercept
output <- array(0,c(n,1)) # collects the results for the n patterns
threshold <- 3 # minimum number of zeros or ones for includion of pattern
for(i in 1:n){
  y <- as.matrix(Y[i,]) 
  sumy <- sum(y) #sum of the pattern
  
  if (sumy < threshold | sumy > (p-threshold)) {
  	output[i] <- 0 #all patterns with only 1 one or zero are assigned score 0.
  	}  else{
    xx <- cbind(ones,X)
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% t(y)
    
    b <- solve(XiCX,XiCY)
    
    h <- t(y) - (xx %*% b)
    
    MSE <- t(h) %*% iC %*% h/(p-2)
    
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2,2])^.5
    output[i] <- b[2]/bSE
  }
}
par(mfrow=c(1,2))
hist(abs(output),breaks=50, main=NULL, xlab='Abs(scores)', ylab='Count')
plot(abs(output), type="l")

# proportion of patterns with scores > 1.96
mean(abs(output) > 1.96)

# patterns with the highest and lowest scores
minY <- t(Y[which.min(output),])
maxY <- t(Y[which.max(output),])
cbind(X, minY, maxY)

# switch direction for minY
cbind(X, 1-minY, maxY)

output.empirical <- output

#############################################################################################################
# I'm going to try this without phylogenetic correlations among species
#############################################################################################################

iC <- diag(p) #inverse of C
ones <- array(1,c(p,1)) # array of ones for the intercept
output <- array(0,c(n,1)) # collects the results for the n patterns
threshold <- 3 # minimum number of zeros or ones for includion of pattern
for(i in 1:n){
  y <- as.matrix(Y[i,]) 
  sumy <- sum(y) #sum of the pattern
  
  if (sumy < threshold | sumy > (p-threshold)) {
  	output[i] <- 0 #all patterns with only 1 one or zero are assigned score 0.
  	}  else{
    xx <- cbind(ones,X)
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% t(y)
    
    b <- solve(XiCX,XiCY)
    
    h <- t(y) - (xx %*% b)
    
    MSE <- t(h) %*% iC %*% h/(p-2)
    
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2,2])^.5
    output[i] <- b[2]/bSE
  }
}
par(mfrow=c(1,2))
hist(abs(output),breaks=50, main=NULL, xlab='Abs(scores)', ylab='Count')
plot(abs(output), type="l")

# proportion of patterns with scores > 1.96
mean(abs(output) > 1.96)

# patterns with the highest and lowest scores
minY <- t(Y[which.min(output),])
maxY <- t(Y[which.max(output),])
cbind(X, minY, maxY)

# switch direction for minY
cbind(X, 1-minY, maxY)

output.star <- output

par(mfrow=c(1,3))
plot(output.empirical, output.VCV, pch=20)
plot(output.empirical, output.star, pch=20)
plot(output.star, output.VCV, pch=20)

