# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# install.packages(c("R.oo","compoisson","R.methodsS3"))

library(readr)
library(ape)
library(geiger)
library(Biostrings)
library(phylolm)
library(dplyr)
library(gdata)

#rank_fh is a function that takes a dataframe and add a rank column according to one column (given the column number)
#such that the numbers above 0 are ranked from high to low as 1,2,3,4 while numbers below 0 are ranked from
#low to high as -1, -2, -3, -4. 0 is ranked as 0.
rank_fh <- function(df, number) {
  df <-
    df_ord <- df[order(df[, number], decreasing = TRUE), ]
  i <- 1
  for (a in df_ord[, number]) {
    if (a > 0) {
      df_ord$rank[i] <- i
    } else if (a == 0) {
      df_ord$rank[i] <- 0
    } else {
      df_ord$rank[i] <- i - nrow(df_ord) - 1
    }
    i <- i + 1
  }
  return(df_ord)
}

k <- 17

# Input tree
phy <- read.tree(file = "Cohen/TBcohen_k17_n1_whole.tre")
p <- Ntip(phy) #number of tips/species 
Vphy <- vcv(phy)[sort(phy$tip.label), sort(phy$tip.label)]
C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier
#Tried to use the distance matrix as the vcv matrix but the results are not good
#dist <- read.table("~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_whole.dist", row.names=1, quote="\"", comment.char="",skip = 1)
#colnames(dist) <- rownames(dist)
#C <- as.matrix(exp(-dist))

iC <- solve(C) #inverse of C
ones <- array(1, c(p, 1)) # array of ones for the intercept

#Introduce true positives: 
# 1. mutation_kmers are kmers surrounding the three drug-resistent mutation sites (V170X, I491X and S493X),
#    including kmers with and without mutations)
# 2. mutation_kmerPatterns are present(1) and absent(0) patterns of the mutation kmers in those 337 TB strains.
mutation_kmers <- read.csv("~/Dropbox/Research/GWAS/Cohen/rpoB_mutations_k17.kmer", header = F, col.names =c('kmer','freq'),sep="", stringsAsFactors = F)
mutation_kmerPatterns <- read.csv("~/Dropbox/Research/GWAS/Cohen/rpoB_mutations_k17.pattern", header = F, col.names =c('kmer','pattern'),sep="", colClasses=c(rep("factor",2)))
# There are fewer kmer in rpoB_mutations_k17.pattern than in rpoB_mutations_k17.kmer because
# kmers that only existed in on strains are not included in the total kmer table, thus not in rpoB_mutations_k17.pattern
nrow(subset(mutation_kmers,!(mutation_kmers$kmer %in% mutation_kmerPatterns$kmer)))
# There are 56 kmers that contain at one of the mutation site but only occurd in one strain.
collapsed_mutationPatterns <- unique(mutation_kmerPatterns$pattern) 
# Only 9 patterns are true positives (excluding those single-strain kmers)

# I will skip the position of mutation in each kmer for now
#kmer_df$position <- rep(c((nchar(kmers) - k + 1):1), p)

# Genetic trait values. 
# These values are imputed from the sequence at those three mutation sites.
# One resistent mutation is sufficient to make a strain resistent.
# Currently there are only 13 resistent strains computed this way, comparing to 220 in
# the cultured results provided by Cohen 2015, the same paper where we got the genome assembly.
w_gen <- read.table(file = 'Cohen/tblastn_hsp_mutations_sequences_genotype.txt',header=T)
# Culture trait values.
w_cul<- read.table(file = 'Cohen/TBCohen_phenotype.txt',header=T)
w <- merge(w_gen,w_cul,by=c('SpecimenID'))  
##############################################
# scoring (Tony's method, rhs being kmer pattern and lhs being the trait)
##############################################

#y(rhs) should be kmer pattern #How did I get the kmer patterns
#X(lhs) should be the trait
#read in the kmer pattern as Y, and calculate score for each y.
Y <- read.fwf(file = "~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_kmerPattern.stats", widths = array(1, c(1, p)), header = F)
colnames(Y) <- sort(phy$tip.label)

#sort the dataframe by tip name, and turn R/S into 1/0
trait_gen <- w$Consensus
trait_gen <- ifelse(trait_gen == 'R',1,0)
names(trait_gen) <- sort(phy$tip.label)
trait_cul <-w$R
trait_cul <- ifelse(trait_cul == 'R',1,0)
names(trait_cul) <- sort(phy$tip.label)

#score based on genetically imputed phenotype
X <- t(t(trait_gen * 1))
xx <- cbind(ones, X)
threshold <- 0

pattern <- array(NA, c(nrow(Y), 1))
output_gen <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scoreGLM=0, scoreLog=0)

for (i in 1:nrow(Y)) {
  y <- t(Y[i, ])
  output_gen$pattern[i] <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  output_gen$sumy[i] <- sumy
  if (sumy >= threshold & sumy <= (p - threshold)) {
    
    ##################
    # LS
    
    XX <- t(xx) %*% xx
    XY <- t(xx) %*% y
    b <- solve(XX, XY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% h/(p - 2)
    iXX <- solve(XX)
    bSE <- (MSE * iXX[2, 2])^0.5
    output_gen$scoreLS[i] <- b[2]/bSE
    
    ##################
    # GLS
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% y
    b <- solve(XiCX, XiCY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% iC %*% h/(p - 2)
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2, 2])^0.5
    output_gen$scoreGLS[i] <- b[2]/bSE
    
    ##################
    # GLM
    # mu = mean(y)
    # B.init = matrix(c(log(mu/(1-mu)),0.0001), ncol=1)
    # show(c(i, 1, sumy, system.time(z.GLM <- binaryPGLM(y ~ X, phy = phy, s2 = 1, B.init = B.init))))
    # if(z.GLM$convergeflag == "converged") {
    #   output$scoreGLM[i] <- z.GLM$B[2]/z.GLM$B.se[2]
    # } else {
    #   show('not converged')
    # }
    ##################
    # GLMM
    
    # show(c(i, 1, sumy, system.time(z.GLMM <- binaryPGLMM(y ~ X, phy=phy, s2.init = 0.001))))
    # if(z.GLMM$convergeflag == "converged") {
    # output$scoreGLM[i] <- z.GLMM$B[2]/z.GLMM$B.se[2]
    # } else {
    # show('not converged')
    # }
    ##################
    # Log
    
    #show(c(i, 2, sumy, system.time(z.Log <- phyloglm(y ~ X, phy=phy))))
    #output$scoreLog[i] <- z.Log$coefficients[2]/z.Log$sd[2]
    
  }
}
#	show(plot(output[abs(output[,3]) < 1000,3:6]))
#	show(plot(output[,4:6]))

ranked.LS.gen <- rank_fh(output_gen, 3)[,c(1,3,7)]
#Get rid of the NAs.
complete.output.gen<- output_gen[complete.cases(output_gen),] 
ranked.GLS.gen <- rank_fh(complete.output.gen, 4)[,c(1,2,7)]
#ranked.GLM <- rank_fh(output, 5)[,c(1,25,7)]
#ranked.Log <- rank_fh(output, 6)[,c(1,6,7)]
colnames(ranked.LS.gen)[3] <- 'rankLS'
colnames(ranked.GLS.gen)[3] <- 'rankGLS'
#colnames(ranked.GLM)[3] <- 'rankGLM'
#colnames(ranked.Log)[3] <- 'rankLog'

###score based on culture determined phenotype###
#################################################
X <- t(t(trait_cul * 1))
xx <- cbind(ones, X)
threshold <- 0

pattern <- array(NA, c(nrow(Y), 1))
output_cul <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scoreGLM=0, scoreLog=0)

for (i in 1:nrow(Y)) {
  y <- t(Y[i, ])
  output_cul$pattern[i] <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  output_cul$sumy[i] <- sumy
  if (sumy >= threshold & sumy <= (p - threshold)) {
    
    ##################
    # LS
    
    XX <- t(xx) %*% xx
    XY <- t(xx) %*% y
    b <- solve(XX, XY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% h/(p - 2)
    iXX <- solve(XX)
    bSE <- (MSE * iXX[2, 2])^0.5
    output_cul$scoreLS[i] <- b[2]/bSE
    
    ##################
    # GLS
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% y
    b <- solve(XiCX, XiCY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% iC %*% h/(p - 2)
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2, 2])^0.5
    output_cul$scoreGLS[i] <- b[2]/bSE
  
  }
}

ranked.LS.cul <- rank_fh(output_cul, 3)[,c(1,3,7)]
complete.output.cul<- output_gen[complete.cases(output_cul),]
ranked.GLS.cul <- rank_fh(complete.output.cul, 4)[,c(1,2,7)]
colnames(ranked.LS.cul)[3] <- 'rankLS'
colnames(ranked.GLS.cul)[3] <- 'rankGLS'

#Save the scores
kmerPattern.stats <- read.table(file = "~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_kmerPattern.stats", colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
total_score.gen <- merge(ranked.LS.gen, kmerPattern.stats, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.GLS.gen, total_score.gen, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.GLM.gen, total_score.gen, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.Log, total_score.gen, by = "pattern", all.x = T)
write.csv(total_score, paste("~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_scores_gen.csv", sep = ""), row.names = FALSE)

total_score.cul <- merge(ranked.LS.cul, kmerPattern.stats, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.GLS.cul, total_score.cul, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.GLM.cul, total_score.cul, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.Log, total_score.cul, by = "pattern", all.x = T)
write.csv(total_score, paste("~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_scores_cul.csv", sep = ""), row.names = FALSE)

#Get scores for true positive patterns
mutation_pattern_score.gen <- merge(mutation_kmerPatterns, total_score.gen, by = "pattern", all.x = T)
mutation_pattern_score.cul <- merge(mutation_kmerPatterns, total_score.cul, by = "pattern", all.x = T)

