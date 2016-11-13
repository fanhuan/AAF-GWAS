# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# install.packages(c("R.oo","compoisson","R.methodsS3"))
# install.packages("/Users/arives/Downloads/phylosim_2.1.1.tar.gz")

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
#read.tree only reads 5 digits after the decimal point, therefore we are going to scale our tree up so there won't be so many edges with a length of 0.
#But we cannot use tree_scaler. I need to write a python script for it.

p <- Ntip(phy) #number of tips/species 

#Use the distance matrix as the vcv matrix.
#Vphy <- vcv(phy)[sort(phy$tip.label), sort(phy$tip.label)]
dist <- read.table("~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_whole.dist", row.names=1, quote="\"", comment.char="",skip = 1)
colnames(dist) <- rownames(dist)
C <- as.matrix(exp(-dist))
#C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier
iC <- solve(C) #inverse of C
ones <- array(1, c(p, 1)) # array of ones for the intercept

#Load in the shared kmer table
#If no header:
#phylokmer <- read.delim(paste("rpoB_GY84_", repn, "/phylokmer.dat", sep = ""), header = FALSE)
#With header:
phylokmer <- read.delim("Cohen/TBcohen_k17_n1_test.dat",header = FALSE,skip=2+p)
colnames(phylokmer) <- c("kmer", sort(phy$tip.label))

#Grep the patterns for kmers invloving S450 (kmer_list)
mutation_kmers <- read.csv("~/Dropbox/Research/GWAS/Cohen/rpoB_mutations_k17.kmer", col.names =c('kmer','freq'),sep="", stringsAsFactors = F)
mutation_kmerPatterns <- read.csv("~/Dropbox/Research/GWAS/Cohen/rpoB_mutations_k17.pattern", col.names =c('kmer','pattern'),sep="", colClasses=c(rep("factor",2)))
collapsed_mutationPatterns <- unique(mutation_kmerPatterns$pattern) 
# I will skip position for now
#kmer_df$position <- rep(c((nchar(kmers) - k + 1):1), p)

# Genetic trait values.
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
output <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scoreGLM=0, scoreLog=0)

for (i in 1:nrow(Y)) {
  y <- t(Y[i, ])
  output$pattern[i] <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  output$sumy[i] <- sumy
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
    output$scoreLS[i] <- b[2]/bSE
    
    ##################
    # GLS
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% y
    b <- solve(XiCX, XiCY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% iC %*% h/(p - 2)
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2, 2])^0.5
    output$scoreGLS[i] <- b[2]/bSE
    
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

ranked.LS <- rank_fh(output, 3)[,c(1,3,7)]
write.table(ranked.LS,'~/Dropbox/Research/GWAS/Cohen/rankedLS.txt',quote = F,row.names = F)
complete_output<- output[complete.cases(output),]
ranked.GLS <- rank_fh(complete_output, 4)[,c(1,2,7)]
#ranked.GLM <- rank_fh(output, 5)[,c(1,25,7)]
#ranked.Log <- rank_fh(output, 6)[,c(1,6,7)]
colnames(ranked.LS)[3] <- 'rankLS'
colnames(ranked.GLS)[3] <- 'rankGLS'
#colnames(ranked.GLM)[3] <- 'rankGLM'
#colnames(ranked.Log)[3] <- 'rankLog'

#read in the kmer pattern again, focus on the frequency this time.
kmerPattern.stats <- read.table(file = "~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_kmerPattern.stats", colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
total_score <- merge(ranked.LS, kmerPattern.stats, by = "pattern", all.x = T)
total_score <- merge(ranked.GLS, total_score, by = "pattern", all.x = T)
total_score <- merge(ranked.GLM, total_score, by = "pattern", all.x = T)
total_score <- merge(ranked.Log, total_score, by = "pattern", all.x = T)
write.csv(total_score, paste("~/Dropbox/Research/GWAS/Cohen/TBcohen_k17_n1_p1_scores.csv", sep = ""), row.names = FALSE)

#Get scores for patterns of kmers in drug resistent gene .
mutation_pattern_score <- merge(mutation_kmerPatterns, ranked.GLS, by = "pattern", all.x = T)
#Note that some of the kmers does not have scores because only patterns with more than two 1 or 0 are scored.
#Merge S450_kmer_score with kmer_df

output_450 <- merge(S450_kmers_score, kmer_df, by = "kmer", all.x = T)
output_450_light <- subset(output_450, select = c("kmer", "position", "pattern", "freq", "scoreLS", "rankLS", "scoreGLS", "rankGLS", "scoreGLM", "rankGLM", "scoreLog", "rankLog"))
write.csv(output_450_light, paste("rpoB_GY84_", as.character(repn), "_450summary.csv", sep = ""), row.names = FALSE)

outputplot <- output[abs(output$scoreLS) < 1000,]
outputplot <- outputplot[order(abs(outputplot$scoreGLS)),]
col450 <- is.element(outputplot$pattern, output_450_light$pattern)

show(plot(outputplot[,3:6], col=(1+col450), pch=(20 - col450), cex=(.5 + 1*col450)))

pdf(paste("rpoB_GY84_", as.character(repn), "_450summary.pdf", sep = ""), width=6, height=6)
par(mfrow=c(2,2))
for(ii in c(6,8,10,12)) hist(output_450_light[,ii], main=colnames(output_450_light)[ii], xlab="rank")
dev.off()
par(mfrow=c(1,1))
}