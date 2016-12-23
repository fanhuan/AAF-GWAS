# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# install.packages(c("R.oo","compoisson","R.methodsS3"))
# install.packages(c("readr","geiger","dplyr","gdata"))

library(readr)
library(ape)
library(geiger)
library(Biostrings)
library(phylolm) #Ho's package
library(logistf) #When it does not converge
library(dplyr)
library(gdata)
library(ggplot2)

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

k <- 19


##############################################
# Input tree
##############################################
phy <- read.tree(file = paste0("Cohen/TBcohen_k",k,"_n1_whole_sorted.tre"))
#Test whether there is zero length branches. If yes, collapse them.
if (min(phy$edge.length) == 0) {
  threshold <- min(phy$edge.length[which(phy$edge.length>0)])
  phy <- di2multi(phy,threshold)
}



p <- Ntip(phy) #number of tips/species 
Vphy <- vcv(phy)[sort(phy$tip.label), sort(phy$tip.label)]
det(Vphy)

# this stabilizes the values so that the determinant isn't zero
Vphy <- Vphy/mean(Vphy) 
det(Vphy)

C <- Vphy
C <- C/det(C)^(1/p) 
#some transformation of C to make the determinant = 1, which should make calculations easier
#Tried to use the distance matrix as the vcv matrix but the results are not good
#dist <- read.table("TBcohen_k17_n1_whole.dist", row.names=1, quote="\"", comment.char="",skip = 1)
#colnames(dist) <- rownames(dist)
#C <- as.matrix(exp(-dist))

iC <- solve(C) #inverse of C
det(C)
ones <- array(1, c(p, 1)) # array of ones for the intercept

##############################################
# input patterns
##############################################

#Introduce true positives: 
# 1. mutation_kmers are kmers surrounding the three drug-resistent mutation sites (V170X, I491X and S493X),
#    including kmers with and without mutations)
# 2. mutation_kmerPatterns are present(1) and absent(0) patterns of the mutation kmers in those 337 TB strains.
mutation_kmers <- read.csv(paste0("Cohen/RRDR_k",k,".kmer"), header = F, col.names =c('kmer','freq'),sep="", stringsAsFactors = F)
mutation_kmerPatterns <- read.csv(paste0("Cohen/RRDR_k",k,".pattern"), header = F, col.names =c('kmer','pattern'),sep="", colClasses=c(rep("factor",2)))

# There might be fewer kmer in pattern than in kmer if p != 1
# because kmers that only existed in on strains are not included in the total kmer table (now fixed)
nrow(subset(mutation_kmers,!(mutation_kmers$kmer %in% mutation_kmerPatterns$kmer)))

# There are 56 kmers that contain at one of the mutation site but only occurred in one strain.
collapsed_mutationPatterns <- unique(mutation_kmerPatterns$pattern) 

# 45(k17) or 50(k19) patterns are true positives (excluding those single-strain kmers)

# I will skip the position of mutation in each kmer for now
#kmer_df$position <- rep(c((nchar(kmers) - k + 1):1), p)

# Genetic trait values. 
# These values are imputed from the sequence at those three mutation sites.
# One resistent mutation is sufficient to make a strain resistent.
# Currently there are only 13 resistent strains computed this way, comparing to 220 in
# the cultured results provided by Cohen 2015, the same paper where we got the genome assembly.
w_gen <- read.table(file = 'Cohen/tblastn_hsp_425_451_genotype.txt',header=F,col.names = c('SpecimenID','Genotype'))

# Culture trait values.
w_cul <- read.table(file = 'Cohen/TBCohen_phenotype.txt',header=T)
w <- merge(w_gen,w_cul,by=c('SpecimenID'),sort = F)  

##############################################
# scoring (Tony's method, rhs being kmer pattern and lhs being the trait)
##############################################

#y(rhs) should be kmer pattern #How did I get the kmer patterns
#X(lhs) should be the trait
#read in the kmer pattern as Y, and calculate score for each y.
Y <- read.fwf(file = paste0("Cohen/TBcohen_k",k,"_n1_p1_kmerPattern.stats"), widths = array(1, c(1, p)), header = F)
colnames(Y) <- sort(phy$tip.label)

#sort the dataframe by tip name, and turn R/S into 0/1
trait_gen <- w$Genotype
trait_gen <- ifelse(trait_gen == 'R',0,1)
names(trait_gen) <- sort(phy$tip.label)
trait_cul <-w$R
trait_cul <- ifelse(trait_cul == 'R',0,1)
names(trait_cul) <- sort(phy$tip.label)

##############################################
# score based on genetically imputed phenotype
X <- t(t(trait_gen * 1))
xx <- cbind(ones, X)
threshold <- 3

pattern <- array(NA, c(nrow(Y), 1))
output_gen <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scoreGLM=0,scorephyloglm=0)


for (i in 1:nrow(Y)) {
  y <- t(Y[i, ])
  output_gen$pattern[i] <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  output_gen$sumy[i] <- sumy
  if (sumy >= threshold & sumy <= (p - threshold)) {
    
    ##################
    # LS

    # XX <- t(xx) %*% xx
    # XY <- t(xx) %*% y
    # b <- solve(XX, XY)
    # h <- y - (xx %*% b)
    # MSE <- t(h) %*% h/(p - 2)
    # iXX <- solve(XX)
    # bSE <- (MSE * iXX[2, 2])^0.5
    # output_gen$scoreLS[i] <- b[2]/bSE
    # 
    # ##################
    # # GLS:with phylogenetic signal
    # 
    # XiCX <- t(xx) %*% iC %*% xx
    # XiCY <- t(xx) %*% iC %*% y
    # b <- solve(XiCX, XiCY)
    # h <- y - (xx %*% b)
    # MSE <- t(h) %*% iC %*% h/(p - 2)
    # iXiCX <- solve(XiCX)
    # bSE <- (MSE * iXiCX[2, 2])^0.5
    # output_gen$scoreGLS[i] <- b[2]/bSE
    # 
    ##################
    # GLM
    fit <- phyloglm_aaf(y ~ X, phy=phy)
    output_gen$scoreGLM[i] <-fit$zB
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
    # z.Log <- phyloglm(y ~ X, phy=phy, method = "logistic_IG10")
    # if (z.Log$convergence == 0) {
    #   count_phyloglm = count_phyloglm + 1
    #   output_gen$scoreLog[i] <- output_gen$scorephyloglm[i] <- z.Log$coefficients[2]/z.Log$sd[2]
    #   }
    # else { 
    #   l.Log <- logistf(y ~ X)
    #   output_gen$scoreLog[i] <- output_gen$scorelogistf[i] <- l.Log$coefficients[2]/sqrt(l.Log$var[2,2])
    # }
    # show(c(count, z.Log$convergence))
    # count = count + 1
  }
}
#Check for NaNs
c(min(output_gen$scoreLS),min(output_gen$scoreGLS),min(output_gen$scoreGLM))
c(max(output_gen$scoreLS),max(output_gen$scoreGLS),max(output_gen$scoreGLM))
#Count the number of NaNs.
sum(is.na(output_gen$scoreGLM))
#If no NA
ranked.LS.gen <- rank_fh(output_gen, 3)[,c(1,3,7)]
ranked.GLS.gen <- rank_fh(output_gen, 4)[,c(1,4,7)]

#If there are NAs.
complete.output.gen <- output_gen[complete.cases(output_gen),] 
ranked.LS.gen <- rank_fh(complete.output.gen, 3)[,c(1,3,7)]
ranked.GLS.gen <- rank_fh(complete.output.gen, 4)[,c(1,4,7)]
ranked.GLM.gen <- rank_fh(complete.output.gen, 5)[,c(1,5,7)]
#ranked.GLM <- rank_fh(output, 5)[,c(1,5,7)]
#ranked.Log <- rank_fh(output, 6)[,c(1,6,7)]
colnames(ranked.LS.gen)[3] <- 'rankLS'
colnames(ranked.GLS.gen)[3] <- 'rankGLS'
colnames(ranked.GLM.gen)[3] <- 'rankGLM'
#colnames(ranked.Log)[3] <- 'rankLog'

##############################################
# score based on culture determined phenotype###
X <- t(t(trait_cul * 1))
xx <- cbind(ones, X)
threshold <- 1

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
#Check for NaNs
c(min(output_cul$scoreLS),min(output_cul$scoreGLS))
c(max(output_cul$scoreLS),max(output_cul$scoreGLS))
#If no NA
ranked.LS.cul <- rank_fh(output_cul, 3)[,c(1,3,7)]
ranked.GLS.cul <- rank_fh(output_cul, 4)[,c(1,4,7)]
#If NA
complete.output.cul <- output_cul[complete.cases(output_cul),]
ranked.LS.cul <- rank_fh(complete.output.cul, 3)[,c(1,3,7)]
ranked.GLS.cul <- rank_fh(complete.output.cul, 4)[,c(1,4,7)]
#add column names
colnames(ranked.LS.cul)[3] <- 'rankLS'
colnames(ranked.GLS.cul)[3] <- 'rankGLS'

#Save the scores
kmerPattern.stats <- read.table(file = paste0("Cohen/TBcohen_k",k,"_n1_p1_kmerPattern.stats"), colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
total_score.gen <- merge(ranked.LS.gen, kmerPattern.stats, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.GLS.gen, total_score.gen, by = "pattern", all.x = T)
#total_score.gen <- merge(ranked.GLM.gen, total_score.gen, by = "pattern", all.x = T)
#total_score.gen <- merge(ranked.Log, total_score.gen, by = "pattern", all.x = T)
write.csv(total_score.gen, paste0("TBcohen_k",k,"_n1_p1_scores_gen.csv"), row.names = FALSE)

total_score.cul <- merge(ranked.LS.cul, kmerPattern.stats, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.GLS.cul, total_score.cul, by = "pattern", all.x = T)
#total_score.cul <- merge(ranked.GLM.cul, total_score.cul, by = "pattern", all.x = T)
#total_score.cul <- merge(ranked.Log, total_score.cul, by = "pattern", all.x = T)
write.csv(total_score.cul, paste0("TBcohen_k",k,"_n1_p1_scores_cul.csv"), row.names = FALSE)

#Get scores for true positive patterns
mutation_pattern_score.gen <- merge(mutation_kmerPatterns, total_score.gen, by = "pattern", all.x = T)
mutation_pattern_score.cul <- merge(mutation_kmerPatterns, total_score.cul, by = "pattern", all.x = T)
write.csv(mutation_pattern_score.gen, paste0("TBcohen_k",k,"_n1_p1_scores_RRDR_gen.csv", sep = ""), row.names = FALSE)
write.csv(mutation_pattern_score.cul, paste("TBcohen_k",k,"_n1_p1_scores_RRDR_cul.csv", sep = ""), row.names = FALSE)
mutation_pattern_score.gen_k17<-mutation_pattern_score.gen
mutation_pattern_score.cul_k17<-mutation_pattern_score.cul
#Drop the first collumn: patterns
#mutation_pattern_score.gen<-mutation_pattern_score.gen[,2:dim(mutation_pattern_score.gen)[2]]
#mutation_pattern_score.cul<-mutation_pattern_score.cul[,2:dim(mutation_pattern_score.cul)[2]]

#Save patterns with top 10 rank (both positive and negative)
TBcohen_n1_p1_scores_gen <- read.csv(paste0("~/Dropbox/Research/GWAS/TBcohen_k",k,"_n1_p1_scores_gen.csv"),colClasses=c("character",rep("numeric",5)))
TBcohen_n1_p1_scores_cul <- read.csv(paste0("~/Dropbox/Research/GWAS/TBcohen_k",k,"_n1_p1_scores_cul.csv"),colClasses=c("character",rep("numeric",5)))
rankGLS_pattern<-TBcohen_n1_p1_scores_gen %>% filter(abs(rankGLS) < 11) %>% filter(abs(rankGLS) > 0) %>% select(pattern)
rankGLS_pattern_cul<-TBcohen_n1_p1_scores_cul %>% filter(abs(rankGLS) < 11) %>% filter(abs(rankGLS) > 0) %>% select(pattern)
rankLS_pattern_cul<-TBcohen_n1_p1_scores_cul %>% filter(abs(rankLS) < 11) %>% filter(abs(rankLS) > 0) %>% select(pattern)

write.table(rankGLS_pattern, paste("rankGLS_k",k,".pattern", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rankGLS_pattern_cul, paste("rankGLS_k",k,"_cul.pattern", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rankLS_pattern_cul, paste("rankLS_k",k,"_cul.pattern", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Get kmers from those top 10 rank patterns
