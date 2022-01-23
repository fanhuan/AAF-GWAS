# This script calculates association scores for the Cohen dataset
##############################################


##############################################
# Load packages
##############################################
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
library(plyr)
library(dplyr)
library(gdata)
library(ggplot2)
library(parallel)
source("~/Dropbox/Research/GWAS/phyloglm_aaf.R")
source("~/Dropbox/Research/GWAS/phyloglm_ML.R")


##############################################
# Self-defined packages
##############################################
#rank_fh is a function that takes a dataframe and add a rank column according to one column (given the column number)
#such that the numbers above 0 are ranked from high to low as 1,2,3,4 while numbers below 0 are ranked from
#low to high as -1, -2, -3, -4. 0 is ranked as 0.
rank_fh <- function(serial) {
  position <- seq(length(serial))
  df <- data.frame(serial, position)
  df_NA <- df %>% filter(is.na(serial))
  if (nrow(df_NA) > 0) {
    df_NA$ranking = NA
  }
  df_Inf <- df %>% filter(is.infinite(serial))
  if (nrow(df_Inf) > 0) {
    df_Inf$ranking = NA
  }
  df_pos <- df %>% filter(serial > 0)
  if (nrow(df_pos) > 0) {
    df_pos <- df_pos %>% arrange(serial)
    df_pos$ranking <- rank(-df_pos$serial, ties.method = 'min')
  }
  df_zero <- df %>% filter(serial==0)
  if (nrow(df_zero) > 0) {
    df_zero$ranking = 0
  }
  df_neg <- df %>% filter(serial < 0)
  if (nrow(df_neg) > 0) {
    df_neg <- df_neg %>% arrange(serial)
    df_neg$ranking <- rank(df_neg$serial, ties.method = 'min')
  }
  df <- rbind(df_NA,df_pos,df_zero,df_neg)
  df <- df[order(df$position),]
  return(df$ranking)
}

LS <- function(xx,y){
  XX <- t(xx) %*% xx
  XY <- t(xx) %*% y
  b <- solve(XX, XY)
  h <- y - (xx %*% b)
  MSE <- t(h) %*% h/(p - 2)
  iXX <- solve(XX)
  bSE <- (MSE * iXX[2, 2])^0.5
  scoreLS <- b[2]/bSE
  return(scoreLS)
}

Logistf <- function(X,y){
  fit_Log <- logistf(y ~ X)
  scoreLog <- sign(fit_Log$coef[2]) * qchisq(fit_Log$prob[2], df=1, lower.tail=F)
  return(scoreLog)
}

GLS <- function(xx,y,iC){
  XiCX <- t(xx) %*% iC %*% xx
  XiCY <- t(xx) %*% iC %*% y
  b <- solve(XiCX, XiCY)
  h <- y - (xx %*% b)
  MSE <- t(h) %*% iC %*% h/(p - 2)
  iXiCX <- solve(XiCX)
  bSE <- (MSE * iXiCX[2, 2])^0.5
  scoreGLS <- b[2]/bSE
  return(scoreGLS)
}

PLog <- function(X,y,phy.extend){
  fit <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=F)
  scorePLog <- fit$zB[2]
  return(scorePLog)
}

PLogF <- function(X,y,phy.extend){
  fit_F <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=T)
  scorePLogF <- fit_F$zB[2]
  return(scorePLogF)
}

PLogML <- function(X,y,phy.extend){
  fit_ML <- try(phyloglm_ML(y ~ X, phy=phy.extend), silent=T)
  fit_ML0 <- try(phyloglm_ML(y ~ 1, phy=phy.extend), silent=T)
  if(is.null(attr(fit_ML,'condition')) & is.null(attr(fit_ML0,'condition'))){
    scorePLogML <- sign(fit_ML$coef[2]) * 2*(fit_ML$logLik - fit_ML0$logLik)
  }else{
    scorePLogML <- NA
  }
  return(scorePLogML)
}

addscore <- function(x) {
  return(paste0('score',x))
}

score_fun <- function(i,Y,X,threshold,phy.extend){
  y <- t(Y[i, ])
  xx <- cbind(ones, X)
  pattern <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  if (sumy >= threshold & sumy <= (p - threshold)) {
    ##################
    # LS
    scoreLS = LS(xx,y)
    ###################
    # Logistf
    scoreLogistf = Logistf(X,y)
    #################
    # GLS
    scoreGLS = GLS(xx,y,iC)
    ##################
    # PLog
    scorePLog = PLog(X,y,phy.extend)
    ##################
    # PLogF
    scorePLogF = PLogF(X,y,phy.extend)
    ##################
    # PLogML
    scorePLogML = PLogML(X,y,phy.extend)
    return(c(pattern,sumy,scoreLS,scoreLogistf,scoreGLS,scorePLog,scorePLogF,scorePLogML))
  }
}

#Parameters
k <- 19
threshold <- 3
#Models available
model_list <- c('LS','Logistf','GLS','PLog','PLogF','PLogML')


##############################################
# Input tree
##############################################

phy <- read.tree(file = paste0("Cohen/TBcohen_k",k,"_n1_whole_sorted.tre"))
#Test whether there is zero length branches. If yes, collapse them.
if (min(phy$edge.length) == 0) {
  threshold <- min(phy$edge.length[which(phy$edge.length>0)])
  phy <- di2multi(phy,threshold)
}
#Extend short tips
phy.extend <- phy
if(any(phy$edge.length[phy$edge[,2] <= Ntip(phy)] < 10^-6)) {
  tip.length <- phy$edge.length[phy$edge[,2] <= Ntip(phy)]
  tip.length <- tip.length[tip.length > 10^-6]
  min_tip <- min(tip.length)
  phy.extend$edge.length[phy.extend$edge[,2] <= Ntip(phy)]<- min_tip + phy.extend$edge.length[phy.extend$edge[,2] <= Ntip(phy)]
}
p <- Ntip(phy.extend) #number of tips/species
Vphy <- vcv(phy.extend)[sort(phy.extend$tip.label), sort(phy.extend$tip.label)]
# this stabilizes the values so that the determinant isn't zero
Vphy <- Vphy/mean(Vphy)

C <- Vphy
C <- C/det(C)^(1/p) 
#some transformation of C to make the determinant = 1, which should make calculations easier
#Tried to use the distance matrix as the vcv matrix but the results are not good
#dist <- read.table("TBcohen_k17_n1_whole.dist", row.names=1, quote="\"", comment.char="",skip = 1)
#colnames(dist) <- rownames(dist)
#C <- as.matrix(exp(-dist))

iC <- solve(C) #inverse of C
ones <- array(1, c(p, 1)) # array of ones for the intercept
##############################################
# scoring (Tony's method, rhs being kmer pattern and lhs being the trait)
##############################################
#X(lhs) should be the trait
# Genetic trait values.
# These values are imputed from the sequence at those three mutation sites.
# One resistent mutation is sufficient to make a strain resistent.
# Currently there are only 13 resistent strains computed this way, comparing to 220 in
# the cultured results provided by Cohen 2015, the same paper where we got the genome assembly.
w_gen <- read.table(file = 'Cohen/tblastn_hsp_425_451_genotype.txt',header=F,col.names = c('SpecimenID','Genotype'))

# Culture trait values.
w_cul <- read.table(file = 'Cohen/TBCohen_phenotype.txt',header=T)
w <- merge(w_gen,w_cul,by=c('SpecimenID'),sort = F)
#sort the dataframe by tip name, and turn R/S into 1/0
trait_gen <- w$Genotype
trait_gen <- ifelse(trait_gen == 'R',1,0)
names(trait_gen) <- sort(phy.extend$tip.label)
trait_cul <-w$R
trait_cul <- ifelse(trait_cul == 'R',1,0)
names(trait_cul) <- sort(phy$tip.label)

##############################################
# input patterns
##############################################
#y(rhs) should be kmer pattern #How did I get the kmer patterns
#read in the kmer pattern as Y, and calculate score for each y.
Y <- read.fwf(file = paste0("Cohen/TBcohen_k",k,"_n1_p1_kmerPattern.stats"), widths = array(1, c(1, p)), header = F)
colnames(Y) <- sort(phy.extend$tip.label)

##############################################
# score based on genetically imputed phenotype
X <- t(t(trait_gen * 1))

reps <- seq(nrow(Y))
r<-sapply(reps, function(x){score_fun(x, Y=Y,X=X,threshold=threshold,
                                      phy.extend=phy.extend)})
output_gen <- data.frame(t(r),stringsAsFactors = F)
asd<-sapply(model_list,addscore,USE.NAMES = F)
output_gen[,2:ncol(output_gen)] <- as.data.frame(sapply(output_gen[,2:ncol(output_gen)], as.numeric))
names(output_gen) = c('pattern','sumy',asd)
output_gen <- output_gen %>% filter(sumy > 1)

reps <- seq(nrow(Y))  
r<-sapply(reps, function(x){score_fun(x, Y=Y,X=X,threshold=threshold,
                                      phy.extend=phy.extend)})
output_cul <- data.frame(t(r),stringsAsFactors=FALSE)

#kmers from RRDR regions
#Introduce true positives:
# 1. mutation_kmers are kmers surrounding the three drug-resistent mutation sites (V170X, I491X and S493X),
#    including kmers with and without mutations)
# 2. mutation_kmerPatterns are present(1) and absent(0) patterns of the mutation kmers in those 337 TB strains.
#mutation_kmers <- read.csv(paste0("Cohen/RRDR_k",k,".kmer"), header = F, col.names =c('kmer','freq'),sep="", stringsAsFactors = F)
mutation_kmerPatterns <- read.csv(paste0("Cohen/RRDR_k",k,".pattern"), header = F, col.names =c('kmer','pattern'),sep="", colClasses=c(rep("factor",2)))
# There might be fewer kmer in pattern than in kmer if p != 1
# because kmers that only existed in on strains are not included in the total kmer table (now fixed)
# nrow(subset(mutation_kmers,!(mutation_kmers$kmer %in% mutation_kmerPatterns$kmer)))
# There are 56 kmers that contain at least one of the mutation site but only occurred in one strain.
RRDR_patterns <- unique(mutation_kmerPatterns$pattern)

#Pair plot
output_gen$correct <- is.element(output_gen$pattern,RRDR_patterns)
pdf(paste0('Cohen_k',k,'_scores.pdf'),width = 10, height = 10)

plot(output_gen[,3:(ncol(output_gen)-1)], col=(1+output_gen$correct), pch=(20 - output_gen$correct), cex=(.5 + 1*output_gen$correct))
dev.off()
# 45(k17) or 50(k19) patterns are true positives (excluding those single-strain kmers)

#Put ranks for each model

output_gen$rankLS <- rank_fh(output_gen$scoreLS)
output_gen$rankLog <- rank_fh(output_gen$scoreLog)
output_gen$rankGLS <- rank_fh(output_gen$scoreGLS)
output_gen$rankPLog <- rank_fh(output_gen$scorePLog)
output_gen$rankPLogF <- rank_fh(output_gen$scorePLogF)
output_gen$rankPLogML <- rank_fh(output_gen$scorePLogML)

# I will skip the position of mutation in each kmer for now
#kmer_df$position <- rep(c((nchar(kmers) - k + 1):1), p)
#Save the scores
output_gen_RRDR <- subset(output_gen,correct==T,-correct)
write.csv(output_gen, paste0("TBcohen_k",k,"_n1_p2_scores_gen.csv", sep = ""), row.names = FALSE)
write.csv(output_gen_RRDR, paste("TBcohen_k",k,"_n1_p1_scores_RRDR_gen.csv", sep = ""), row.names = FALSE)

##############################################
# score based on culture determined phenotype###
X <- t(t(trait_cul * 1))

#output_cul = data.frame(matrix(vector(),nrow(Y),2+length(model_list)))
#names(output_cul) = c('pattern','sumy',asd)
reps <- seq(nrow(Y))  
r<-sapply(reps, function(x){score_fun(x, Y=Y,X=X,threshold=threshold,
                                   phy.extend=phy.extend)})
output_cul <- data.frame(t(r),stringsAsFactors=FALSE)
output_cul[,2:ncol(output_cul)] <- as.data.frame(sapply(output_cul[,2:ncol(output_cul)], as.numeric))
names(output_cul) = c('pattern','sumy',asd)
output_cul <- output_cul %>% filter(sumy > 1)

#Add rank    
output_cul$rankLS <- rank_fh(output_cul$scoreLS)
output_cul$rankLog <- rank_fh(output_cul$scoreLog)
output_cul$rankGLS <- rank_fh(output_cul$scoreGLS)
output_cul$rankPLog <- rank_fh(output_cul$scorePLog)
output_cul$rankPLogF <- rank_fh(output_cul$scorePLogF)
output_cul$rankPLogML <- rank_fh(output_cul$scorePLogML)
#output_cul <- output_cul %>% dplyr::select(-id)
#Save the scores
write.csv(output_cul, paste0("TBcohen_k",k,"_n1_p2_scores_cul.csv"), row.names = FALSE)

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
