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

#Parameters
k <- 17
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
xx <- cbind(ones, X)
threshold <- 3

pattern <- array(NA, c(nrow(Y), 1))
output_gen = data.frame(matrix(vector(),nrow(Y),2+length(model_list)))
asd<-sapply(model_list,addscore,USE.NAMES = F)
names(output_gen) = c('pattern','sumy',asd)
output_gen$pattern<-pattern

for (i in 1:nrow(Y)) {
  y <- t(Y[i, ])
  output_gen$pattern[i] <- paste(y, collapse = "")
  sumy <- sum(y) #sum of the pattern
  output_gen$sumy[i] <- sumy
  if (sumy >= threshold & sumy <= (p - threshold)) {
    
    ##################
    # LS
	
    output_gen$scoreLS[i] <- LS(xx,y)
	
    ###################
    # Logistf

    output_gen$scoreLogistf[i] <- Logistf(X,y)
	
    #################
    # GLS
	
    output_gen$scoreGLS[i] <- GLS(xx,y,iC)
	
	##################
    # PLog
	
	output_gen$scorePLog[i] <- PLog(X,y,phy.extend)
	##################
    # PLogF
	
	output_gen$scorePLogF[i] <- PLogF(X,y,phy.extend)
	##################
	# PLogML
	
	output_gen$scorePLogML[i] <- PLogML(X,y,phy.extend)
	
	}
}

complete.output.gen <- output_gen[is.na(output_gen)] <- 0
ranked.LS.gen <- rank_fh(complete.output.gen, 3)[,c(1,3,6)]
ranked.GLS.gen <- rank_fh(complete.output.gen, 4)[,c(1,4,6)]
ranked.PhyloGLM.gen <- rank_fh(complete.output.gen, 5)[,c(1,5,6)]
#ranked.GLM <- rank_fh(output, 5)[,c(1,5,7)]
#ranked.Log <- rank_fh(output, 6)[,c(1,6,7)]
colnames(ranked.LS.gen)[3] <- 'rankLS'
colnames(ranked.GLS.gen)[3] <- 'rankGLS'
colnames(ranked.PhyloGLM.gen)[3] <- 'rankPhyloGLM'
#colnames(ranked.Log)[3] <- 'rankLog'


# # basic plotting
# if(F){
# 	pairs(~ scoreLS + scoreGLS + scorePhyloGLM, data=output_gen)
# 	
# 	par(mfrow=c(3,1))
# 	hist(output_gen$scoreLS, breaks=60)
# 	hist(output_gen$scoreGLS, breaks=60)
# 	hist(output_gen$scorePhyloGLM, breaks=60)
# }

# #Huan's fancy ploting
# PhyloGLM20<-ranked.PhyloGLM.gen %>% filter(rankPhyloGLM == 20) %>% select(scorePhyloGLM)
# PhyloGLM_20<-ranked.PhyloGLM.gen %>% filter(rankPhyloGLM == -20) %>% select(scorePhyloGLM)
# PhyloGLM10<-ranked.PhyloGLM.gen %>% filter(rankPhyloGLM == 10) %>% select(scorePhyloGLM)
# PhyloGLM_10<-ranked.PhyloGLM.gen %>% filter(rankPhyloGLM == -10) %>% select(scorePhyloGLM)
# RRDR_gen <- merge(mutation_kmerPatterns, ranked.PhyloGLM.gen, by = "pattern", all.x = T)
# RRDR_gen <- unique(subset(RRDR_gen,select=c(1,3,4)))
# #ranked.PhyloGLM.gen %>% filter(scorePhyloGLM == 0) %>% tally
# score_hist <- function(k) {ggplot(ranked.PhyloGLM.gen,aes(scorePhyloGLM)) +
#     geom_histogram(data = RRDR_gen, fill = "green", alpha = 0.3) +
#     geom_histogram(data = ranked.PhyloGLM.gen, color = 'green',fill = 'white',alpha = 0.01) +
#     coord_trans(y='log1p') +
#     geom_vline(aes(xintercept=PhyloGLM20$scorePhyloGLM),linetype = 2,size = 0.5, color = 'green') +
#     geom_vline(aes(xintercept=PhyloGLM10$scorePhyloGLM),linetype = 1,size = 0.5, color = 'green') +
#     geom_vline(aes(xintercept=PhyloGLM_20$scorePhyloGLM),linetype = 2,size = 0.5, color = 'green') +
#     geom_vline(aes(xintercept=PhyloGLM_10$scorePhyloGLM),linetype = 1,size = 0.5, color = 'green')+
#     ylab('Number of patterns') +
#     xlab('Association score') +
#     annotate("text",  x=Inf, y = Inf, label = paste0("k=",k), vjust=2, hjust=2)
# } 
# score_hist(k)
##############################################
# score based on culture determined phenotype###
X <- t(t(trait_cul * 1))
xx <- cbind(ones, X)
threshold <- 1

pattern <- array(NA, c(nrow(Y), 1))
output_cul <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scorePhyloGLM=0)

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
    
    ##################
    # PhyloGLM
    fit <- phyloglm_aaf(y ~ X, phy=phy)
    output_cul$scorePhyloGLM[i] <-fit$zB[2]
  }
}


#If NA
complete.output.cul <- output_cul[is.na(output_cul)] <- 0
ranked.LS.cul <- rank_fh(complete.output.cul, 3)[,c(1,3,6)]
ranked.GLS.cul <- rank_fh(complete.output.cul, 4)[,c(1,4,6)]
ranked.PhyloGLM.cul <- rank_fh(complete.output.cul, 5)[,c(1,5,6)]
#add column names
colnames(ranked.LS.cul)[3] <- 'rankLS'
colnames(ranked.GLS.cul)[3] <- 'rankGLS'
colnames(ranked.PhyloGLM.cul)[3] <- 'rankPhyloGLM'

#Save the scores
kmerPattern.stats <- read.table(file = paste0("Cohen/TBcohen_k",k,"_n1_p1_kmerPattern.stats"), colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
total_score.gen <- merge(ranked.LS.gen, kmerPattern.stats, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.GLS.gen, total_score.gen, by = "pattern", all.x = T)
total_score.gen <- merge(ranked.PhyloGLM.gen, total_score.gen, by = "pattern", all.x = T)
#total_score.gen <- merge(ranked.Log, total_score.gen, by = "pattern", all.x = T)
write.csv(total_score.gen, paste0("TBcohen_k",k,"_n1_p1_scores_gen.csv"), row.names = FALSE)

total_score.cul <- merge(ranked.LS.cul, kmerPattern.stats, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.GLS.cul, total_score.cul, by = "pattern", all.x = T)
total_score.cul <- merge(ranked.PhyloGLM.cul, total_score.cul, by = "pattern", all.x = T)
#total_score.cul <- merge(ranked.Log, total_score.cul, by = "pattern", all.x = T)
write.csv(total_score.cul, paste0("TBcohen_k",k,"_n1_p1_scores_cul.csv"), row.names = FALSE)

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

#Get scores for true positive patterns
mutation_pattern_score.gen <- merge(mutation_kmerPatterns, total_score.gen, by = "pattern", all.x = T)
mutation_pattern_score.cul <- merge(mutation_kmerPatterns, total_score.cul, by = "pattern", all.x = T)
write.csv(mutation_pattern_score.gen, paste0("TBcohen_k",k,"_n1_p1_scores_RRDR_gen.csv", sep = ""), row.names = FALSE)
write.csv(mutation_pattern_score.cul, paste("TBcohen_k",k,"_n1_p1_scores_RRDR_cul.csv", sep = ""), row.names = FALSE)

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
