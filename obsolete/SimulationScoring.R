# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# install.packages(c("R.oo","compoisson","R.methodsS3"))
# install.packages("/Users/arives/Downloads/phylosim_2.1.1.tar.gz")

library(readr)
library(ape)
library(Biostrings)
library(phylolm)
source("~/Dropbox/Research/GWAS/phyloglm_aaf.R")
source("~/Dropbox/Research/GWAS/phyloglm_ML.R")
library(plyr)
library(dplyr)
library(logistf)
library(ggplot2)

#User defined functions
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


#Parameters
k <- 9
prefix_base <- 'phylosim_sp?d?_'
repn <- ?
threshold <- 3

for (iter in 0:(repn-1))
	prefix <- paste0(prefix_base,iter)
	phy <- read.tree(file = paste0(prefix,"/",prefix,".tre"))
	#Test whether there is zero length branches. If yes, collapse them.
	if (min(phy$edge.length) == 0) {
		min_edge <- min(phy$edge.length[which(phy$edge.length>0)])
		phy <- di2multi(phy,min_edge)
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
	C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier
	iC <- solve(C) #inverse of C
	ones <- array(1, c(p, 1)) # array of ones for the intercept

	w <- read.csv(file = paste0(prefix,"_trait.csv"))

	##############################################
	# scoring
	##############################################

	#y should be kmer pattern
	#X should be the trait
	#read in the kmer pattern as Y, and calculate score for each y.
	Y <- read.fwf(file = paste0(prefix,"_kmerPattern.stats"), widths = array(1, c(1, p)), header = F)
	colnames(Y) <- sort(phy$tip.label)

	#sort the dataframe by tip name, get the serine column with Trues and falses and convert it into 0/1(by *1 or +0)
	# Resistant (not serine) = 1, Susceptible (serine) = 0
	trait <- w$serine[order(w$tip)]
	names(trait) <- w$tip[order(w$tip)]

	X <- t(t(trait * -1 +1))
	xx <- cbind(ones, X)
	
	pattern <- array(NA, c(nrow(Y), 1))
	output <- data.frame(pattern, sumy=0, scoreLS=0, scoreLog=0, scoreGLS=0, scorePLog=0, scorePLogF=0, scorePLogML=0)
	nModel <- ncol(output)-2
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
			scoreLS <- b[2]/bSE
			output$scoreLS[i] <- scoreLS
			
			##################
			# logistf
			fit_Log <- logistf(y ~ X)
			scoreLog <- sign(fit_Log$coef[2]) * qchisq(fit_Log$prob[2], df=1, lower.tail=F)
			output$scoreLog[i] <- scoreLog
			
			##################
			# GLS
			XiCX <- t(xx) %*% iC %*% xx
			XiCY <- t(xx) %*% iC %*% y
			b <- solve(XiCX, XiCY)
			h <- y - (xx %*% b)
			MSE <- t(h) %*% iC %*% h/(p - 2)
			iXiCX <- solve(XiCX)
			bSE <- (MSE * iXiCX[2, 2])^0.5
			scoreGLS <- b[2]/bSE
			output$scoreGLS[i] <- scoreGLS
			
			##################
			# PLog
			fit <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=F)
			scorePLog <- fit$zB[2]
			output$scorePLog[i] <- scorePLog
			
			##################
			# PLogF
			fit_F <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=T)
			scorePLogF <- fit_F$zB[2]
			output$scorePLogF[i] <- scorePLogF
			
			##################
			# PLogML
			fit_ML <- try(phyloglm_ML(y ~ X, phy=phy.extend), silent=T)
			fit_ML0 <- try(phyloglm_ML(y ~ 1, phy=phy.extend), silent=T)
			if(is.null(attr(fit_ML,'condition')) & is.null(attr(fit_ML0,'condition'))){
				scorePLogML <- sign(fit_ML$coef[2]) * 2*(fit_ML$logLik - fit_ML0$logLik)
			}else{
				scorePLogML <- NA
			}
			output$scorePLogML[i] <- scorePLogML
			
			show(c(i, scoreLS, scoreLog, scoreGLS, scorePLog, scorePLogF, scorePLogML))
			
		}
	}

	#Put ranks for each model

	output$rankLS <- rank_fh(output$scoreLS)
	output$rankLog <- rank_fh(output$scoreLog)
	output$rankGLS <- rank_fh(output$scoreGLS)
	output$rankPLog <- rank_fh(output$scorePLog)
	output$rankPLogF <- rank_fh(output$scorePLogF)
	output$rankPLogML <- rank_fh(output$scorePLogML)

	#read in the kmer pattern again, focus on the frequency this time.
	kmerPattern.stats <- read.table(file = paste0(prefix, "_kmerPattern.stats"), colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
	total_score <- merge(output, kmerPattern.stats, by = "pattern", all.x = T)

	#Load kmers including S450 (we do need to run aaf_phylosim.py first)
	S450_kmers <- read.csv(paste0(prefix,"_S450.kmer"))
	#Note that less than half of the kmers in kmer_list_rc ends up in S450_kmers because
	#1. only the original OR the rc kmer is in phylokmer
	#2. only kmers shared at least by two species are in phylokmer

	#Get scores for S450 containing kmers.
	kft <- 1 #kmer frequency threshold
	pattern450 <- array(0, nrow(S450_kmers))
	for (ii in 1:nrow(S450_kmers)) {
		pattern450_j <- array(1, ncol(S450_kmers) - 1)
		for (jj in 2:ncol(S450_kmers)) {
			if (S450_kmers[ii, ][jj] < kft) {
				pattern450_j[jj - 1] <- 0
			}
		}
		pattern450[ii] <- paste(pattern450_j, collapse = "")
	}
	S450_kmers$pattern <- pattern450

	total_score <- merge(S450_kmers, total_score, by = "pattern", all.y = T)
	total_score$correct <- !is.na(total_score$kmer)
	total_score <- total_score[,-grep("^t",colnames(total_score))]
	total_score <- unique(subset(total_score,select=-kmer))
	write.csv(total_score, paste0(prefix, "_scores.csv"), row.names = FALSE)

	#Note that some of the kmers does not have scores because only patterns with more than two 1 or 0 are scored.
	#Merge S450_kmer_score with kmer_df

	output_450 <- total_score %>% filter(correct==T)
	outputplot <- total_score[abs(total_score$scoreLS) < 1000,]
	outputplot <- outputplot[order(abs(outputplot$scoreGLS)),]
	col450 <- is.element(outputplot$pattern, output_450$pattern)
	show(plot(outputplot[,3:ncol(outputPlot)], col=(1+col450), pch=(20 - col450), cex=(.5 + 1*col450)))


	#reshpe to long form
	output_450_long <- gather(output_450, Model, Score, seq(2,2*nModel,by=2))
	output_450_long1 <-gather(output_450, Model, Rank, seq(3,2*nModel+1,by=2))
	output_450_long$Rank <-output_450_long1$Rank
	output_450_long<-output_450_long[, -grep("^rank", colnames(output_450_long))]
	output_450_long<- output_450_long %>% mutate(Model = as.factor(sub("score", "", Model)),Rank_450=NA)%>%
	mutate(type = ifelse(is.na(Score), 'notappilicable',
	ifelse(Score > 0, 'pos',
	ifelse(Score == 0, 0, 'neg'))))
	#Copy the dataFRAME without copying the numbers
	output_450_rank <- output_450_long[0,]
	for (md in levels(output_450_long$Model)) {
		subset <- output_450_long %>% filter(Model == md)
		subset$Rank_450 <- rank_fh(subset$Score)
		output_450_rank <- rbind(output_450_rank,subset)
	}

	rank_plot<-output_450_rank %>% filter(type == 'pos'|type == 'neg')%>% ggplot(aes(Rank_450, abs(Rank))) +
	geom_point(aes(color=Model)) +
	geom_line(aes(color=Model)) +
	scale_y_log10(aes(abs(Rank))) +
	facet_grid(type~.) +
	annotate("text",  x=Inf, y = Inf, label = paste0(prefix), vjust=2, hjust=2)

	pdf(rank_plot, paste0(prefix, "_ranks.pdf"), width=10, height=10)
	dev.off()
	write.csv(output_450_rank,paste0(prefix, "_450summary.csv"), row.names = FALSE)
}