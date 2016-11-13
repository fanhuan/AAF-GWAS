# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# install.packages(c("R.oo","compoisson","R.methodsS3"))
# install.packages("/Users/arives/Downloads/phylosim_2.1.1.tar.gz")

library(readr)
library(phylosim)
library(ape)
library(Biostrings)
library(phylolm)

#rank_fh is a function that takes a dataframe and add a rank column according to one column (given the column number)
#such that the numbers above 0 are ranked from high to low as 1,2,3,4 while numbers below 0 are ranked from
#low to high as -1, -2, -3, -4. 0 is ranked as 0.
rank_fh <- function(df, number) {
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

setwd("~/Box Sync/Traveller Box/Huan GWAS/rpo_GY84_simulations_folder")

k <- 9

# Input tree
# cat("(((((MAF_11821_03:0.412,(MAF_GM_0981:0.582,MTB_95_0545:0.550):0.184):0.212,((MTB_K21:0.242,(MTB_K67:0.280,MTB_K93:0.234):0.096):0.108,(MTB_T17:0.624,MTB_T92:0.584):0.150):0.114):0.096,(((((MTB_00_1695:0.510,MTB_T67:0.518):0.182,MTB_T85:0.278):0.154,(MTB_98_1833:0.392,MTB_M4100A:0.072):0.142):0.134,(MTB_4783_04:0.574,MTB_GM_1503:0.596):0.166):0.118,MTB_91_0079:0.576):0.094):0.074,(MTB_K37:0.518,MTB_K49:0.288):0.070):0.138,MTB_H37Rv:0.256);", 
# file = "TBSimulation.nwk")

phy <- read.tree(file = "TBSimulation.nwk")
p <- Ntip(phy) #number of tips/species 

# sort(phy$tip.label)
# This matches the species list from kmercount
# MAF_11821_03
# MAF_GM_0981
# MTB_00_1695
# MTB_4783_04
# MTB_91_0079
# MTB_95_0545
# MTB_98_1833
# MTB_GM_1503
# MTB_H37Rv
# MTB_K21
# MTB_K37
# MTB_K49
# MTB_K67
# MTB_K93
# MTB_M4100A
# MTB_T17
# MTB_T67
# MTB_T85
# MTB_T92

Vphy <- vcv(phy)[sort(phy$tip.label), sort(phy$tip.label)]

C <- as.matrix(Vphy)
C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier
iC <- solve(C) #inverse of C
ones <- array(1, c(p, 1)) # array of ones for the intercept

for (repn in 50) {

	sim <- readRDS(file = paste("sim_", repn, ".rds", sep = ""))
	alignment.names <- names(sim$alignment[1, ])
	index <- (1:1172)[alignment.names == 450 & !is.na(alignment.names)]

	#Load in the shared kmer table
	phylokmer <- read.delim(paste("rpoB_GY84_", repn, "/phylokmer.dat", sep = ""), header = FALSE)
	colnames(phylokmer) <- c("kmer", sort(phy$tip.label))
	#Grep the patterns for kmers invloving S450 (kmer_list)
	
	#generate kmers including S450 (we do need to run aaf_phylosim.py first)
	kmer <- sim$alignment[, (index - k/3 + 1):(index + k/3 - 1)]
	kmer_list <- NULL #This does not consider possible deletions ('NA' in some codon)
	for (x in 1:p) {
		kmers <- paste(kmer[x, ], collapse = "")
		for (j in 1:(nchar(kmers) - k + 1)) {
			kmer_list <- c(kmer_list, substr(kmers, j, j + k - 1))
		}
	}
	kmer_df <- data.frame(kmer_list)
	kmer_df <- rename(kmer_df, c(kmer_list = "kmer"))
	kmer_df$position <- rep(c((nchar(kmers) - k + 1):1), p)
	#one kmer won't have different positions, too short
	kmer_df <- unique(kmer_df)
	
	#add their reverse compliment conterparts in (kmer_count only uses whoever that is alphabetically
	#first, the original kmer or it's rc.)
	rc <- NULL
	for (i in 1:nrow(kmer_df)) {
		kmer_df <- rbind(kmer_df, data.frame(kmer = as.character(reverseComplement(DNAString(kmer_df$kmer[i]))), position = -kmer_df$position[i]))
	}
	kmer_df <- unique(kmer_df)
	kmer_list_rc <- kmer_df$kmer
	S450_kmers <- phylokmer[phylokmer$kmer %in% kmer_list_rc, ]
	#Note that less than half of the kmers in kmer_list_rc ends up in S450_kmers because
	#1. only the original OR the rc kmer is in phylokmer
	#2. only kmers shared at least by two species are in phylokmer

	w <- read.csv(file = paste("rpoB_GY84_trait_", repn, ".csv", sep = ""))

	##############################################
	# scoring
	##############################################

	#y should be kmer pattern
	#X should be the trait
	#read in the kmer pattern as Y, and calculate score for each y.
	Y <- read.fwf(file = paste("rpoB_GY84_", repn, "_kmerPattern.stats", sep = ""), widths = array(1, c(1, p)), header = F)
	colnames(Y) <- sort(phy$tip.label)

	#sort the dataframe by tip name, get the serine column with Trues and falses and convert it into 0/1(by *1 or +0)
	trait <- w$serine[order(w$tip)]
	names(trait) <- w$tip[order(w$tip)]

	X <- t(t(trait * 1))
	xx <- cbind(ones, X)
	threshold <- 3
	
	pattern <- array(NA, c(nrow(Y), 1))
	output <- data.frame(pattern, sumy=0, scoreLS=0, scoreGLS=0, scoreGLM=0, scoreLog=0)
	for (i in 1:nrow(Y)) {
	#for (i in 1:500) {
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
			MSE <- t(h) %*% iC %*% h/(p - 2)
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
			mu = mean(y)
			B.init = matrix(c(log(mu/(1-mu)),0.0001), ncol=1)
			show(c(i, 1, sumy, system.time(z.GLM <- binaryPGLM(y ~ X, phy = phy, s2 = 1, B.init = B.init))))
			if(z.GLM$convergeflag == "converged") {
				output$scoreGLM[i] <- z.GLM$B[2]/z.GLM$B.se[2]
			} else {
				show('not converged')
			}
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
			
			show(c(i, 2, sumy, system.time(z.Log <- phyloglm(y ~ X, phy=phy))))
			output$scoreLog[i] <- z.Log$coefficients[2]/z.Log$sd[2]
			
		}
	}
	output <- output[1:i,]
#	show(plot(output[abs(output[,3]) < 1000,3:6]))
#	show(plot(output[,4:6]))
	ranked.LS <- rank_fh(output, 3)[,c(1,3,7)]
	ranked.GLS <- rank_fh(output, 4)[,c(1,4,7)]
	ranked.GLM <- rank_fh(output, 5)[,c(1,5,7)]
	ranked.Log <- rank_fh(output, 6)[,c(1,6,7)]
	colnames(ranked.LS)[3] <- 'rankLS'
	colnames(ranked.GLS)[3] <- 'rankGLS'
	colnames(ranked.GLM)[3] <- 'rankGLM'
	colnames(ranked.Log)[3] <- 'rankLog'
	
	#read in the kmer pattern again, focus on the frequency this time.
	kmerPattern.stats <- read.table(file = paste("rpoB_GY84_", repn, "_kmerPattern.stats", sep = ""), colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
	total_score <- merge(ranked.LS, kmerPattern.stats, by = "pattern", all.x = T)
	total_score <- merge(ranked.GLS, total_score, by = "pattern", all.x = T)
	total_score <- merge(ranked.GLM, total_score, by = "pattern", all.x = T)
	total_score <- merge(ranked.Log, total_score, by = "pattern", all.x = T)
	write.csv(total_score, paste("rpoB_GY84_", as.character(repn), "_scores.csv", sep = ""), row.names = FALSE, )

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
	S450_kmers_score <- merge(S450_kmers, total_score, by = "pattern", all.x = T)
	#Note that some of the kmers does not have scores because only patterns with more than two 1 or 0 are scored.
	#Merge S450_kmer_score with kmer_df
	
	output_450 <- merge(S450_kmers_score, kmer_df, by = "kmer", all.x = T)
	output_450_light <- subset(output_450, select = c("kmer", "position", "pattern", "freq", "scoreLS", "rankLS", "scoreGLS", "rankGLS", "scoreGLM", "rankGLM", "scoreLog", "rankLog"))
	write.csv(output_450_light, paste("rpoB_GY84_", as.character(repn), "_450summary.csv", sep = ""), row.names = FALSE, )

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