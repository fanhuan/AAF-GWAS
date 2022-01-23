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

LogF <- function(X,y){
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

PhyLog <- function(X,y,phy.extend){
    fit <- phyloglm(y ~ X, phy=phy.extend,method = "logistic_MPLE")
    scorePhyLog <- fit$coef[2]/fit$sd[2]
    return(scorePhyLog)
}

PLog <- function(X,y,phy.extend){
    fit <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=F)
    scorePLog <- fit$zB[2]
    return(scorePLog)
}

PLogF <- function(X,y,phy.extend){
    if(any(X == y) & any(X != y)) {
        fit_F <- try(phyloglm_aaf(y ~ X, phy=phy.extend, Firth=T, ultrametric=F), silent=T)
        if(is.null(attr(fit_F,'condition'))){
            scorePLogF <- fit_F$zB[2]
            # output$convergencePLogF[i]<- fit_F$convergence
        }else{
            scorePLogF <- 0
            # output$convergencePLogF[i] <- 100
        }
    }else{
        if(any(X == y)){
            scorePLogF <- Inf
        }else{
            scorePLogF <- -Inf
        }
    }
    return(scorePLogF)
}

PLogML <- function(X,y,phy.extend){
    if(any(X == y) & any(X != y)) {
        fit_ML <- try(phyloglm_ML(y ~ X, phy=phy.extend, ultrametric=T), silent=T)
        fit_ML0 <- try(phyloglm_ML(y ~ 1, phy=phy.extend, ultrametric=T), silent=T)
        if(is.null(attr(fit_ML,'condition')) & is.null(attr(fit_ML0,'condition'))){
            scorePLogML <- fit_ML$coef[2]/fit_ML$sd[2]
            LR <- fit_ML$logLik - fit_ML0$logLik
            if(LR > 0){
                scorePLogML_LR <- sign(fit_ML$coef[2]) * 2 * LR
            }else{
                scorePLogML_LR <- 0
            }
        }else{
            scorePLogML <- 0
            scorePLogML_LR <- 0
        }
    }else{
        if(any(X == y)){
            scorePLogML <- Inf
            scorePLogML_LR <- Inf
        }else{
            scorePLogML <- -Inf
            scorePLogML_LR <- -Inf
        }
    }
    #R functions don't return multiple objects in the strict sense.
    return(c(scorePLogML,scorePLogML_LR))
}


#Parameters
k <- 9
prefix_base <- 'phylosim_sp100d01_'
repn <- 115
threshold <- 3 #While in simulation this threshold is 4 therefore 3 does not exclude further
output.scores <- NULL
for (irep in 67:(repn-1)) {
    prefix <- paste0(prefix_base,irep)
    phy <- read.tree(file = paste0('rpoBsimulation/',prefix,".tre"))
    #Test whether there is zero length branches. If yes, collapse them.
    if (min(phy$edge.length) == 0) {
        min_edge <- min(phy$edge.length[which(phy$edge.length>0)])
        phy <- di2multi(phy,min_edge)
    }
    
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
    
    w <- read.csv(file = paste0('rpoBsimulation/',prefix,'/',prefix,"_trait.csv"))
    
    
    ##############################################
    # scoring
    ##############################################
    #y should be kmer pattern
    #X should be the trait
    #read in the kmer pattern as Y, and calculate score for each y.
    Y <- read.fwf(file = paste0('rpoBsimulation/',prefix,'/',prefix,"_kmerPattern.stats"), widths = array(1, c(1, p)), header = F)
    colnames(Y) <- sort(phy$tip.label)
    
    #sort the dataframe by tip name, get the serine column with Trues and falses and convert it into 0/1(by *1 or +0)
    # Resistant (not serine) = 1, Susceptible (serine) = 0
    trait <- w$serine[order(w$tip)]
    names(trait) <- w$tip[order(w$tip)]
    
    X <- t(t(trait * -1 +1))
    xx <- cbind(ones, X)
    
    pattern <- array(NA, c(nrow(Y), 1))
    output <- data.frame(pattern, sumy=0, scoreGLS=0, scoreLS=0, scorePhyLog=0, scoreLogF=0)
    
    for (i in 1:nrow(Y)) {
        y <- t(Y[i, ])
        output$pattern[i] <- paste(y, collapse = "")
        sumy <- sum(y) #sum of the pattern
        output$sumy[i] <- sumy
        if (sumy >= threshold & sumy <= (p - threshold)) {
       ##################
        # LS
        output$scoreLS[i] <- LS(xx,y)
        
        ##################
        # LogF
        output$scoreLogF[i] <- LogF(X,y)
        
        ##################
        # GLS
        output$scoreGLS[i]  <- GLS(xx,y,iC)
        
        ##################
        # PhyLog
        output$scorePhyLog[i] <- PhyLog(X,y,phy.extend)

        }
    }
    #fill NA with 0
    output$rankLS <- rank_fh(output$scoreLS)
    output$rankLogF <- rank_fh(output$scoreLogF)
    output$rankGLS <- rank_fh(output$scoreGLS)
    output$rankPhyLog <- rank_fh(output$scorePhyLog)
    
    #read in the kmer pattern again, focus on the frequency this time.
    kmerPattern.stats <- read.table(file = paste0('rpoBsimulation/',prefix,'/',prefix,"_kmerPattern.stats"), colClasses = c("character", "integer"), col.names = c("pattern", "freq"))
    total_score <- merge(output, kmerPattern.stats, by = "pattern", all.x = T)
    write.csv(total_score, paste0(prefix, "_scores_010Feb17.csv"), row.names = FALSE)
    
    #Grep the patterns for kmers invloving S450 (kmer_list)
    phylokmer <- read.delim(paste0('rpoBsimulation/',prefix, "/phylokmer.dat"), header=FALSE)
    intree <- read.tree(paste0('rpoBsimulation/',prefix,'.tre'))
    colnames(phylokmer)<-c("kmer",sort(intree$tip.label))
    
    #generate kmers including S450 (we do need to run aaf_phylosim.py first)
    kmer_df <- read.table(paste0('~/Dropbox/Research/GWAS/rpoBsimulation/',prefix,'/S450_',irep,'.kmer'),col.names = 'kmer')
    #add their reverse compliment conterparts in (kmer_count only uses whoever that is alphabetically
    #first, the original kmer or it's rc.)
    rc<-NULL
    for (i in 1:nrow(kmer_df)) {
        kmer_df<-rbind(kmer_df,data.frame(kmer=as.character(
            reverseComplement(DNAString(kmer_df$kmer[i])))))
    }
    kmer_list_rc<-kmer_df$kmer
    S450_kmers <-phylokmer[phylokmer$kmer %in% kmer_list_rc,]
    
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
    output_450_light <- subset(S450_kmers_score, select = c("kmer", "pattern", "sumy", "scoreGLS", "rankGLS", "scoreLS", "rankLS", "scorePhyLog", "rankPhyLog", "scoreLogF", "rankLogF"))
    write.csv(output_450_light, paste0('rpoBsimulation/',prefix,'/',prefix, "_450summary_10Feb17.csv"), row.names = FALSE)
    # reshape
    output_450_plot <- unique(output_450_light[,2:ncol(output_450_light)])
    output_450_plot <- output_450_plot[complete.cases(output_450_plot),]
    output_450_plot_long <- output_450_plot %>% reshape(varying=3:ncol(output_450_plot),
                                                        timevar="Model",
                                                        times = c("GLS","LS","PhyLog","LogF"),
                                                        v.names = c("Score","Rank"),
                                                        direction = "long",
                                                        new.row.names=1:10000) %>%
        mutate(Score1=Rank, Rank=Score, Score=Score1,Model = as.factor(Model),Rank_450=NA) %>%
        dplyr::select(-c(Score1,id)) %>% mutate(type = ifelse(is.na(Score), 'notappilicable',
                                                              ifelse(Score > 0, 'Positive',
                                                                     ifelse(Score == 0, 0, 'Negative'))))
    
    output_450_rank <- output_450_plot_long[0,]
    for (md in levels(output_450_plot_long$Model)) {
        subset <- output_450_plot_long %>% filter(Model == md)
        subset$Rank_450 <- rank_fh(subset$Score)
        output_450_rank <- rbind(output_450_rank,subset)
    }
    show(irep)
    output.scores <- rbind(output.scores, cbind(array(irep, dim=dim(output_450_rank)[1]),output_450_rank[,2:7]))
    
}
names(output.scores)[1] <- 'irep'

write.csv(file = 'rpoBsimulation_output_scores_10Feb17.csv', output.scores)

d <- output.scores[is.element(output.scores$type, c('pos','neg')),] 
d <- d[!is.infinite(d$Score),]
summary(d)

pdf(file=paste0('rpoBsimulation_rank_hist_pos_10Feb17.pdf'), width=10, height=10)
par(mfrow=c(3,2))
for(i in levels(d$Model)){
    dd <- d[d$type == 'pos' & d$Rank_450 == 1 & d$Model == i,]
    hist(dd$Rank, main = i, breaks = 0.5 + 0:max(dd$Rank))
}
dev.off()

pdf(file=paste0('rpoBsimulation_rank_hist_neg_10Feb17.pdf'), width=10, height=10)
par(mfrow=c(3,2))
for(i in levels(d$Model)){
    dd <- d[d$type == 'neg' & d$Rank_450 == 1 & d$Model == i,]
    hist(dd$Rank, main = i, breaks = 0.5 + 0:max(dd$Rank))
}
dev.off()

ranks <- 5
x <- data.frame(MeanRank = array(0,dim=length(levels(d$Model))*2*ranks), Model=levels(d$Model), Rank_450=0, type=c('pos','neg'))
i <- 0
for(iModel in levels(d$Model)) for(itype in c('pos','neg')) for(irank in 1:ranks) {
    i <- i + 1
    dd <- d[d$Model == iModel & d$type == itype & d$Rank_450 == irank,]
    x$MeanRank[i] <- mean(dd$Rank, na.rm=T)
    x$Model[i] <- iModel
    x$Rank_450[i] <- irank
    x$type[i] <- itype
}

pdf(file=paste0('rpoBsimulation_MeanRank_10Feb17.pdf'), width=10, height=10)
x %>% filter(type == 'pos'|type == 'neg')%>% ggplot(aes(Rank_450, abs(MeanRank))) +
    geom_point(aes(color=Model)) +
    geom_line(aes(color=Model)) +
    scale_y_log10(aes(abs(MeanRank))) +
    facet_grid(type~.)
dev.off()
