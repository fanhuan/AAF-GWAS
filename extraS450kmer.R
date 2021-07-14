library(ape)
library(Biostrings)
library(plyr)
library(dplyr)
library(ggplot2)

rank_fh <- function(serial) {
  position <- seq(length(serial))
  #ranking <- is.na(serial)
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

sp <- 100
depth <- 0.1
output.scores <- NULL
for (irep in 0:66) {
  prefix = paste0("phylosim_","sp",sp,"d0",strsplit(as.character(depth),split = '.',fixed = TRUE)[[1]][2],'_',irep)
  
  phylokmer <- read.delim(paste0('rpoBsimulation/',prefix, "/phylokmer.dat"), header=FALSE)
  intree <- read.tree(paste0('rpoBsimulation/',prefix,'.tre'))
  colnames(phylokmer)<-c("kmer",sort(intree$tip.label))
  #Grep the patterns for kmers invloving S450 (kmer_list)
  
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
  total_score <- read.csv(paste0('rpoBsimulation/',prefix,'/',prefix, "_scores_15Nov16.csv"),colClasses=c("character",rep("numeric",14)))
  S450_kmers_score <- merge(S450_kmers, total_score, by = "pattern", all.x = T)
  #Note that some of the kmers does not have scores because only patterns with more than two 1 or 0 are scored.
  #Merge S450_kmer_score with kmer_df
  output_450_light <- subset(S450_kmers_score, select = c("kmer", "pattern", "sumy", "scoreGLS", "rankGLS", "scoreLS", "rankLS", "scorePLogML_LR", "rankPLogML_LR", "scorePLogF", "rankPLogF", "scorePLogML", "rankPLogML"))
  write.csv(output_450_light, paste0('rpoBsimulation/',prefix,'/',prefix, "_450summary_15Nov16_7Feb17.csv"), row.names = FALSE)
  # reshape
  output_450_plot <- unique(output_450_light[,2:ncol(output_450_light)])
  output_450_plot <- output_450_plot[complete.cases(output_450_plot),]
  output_450_plot_long <- output_450_plot %>% reshape(varying=3:ncol(output_450_plot),
                                                      timevar="Model",
                                                      times = c("GLS","LS","PLogML_LR","PLogF","PLogML"),
                                                      v.names = c("Score","Rank"),
                                                      direction = "long",
                                                      new.row.names=1:10000) %>%
    mutate(Score1=Rank, Rank=Score, Score=Score1,Model = as.factor(Model),Rank_450=NA) %>%
    dplyr::select(-c(Score1,id)) %>% mutate(type = ifelse(is.na(Score), 'notappilicable',
                                                          ifelse(Score > 0, 'pos',
                                                                 ifelse(Score == 0, 0, 'neg'))))
  
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

write.csv(file = 'rpoBsimulation_output_scores_15Nov16_7Feb17.csv', output.scores)

d <- output.scores[is.element(output.scores$type, c('pos','neg')),] 
d <- d[!is.infinite(d$Score),]
summary(d)

pdf(file=paste0('rpoBsimulation_rank_hist_pos_15Nov16_7Feb17.pdf'), width=10, height=10)
par(mfrow=c(3,2))
for(i in levels(d$Model)){
  dd <- d[d$type == 'pos' & d$Rank_450 == 1 & d$Model == i,]
  hist(dd$Rank, main = i, breaks = 0.5 + 0:max(dd$Rank))
}
dev.off()

pdf(file=paste0('rpoBsimulation_rank_hist_neg_15Nov16_7Feb17.pdf'), width=10, height=10)
par(mfrow=c(3,2))
for(i in levels(d$Model)){
  dd <- d[d$type == 'neg' & d$Rank_450 == 1 & d$Model == i,]
  hist(dd$Rank, main = i, breaks = 0.5 + 0:max(dd$Rank))
}
dev.off()

x <- data.frame(MeanRank = array(0,dim=5*2*4), Model=levels(d$Model), Rank_450=0, type=c('pos','neg'))
i <- 0
for(iModel in levels(d$Model)) for(itype in c('pos','neg')) for(irank in 1:4) {
  i <- i + 1
  dd <- d[d$Model == iModel & d$type == itype & d$Rank_450 == irank,]
  x$MeanRank[i] <- mean(dd$Rank, na.rm=T)
  x$Model[i] <- iModel
  x$Rank_450[i] <- irank
  x$type[i] <- itype
}

pdf(file=paste0('rpoBsimulation_MeanRank_15Nov16_7Feb17.pdf'), width=10, height=10)
x %>% filter(type == 'pos'|type == 'neg')%>% ggplot(aes(Rank_450, abs(MeanRank))) +
  geom_point(aes(color=Model)) +
  geom_line(aes(color=Model)) +
  scale_y_log10(aes(abs(MeanRank))) +
  facet_grid(type~.)
dev.off()
