library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(hash)
library(parallel)
# Figure 1: 100 rpoBsimulation for all 11 scoring methods.
output.scores <- read_csv("rpoBsimulation/sp100d01/S450_summary/rpoBsimulation_output_scores_1Dec21.csv")
output.scores <- read_csv("rpoBsimulation/sp100d01/S450_summary/rpoBsimulation_output_scores_thr10_1Dec21.csv")
# Rename the models
old_model <- c('LS','GLS','EGLS','LogF','PhyLog','PLogML_LR_a',
               'PLog','PLogML','PLogML_LR','PLogML_LR_p','phyloglmA')
new_model <- c('LS','GLS','EGLS','LogF','PhylogF_EstA','PhylogML_EstA',
               'Phylog','PhylogML','PhylogML_LR','PhylogML_LR_p','PhylogF')
final_model <- c('LS','GLS','EGLS','LogF','Phylog_EstA','Phylog_EstA_LR',
                'Remove','Phylog','Phylog_LR','Phylog_pLR','PhylogIG')
h <- hash(final_model,old_model)
output.scores.1 <- output.scores %>% 
  left_join(stack(h), by = c('Model' = 'values')) %>% 
  rename(New_Model = ind)

# Remove 0s 
d <- output.scores.1[is.element(output.scores.1$Type, c('Positive','Negative')),] 
# Remove infinites
d <- d[!is.infinite(d$Score),]
# Remove removed models and inbalanced data
d <- d %>% filter(New_Model != 'Remove', sumy > 4, sumy < 96, irep < 101)
summary(d)
d$Rank <- as.integer(d$Rank)
d$model_f = factor(d$New_Model, levels=final_model)
##############################################################
###############PLOTTING#######################################
# Only plot positives.
# All ranks and top ranks.
################################################################
# Using ggplot so the x axis are comparable
pdf(file=paste0('Fig1.rpoBsimulation_positive_allrank_toprank_hist_15Feb22.pdf'), width=15, height=5)
p1 = d %>% filter(Type == 'Positive') %>% ggplot(mapping = aes(x = Rank)) + 
  geom_histogram(position="identity") + 
  xlab("All Ranks") + 
  facet_grid(~model_f)
p2 = d %>% filter(Rank_450 == 1, Type == 'Positive') %>% ggplot(mapping = aes(x = Rank)) + 
  geom_histogram(position="identity") + 
  xlab("Top Ranks") +
  facet_grid(~model_f)
grid.arrange(p2,p1,nrow = 2)
dev.off()

pdf(file=paste0('FigS1.rpoBsimulation_negative_allrank_toprank_hist_15Feb.pdf'), width=15, height=5)
p1 = d %>% filter(Type == 'Negative') %>% ggplot(mapping = aes(x = Rank)) + 
  geom_histogram(position="identity") + 
  xlab("All Ranks") +
  facet_grid(~model_f)
p2 = d %>% filter(Rank_450 == 1, Type == 'Negative') %>% ggplot(mapping = aes(x = Rank)) + 
  geom_histogram(position="identity") + 
  xlab("Top Ranks") +
  facet_grid(~model_f)
grid.arrange(p2,p1,nrow = 2)
dev.off()

# We can see that:
# 1. Sensitivity: ability to pick up true positives: EGLS and GLS performs the best in terms of the height of the first bar (higher percentage in the high ranks)
# 2. Specificity: ability to not pick up true negatives: PLogML_LR, the threshold is much smaller (<400 vs ~600 for others) to include all true positives. 

# Mean rank (for the top5-ranked patterns)
ranks <- 10 # top n, *2 is for negative and positve. The highest Rank_450 is 13(at most 13 pos/neg patterns in our 115 simulations )
x <- data.frame(MeanRank = array(0,dim=length(unique(d$model_f))*2*ranks),
                MedianRank = array(0,dim=length(unique(d$model_f))*2*ranks),
                Model=unique(d$model_f), Rank_450=0, Type=c('Positive','Negative'))
i <- 0
for(iModel in unique(d$model_f)) for(itype in c('Positive','Negative')) for(irank in 1:ranks) {
  i <- i + 1
  dd <- d[d$model_f == iModel & d$Type == itype & d$Rank_450 == irank,]
  x$MeanRank[i] <- mean(dd$Rank, na.rm=T)
  x$MedianRank[i] <- median(dd$Rank, na.rm=T)
  x$Model[i] <- iModel
  x$Rank_450[i] <- irank
  x$Type[i] <- itype
}
date = '15Feb22'
pdf(file=paste0('Fig2.rpoBsimulation_MeanRank', ranks, '_', date,'.pdf'), width=10, height=10)
x %>% filter(Type == 'Positive'|Type == 'Negative')%>% ggplot(aes(Rank_450, abs(MeanRank))) +
  geom_point(aes(color=Model)) +
  geom_line(aes(color=Model)) +
  scale_y_log10(aes(abs(MeanRank))) +
  facet_grid(Type~.)
dev.off()
pdf(file=paste0('Fig2.rpoBsimulation_MedianRank', ranks, '_', date,'.pdf'), width=10, height=10)
x %>% filter(Type == 'Positive'|Type == 'Negative')%>% ggplot(aes(Rank_450, abs(MedianRank))) +
  geom_point(aes(color=Model)) +
  geom_line(aes(color=Model)) +
  scale_y_log10(aes(abs(MedianRank))) +
  facet_grid(Type~.)
dev.off()

# This is very obvious that PLogML_LR is the best (stays the lowest). 
# But Tony thinks who is the lowest at Rank_450=1 is more important.
# But why it has less Rank450 than others? Does it give the same scores to different patterns?
irep113 <- output.scores %>% filter(irep == 113)
irep113 %>% group_by(Model) %>% tally
irep113 %>% filter(Rank_450 == 13)
irep113 %>% filter(Model == 'LogF') %>% arrange(Rank_450)
# Because not all patterns gets the same type! The ones around 0 can sometimes be positive and sometimes negative!
# I could look into why but the ones around 0 is trivial. 
# What I worry more is why in the mean plot the second rank dipped 
# 1. inbalanced kmer pattern (too many 1s or too many 0s) resulted in very bad Rank for Rank_450=1 (fixed by increasing threshold to 4).
# But we cannot just remove them. They are true positive even if the complexity is low.
# A more reasonable way is to remove simulations with inbalanced trait data.
# Currently, simulations with at least 3 serines and 3 non-serines are saved. 
# Let's see whether increase the threshold would help.
s <- phylosim_sp100d01_49_trait[phylosim_sp100d01_49_trait$serine == 'TRUE',]
dim(s)[1]
c <- rep(0,115)
for (i in 0:114){
  trait <- read_csv(paste0("rpoBsimulation/sp100d01/trait/phylosim_sp100d01_", i, 
                           "_trait.csv"))
  s <- trait[trait$serine == 'FALSE',]
  c[i+1] <- dim(s)[1]
}
hist(c)
# The lowest mutations are 4 (too litter info?) and the highest is 93 (also pretty
# in balanced). But is there a correlation of the imbalance of the trait data 
# versus the imbalance in the true positive k-mer patterns?

da <- read_csv("rpoBsimulation/sp100d01/S450_summary/rpoBsimulation_output_scores_1Dec21.csv")
da <- read_csv("rpoBsimulation/sp100d01/S450_summary/rpoBsimulation_output_scores_thr10_1Dec21.csv")
da %>% group_by(Model) %>% tally()
thr = 10
da_LS <- da %>% filter(Model == 'LS') %>% 
  mutate(Group = ifelse(sumy < thr, 'imbalanced', 
                        ifelse(sumy > 100-thr, 'imbalanced', 'balanced'))) %>% 
  group_by(irep, Group) %>% tally() %>% group_by(irep) %>% mutate(per_kmer = n/sum(n))
irep <- seq(0,114)
ids <- cbind(c, irep)
da <- merge(da, ids, by = 'irep')
da_LS <- merge(da_LS, ids, by = 'irep')
da_LS_im <- da_LS %>% filter(Group == 'imbalanced')
da_zero <- da %>% filter(Rank == 0) %>% group_by(Model) %>% tally()
da_one <- da %>% filter(Rank == 1) %>% group_by(Model) %>% tally() 
# 536 true positives in our 115 simulations, 313 was able to rank it 1st among all k-mers.
da_one_irep <- da %>% filter(Rank == 1) %>% group_by(Model, irep, Type) %>% tally()
length(unique(da_one_irep$irep))
# So among 115 simulations, there are 101 ended up having the first ranked k-mer being true positives.
da_one_type <- da_one_irep %>% group_by(Model,Type) %>% tally()
# If we don't distinguish positive/negative types, among those 101 ireps, PLogML_LR_a performed the best(100 ireps) followed by GLS, 
# PLogML_LRand PLogML_LR_p (ireps 99) 
da_one_type <- da_one_irep %>% group_by(Model,Type) %>% tally()
# But if we do, GLS performs the best in terms of the percentage of ireps with the highest positive scores being a true positive kmer pattern.
# Interestingly, the negative patterns peformed way better (6 methods had 95 ireps among 101)
# So now two things to check:
# 1. What is the biological meaning of negative k-mer patterns. Are their counterparts (1101 vs. 0010)
# instantly scores the same with a different sign?

# How do we utilize the negative info?

# 2. If we combine different methods, would it increase the number of ireps with the highest positve score being true (now highest is GLS 81)?

plot(da_LS_im$c, da_LS_im$per) %>% group_by(irModel) %>% tally
plot(da_LS)
# how to summarize imbalance?
hist(da$sumy)
# It seems like the true positive k-mer patterns tend to be distributed at the two ends of sumy.
da %>% group_by(Model) %>% tally() %>% arrange(n)
# Every truo positive k-mer pattern was scored using all methods. Let's only use rows under LS as a way to dereplicate.

d_n2 <- d %>% filter(Rank_450 < 3)
tmp <- d_n2 %>% filter(model_f == "Phylog_EstA_LR",Type == 'Negative')
tmp <- d %>% filter(irep == '20', Type == 'Negative')
tmp <- d %>% filter(irep == '20', Model == 'phyloglmA')

# The result for imbalanced k-mer patterns are less trustful. Therefore we will 
# filter out all imbalanced k-mers and re-rank. Where is the original data?
prefix_base <- 'phylosim_sp100d01_'
repn <- 115
date <- '1Dec21'
thr = 10 
model_list = c('LS','LogF','GLS','EGLS','PhyLog','PLogML_LR_a','PLog','phyloglmA','PLogML_all')
if ('PLogML_all' %in% model_list){
  index <- match('PLogML_all',model_list)
  model_name_list <- c(model_list[1:index-1],"PLogML", "PLogML_LR", "PLogML_LR_p", 
                       model_list[index+1:length(model_list)])[1:(length(model_list)+2)]
}else{
  model_name_list <- model_list            
}
data_dir <- "rpoBsimulation/sp100d01/"
prefix_base <- 'phylosim_sp100d01_'
output.scores <- NULL #saving all S450 scores
date <- '1Dec21'
for (irep in 0:114) { # count from 0
  print(irep)
  prefix <- paste0(prefix_base,irep)
  score <- as.data.frame(read_csv(paste0("rpoBsimulation/sp100d01/score/phylosim_sp100d01_", 
                           irep,"_scores_1Dec21.csv")))
  score.ba <- score %>% filter(sumy > thr, sumy < 100-thr)
  score.ba <- score.ba[ , !grepl( "rank" , names( score.ba ) ) ]
  reps_m <- seq(length(model_name_list))
  # x+3 because now there is an extra column of freq.
  r_m <- mclapply(reps_m, function(x){rank_fh(as.matrix(score.ba[,x+3]))},
                  mc.cores = 4)
  output_rank <- as.data.frame(do.call(cbind, r_m))
  names(output_rank)<-c(sapply('rank',paste0,model_name_list,USE.NAMES = F))
  output <- cbind(score.ba,output_rank)
  S450_kmer <- read.csv(paste0(data_dir, 'S450_kmer_pattern/', prefix, '_S450.kmer'),colClasses = c("character", "character"))
  S450_score <- merge(S450_kmer, output, by = "pattern")
  if (dim(S450_score)[1] > 0){
    # Summarize into one excel sheet like rpoBsimulation_output_scores_20Aug2021.csv
    output_450_plot <- S450_score %>% select(-kmer)
    output_450_plot <- unique(output_450_plot) 
    # rearrange the data set so score and rank for one method is together
    order_c <- rep(0,length(model_name_list) * 2)
    for (i in 1:length(model_name_list)){
      # j is the position in order_c
      j = i*2
      # position of the score
      order_c[j-1] = 3 + i
      # position of the rank
      order_c[j] = 3 + length(model_name_list) + i
    }
    output_450_plot <- output_450_plot[,c(1:3,order_c)]
    # reshape into long form
    output_450_plot_long <- output_450_plot %>% reshape(varying=4:ncol(output_450_plot),
                                                        timevar="Model",
                                                        times = model_name_list,
                                                        v.names = c("Rank","Score"),
                                                        direction = "long",
                                                        new.row.names=1:10000) %>%
      # for some reason the rank and score is swapped
      mutate(Score1=Rank, Rank=Score, Score=Score1,Model = as.factor(Model),Rank_450=NA) %>%
      dplyr::select(-c(Score1,id)) %>% mutate(Type = ifelse(is.na(Score), 'notappilicable',
                                                            ifelse(Score > 0, 'Positive',
                                                                   ifelse(Score == 0, 0, 'Negative'))))
    output_450_rank <- output_450_plot_long[0,]
    for (md in levels(output_450_plot_long$Model)) {
      subset <- output_450_plot_long %>% filter(Model == md)
      subset$Rank_450 <- rank_fh(subset$Score)
      output_450_rank <- rbind(output_450_rank,subset)
    }
    output.scores <- rbind(output.scores, cbind(array(irep, dim=dim(output_450_rank)[1]),output_450_rank[,2:ncol(output_450_rank)]))
  }else{
    print(paste('all S450 patterns are imbalanced from irep', irep))
  }
}
names(output.scores)[1] <- 'irep'
write.csv(file = paste0(data_dir, "S450_summary/rpoBsimulation_output_scores_thr", thr, '_',date, ".csv"), output.scores, row.names = F)

