library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(hash)
# Figure 1: 100 rpoBsimulation for all 11 scoring methods.
output.scores <- read_csv("rpoBsimulation/sp100d01/S450_summary/rpoBsimulation_output_scores_1Dec21.csv")
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
# 1. inbalanced data resulted in very bad Rank for Rank_450=1 (fixed by increasing threshold to 4).
d_n2 <- d %>% filter(Rank_450 < 3)
tmp <- d_n2 %>% filter(model_f == "Phylog_EstA_LR",Type == 'Negative')
tmp <- d %>% filter(irep == '20', Type == 'Negative')
tmp <- d %>% filter(irep == '20', Model == 'phyloglmA')

