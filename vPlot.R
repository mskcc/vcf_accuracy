#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript


library(Cairo)
library(data.table)
library(ggplot2)
library(reshape)
library(RColorBrewer)

mgrep <- function(patterns, x, exact=F, v=F, ...) {
  gpat <- patterns
  if (exact) {
    gpat = paste("^", patterns, "$", sep="")
  }
  if (!v) {
    out <- unlist(lapply(gpat, function(p,y=x) grep(p,y)))
  }
  else {
    out <- x[unlist(lapply(gpat, function(p,y=x) grep(p,y)))]
  }
  out
}

get.percentages = function(Missed_SNP, Novel_SNP, Correct_SNP_Genotype, Incorrect_SNP_Genotype, 
                           Missed_INDEL,  Novel_INDEL, Correct_INDEL_Genotype, Incorrect_INDEL_Genotype,  
                           Union_SNPs, Union_INDELs){
  
  if(Union_SNPs == 0){   
    PC_Missed_SNP = 0
    PC_Correct_SNP_Genotype = 0
    PC_Incorrect_SNP_Genotype = 0
  }
  
  if(Union_INDELs == 0){  
    PC_Missed_INDEL = 0
    PC_Correct_INDEL_Genotype = 0
    PC_Incorrect_INDEL_Genotype = 0
  }
  
  if(Union_SNPs == 0){PC_Novel_SNP = 0}
  if(Union_INDELs == 0){PC_Novel_INDEL = 0}
  
  if(Union_SNPs > 0){
    PC_Missed_SNP = (Missed_SNP / Union_SNPs) * 100
    PC_Correct_SNP_Genotype = (Correct_SNP_Genotype / Union_SNPs) * 100
    PC_Incorrect_SNP_Genotype = (Incorrect_SNP_Genotype / Union_SNPs) * 100
  }
  
  if(Union_INDELs > 0){
    PC_Missed_INDEL = (Missed_INDEL / Union_INDELs) * 100
    PC_Correct_INDEL_Genotype = (Correct_INDEL_Genotype / Union_INDELs) * 100
    PC_Incorrect_INDEL_Genotype =  (Incorrect_INDEL_Genotype / Union_INDELs) * 100
  }
  
  ###################################################################################################################
  
  if(Union_SNPs > 0){
    PC_Novel_SNP = (Novel_SNP / Union_SNPs) *100 }
  
  if(Union_INDELs > 0){
    PC_Novel_INDEL = (Novel_INDEL / Union_INDELs) *100 }
  
  ###################################################################################################################
  
  list(PC_Missed_SNP,
       PC_Novel_SNP,
       PC_Correct_SNP_Genotype,
       PC_Incorrect_SNP_Genotype,
       PC_Missed_INDEL,
       PC_Novel_INDEL,
       PC_Correct_INDEL_Genotype,
       PC_Incorrect_INDEL_Genotype) 
}

#####################################################################################################################
#####################################################################################################################

args = commandArgs(trailingOnly=T)
breakdown = args[1]
details = args[2]
prefix = args[3]
samples = fread(breakdown)
samples.details = fread(details)

samples[,c('Missed SNP', 
           'Novel SNP', 
           'Correct SNP Genotype', 
           'Incorrect SNP Genotype',
           'Missed INDEL', 
           'Novel INDEL',
           'Correct INDEL Genotype' ,
           'Incorrect INDEL Genotype'):=get.percentages(Missed_SNP, Novel_SNP, Correct_SNP_Genotype, Incorrect_SNP_Genotype, 
                                                        Missed_INDEL,  Novel_INDEL, Correct_INDEL_Genotype, Incorrect_INDEL_Genotype,  
                                                        Union_SNPs, Union_INDELs), by=1:nrow(samples)]

# NORMAL_TAGS <- c('NORMAL', 'Normal', '_N[0-9]*$')
# samples <- samples[-mgrep(NORMAL_TAGS, sample), ]

sample = melt(samples, id.vars=c('sample'), measure.vars=c('Missed SNP', 'Novel SNP', 'Correct SNP Genotype', 'Incorrect SNP Genotype',  
                                                           'Missed INDEL', 'Novel INDEL', 'Correct INDEL Genotype', 'Incorrect INDEL Genotype'))

breakdown.plot = ggplot(sample, aes(y=value, x =variable, fill =sample)) +  
  geom_bar(stat = 'identity', position = 'dodge')  +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=.6, size=12), 
        axis.text.y  = element_text(size=12), 
        axis.title   = element_text(size=15),
        legend.text  = element_text(size=12), 
        legend.title = element_text(size=15, face='bold'), 
        legend.position='none') +
  labs(x='', y='Percentage of Variants')

CairoPNG(paste(prefix,'_percentages.png', sep=''),units='px',w=800,h=600)
breakdown.plot
dev.off()

setnames(samples.details, c('Type', 'Coordinates', 'Class', 'Sample', 'DP', 'AD'))
samples.details$DP <- as.numeric(samples.details$DP)
samples.details$AD <- as.numeric(samples.details$AD)
samples.details[, AF := (AD / DP)]

detail = melt(samples.details, id.vars='Class', measure.vars='AF')

detail.plot = ggplot(detail, aes(y=value, x =Class)) +  
  geom_point(position = position_jitter(width=0.1)) +
  geom_boxplot(fill=NA, aes(color=Class), outlier.shape=NA)  +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=.6, size=12), 
        axis.text.y  = element_text(size=12), 
        axis.title   = element_text(size=15),
        legend.text  = element_text(size=12), 
        legend.title = element_text(size=15, face='bold'), 
        legend.position='none') +
  labs(y='Allelic Fraction')

CairoPNG(paste(prefix,'_allelicFractions.png', sep=''),units='px',w=800,h=600)
detail.plot
dev.off()

total.SNPs <- sum(samples$Union_SNPs)
total.Indels <- sum(samples$Union_INDELs)
concordant.SNPs <- sum(samples$Correct_SNP_Genotype)
discordant.SNPs <- total.SNPs - concordant.SNPs
concordant.Indels <- sum(samples$Correct_INDEL_Genotype)
discordant.Indels <- total.Indels - concordant.Indels

confusion <- matrix(c(concordant.SNPs, concordant.Indels, discordant.SNPs, discordant.Indels), nrow=2, byrow = T)

colnames(confusion) <- c('SNPs', 'Indels')
rownames(confusion) <- c('Concordant', 'Discordant')
confusion.norm <- (confusion / colSums(confusion)[col(confusion)])
confusion.norm <- as.data.frame(as.table(confusion.norm))
confusion <- as.data.frame(as.table(confusion))

pal <- colorRampPalette(c("white", "red"))
confusion.plot <- ggplot(confusion.norm, aes(as.factor(Var1), Var2, group=Var2)) + 
  geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + 
  geom_text(aes(fill =confusion.norm$Freq, label=round(confusion$Freq, 2)), size=12, fontface='bold') +
  scale_x_discrete() + 
  scale_y_discrete() + 
  scale_fill_gradientn(colours=pal(100), limits=c(0,1), breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1)) + 
  theme_bw() +
  theme(axis.text.x  = element_text(size=12), 
        axis.text.y  = element_text(size=12), 
        axis.title   = element_text(size=15),
        legend.text  = element_text(size=12), 
        legend.title = element_text(size=15, face='bold')) +
  labs(x='', y='')

CairoPNG(paste(prefix,'_confusion.png', sep=''),units='px',w=800,h=600)
confusion.plot
dev.off()
