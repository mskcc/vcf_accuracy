#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript


library(Cairo)
library(data.table)
library(ggplot2)
library(reshape)


get.percentages = function(Missed_SNP, Novel_SNP, Correct_SNP_Genotype, Incorrect_SNP_Genotype, 
                           Missed_INDEL,  Novel_INDEL, Correct_INDEL_Genotype, Incorrect_INDEL_Genotype,  
                           Total_Truth_SNPs,	Total_Truth_INDELs, Total_Test_SNPs, Total_Test_INDELs){
  
  if(Total_Truth_SNPs == 0){
    
    PC_Missed_SNP =0
    PC_Correct_SNP_Genotype =0
    PC_Incorrect_SNP_Genotype =0
  }
  
  if(Total_Truth_INDELs == 0){
    
    PC_Missed_INDEL =0
    PC_Correct_INDEL_Genotype =0
    PC_Incorrect_INDEL_Genotype =0
  }
  
  if(Total_Test_SNPs == 0){
    
    PC_Novel_SNP =0
  }
  
  if(Total_Test_INDELs == 0){
    
    PC_Novel_INDEL =0
  }
  
  if(Total_Truth_SNPs > 0){
    PC_Missed_SNP = (Missed_SNP / Total_Truth_SNPs) *100
    PC_Correct_SNP_Genotype = (Correct_SNP_Genotype / Total_Truth_SNPs) *100
    PC_Incorrect_SNP_Genotype = (Incorrect_SNP_Genotype / Total_Truth_SNPs) *100
  }
  
  if(Total_Truth_INDELs > 0){
    PC_Missed_INDEL = (Missed_INDEL / Total_Truth_INDELs) *100
    PC_Correct_INDEL_Genotype = (Correct_INDEL_Genotype / Total_Truth_INDELs) *100
    PC_Incorrect_INDEL_Genotype =  (Incorrect_INDEL_Genotype / Total_Truth_INDELs) *100
  }
  
  ###################################################################################################################
  ###################################################################################################################
  
  if(Total_Test_SNPs > 0){
    PC_Novel_SNP = (Novel_SNP / Total_Test_SNPs) *100
  }
  
  if(Total_Test_INDELs > 0){
    PC_Novel_INDEL = (Novel_INDEL / Total_Test_INDELs) *100
  }
  
  ###################################################################################################################
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
filename = args[1]
prefix = args[2]
samples = fread(filename)

samples[,c('Missed SNP', 
           'Novel SNP', 
           'Correct SNP Genotype', 
           'Incorrect SNP Genotype',
           'Missed INDEL', 
           'Novel INDEL',
           'Correct INDEL Genotype' ,
           'Incorrect INDEL Genotype'):=get.percentages(Missed_SNP, Novel_SNP, Correct_SNP_Genotype, Incorrect_SNP_Genotype, 
                                                        Missed_INDEL,  Novel_INDEL, Correct_INDEL_Genotype, Incorrect_INDEL_Genotype,	
                                                        Total_Truth_SNPs,	Total_Truth_INDELs, Total_Test_SNPs, Total_Test_INDELs), by=1:nrow(samples)]

sample = melt(samples, id.vars=c('sample'), measure.vars=c('Missed SNP', 'Novel SNP', 'Correct SNP Genotype', 'Incorrect SNP Genotype',  
                                                           'Missed INDEL', 'Correct INDEL Genotype', 'Incorrect INDEL Genotype'))


my.plot = ggplot(sample, aes(y=value, x =variable, fill =sample)) +  
  geom_bar(stat = 'identity', position = 'dodge')  +
  theme(axis.text.x  = element_text(angle=90, vjust=.6, size=8), legend.position='right')  +
  labs(x='', y='Percentage of variants')

CairoPNG(paste(prefix,'_percentages.png'), sep='',units='px',w=800,h=600)
my.plot
dev.off()

