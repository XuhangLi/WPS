

###
cutoff_all = c(1, 0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005, 1e-04, 1e-5, 1e-10, 1e-30, 0) # better revise this in the final run 

strvector = c('1','0.5','0.2','0.1','0.075','0.05', '0.025','0.01', '0.0075','0.005','0.0025', '0.001','0.00075','5e-04','1e-04', '1e-05', '1e-10','1e-30', '0')
FCtypes = c("log2FoldChange","log2FoldChange_raw")

for (FCtype in FCtypes){
#str = '0.005'
  for (str in strvector){
  library(DESeq2)
  ############# make the pvalue and fc matrix ######
  plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
  allConditions = c()
  allgenes = c()
  for (plateID in plateIDs){
    # load the FC matrix from DE
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
    allgenes = union(allgenes, rownames(input))
  }
  RNAiName = c()
  batchID = c()
  for (i in 1:length(allConditions)){
    RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
    batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])
    
  }
  # remove vectors
  batchID = batchID[RNAiName != 'x.vector']
  RNAiName = RNAiName[RNAiName != 'x.vector']
  

  
  #################### re-cutoff optimal threshold: using the vectorlike samples to find best FC and FDR cutoff for DE calling  ######################
  library(stringr)
  plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
  allConditions = c()
  for (plateID in plateIDs){
    # load the FC matrix from DE
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
  }
  RNAiName = c()
  batchID = c()
  for (i in 1:length(allConditions)){
    RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
    batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])
    
  }
  

  # load a set of manually curated vectorlike conditions
  vectorlikeTbl = read.csv('./../../MetabolicLibrary/input_data/metaData/manually_curated_vectorlike_conditions.csv')
  vectorlikeTbl$vectorlike_RNAi = str_remove(vectorlikeTbl$NTP_condID,' met[0-9]+_lib.$')
  vectorlikes = RNAiName[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                         | str_detect(RNAiName, '^x.REALVECTOR_')]
  batch_vectorlikes = batchID[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                              | str_detect(RNAiName, '^x.REALVECTOR_')]
  # some of the realvectors are technical replicates, we remove them
  rmInd = str_detect(vectorlikes, '^x.REALVECTOR_') & batch_vectorlikes %in% c('met11_lib1','met11_lib2','met11_lib3','')
  vectorlikes = vectorlikes[!rmInd]
  batch_vectorlikes = batch_vectorlikes[!rmInd]
  
  
  # aggregate the raw tables for cutoff titrating
  mergeTbl = data.frame()
  for (i in 1:length(vectorlikes)){
    tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',vectorlikes[i],'_',batch_vectorlikes[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))  
    tbl$RNAi = vectorlikes[i]
    tbl$batchID = batch_vectorlikes[i]
    mergeTbl = rbind(mergeTbl, tbl)
  }
  maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
  # N_DE_list = list()
  fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
  fc_seq = seq(1,3,0.1)
  
  N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)), 
                     c(length(fdr_seq), length(fc_seq)))
  for (FDR in fdr_seq){
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < FDR),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = quantile(N_DE,0.9) #
    }
  }
  
  library(pheatmap)
  #dev.off()
  pdf(paste('figures/0_DE_QA/2d_cutoff_titration_71NTP_p',str,'_',FCtype,'.pdf',sep = ''),height = 10,width = 10) 
  pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
           breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]
  
  # evaluate FC cutoff at a fixed FDR
  # we end up with fdr and eFDR 0.05
  N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
  colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
  for (FC in fc_seq){
    passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < 0.1),]
    passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
    N_DE = table(passTbl$ID)
    zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
    tmp = rep(0, length(zeros))
    names(tmp) = zeros
    N_DE = c(N_DE, tmp)
    N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
  }
  plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
  for (i in 2:ncol(N_DE_fdr01)){
    points(fc_seq, N_DE_fdr01[,i],type = 'o')
  }
  abline(h=5, col = 'red')
  dev.off()
  
  # do this evaluation with real vectors
  realvectors = RNAiName[str_detect(RNAiName, '^x.REALVECTOR_')]
  batch_realvectors = batchID[str_detect(RNAiName, '^x.REALVECTOR_')]
  rmInd = str_detect(realvectors, '^x.REALVECTOR_') & batch_realvectors %in% c('met11_lib1','met11_lib2','met11_lib3','')
  realvectors = realvectors[!rmInd]
  batch_realvectors = batch_realvectors[!rmInd]
  vectorlikes = realvectors
  batch_vectorlikes = batch_realvectors
  
  # aggregate the raw tables for cutoff titrating
  mergeTbl = data.frame()
  for (i in 1:length(vectorlikes)){
    tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',vectorlikes[i],'_',batch_vectorlikes[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))  
    tbl$RNAi = vectorlikes[i]
    tbl$batchID = batch_vectorlikes[i]
    mergeTbl = rbind(mergeTbl, tbl)
  }
  maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
  # N_DE_list = list()
  fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
  fc_seq = seq(1,3,0.1)
  
  N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)), 
                     c(length(fdr_seq), length(fc_seq)))
  for (FDR in fdr_seq){
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < FDR),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = mean(N_DE) #
    }
  }
  
  library(pheatmap)
  #dev.off()
  pdf(paste('figures/0_DE_QA/2d_cutoff_titration_realVector_techRep_removed_p',str,'_',FCtype,'.pdf',sep = ''),height = 10,width = 10) 
  pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
           breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]
  
  # evaluate FC cutoff at a fixed FDR
  # we end up with fdr and eFDR 0.05
  N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
  colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
  for (FC in fc_seq){
    passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < 0.1),]
    passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
    N_DE = table(passTbl$ID)
    zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
    tmp = rep(0, length(zeros))
    names(tmp) = zeros
    N_DE = c(N_DE, tmp)
    N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
  }
  plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
  for (i in 2:ncol(N_DE_fdr01)){
    points(fc_seq, N_DE_fdr01[,i],type = 'o')
  }
  abline(h=5, col = 'red')
  dev.off()
}
}

# plot that for the uncleaned data 
rm(list = ls())
plotControl = TRUE
FCtypes = c("log2FoldChange","log2FoldChange_raw")
for (FCtype in FCtypes){
  if (plotControl){ # vectorlike and realvector titration
  library(stringr)
  plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
  allConditions = c()
  for (plateID in plateIDs){
    # load the FC matrix from DE
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
  }
  RNAiName = c()
  batchID = c()
  for (i in 1:length(allConditions)){
    RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
    batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])
    
  }
  
  # load a set of manually curated vectorlike conditions
  vectorlikeTbl = read.csv('./../../MetabolicLibrary/input_data/metaData/manually_curated_vectorlike_conditions.csv')
  vectorlikeTbl$vectorlike_RNAi = str_remove(vectorlikeTbl$NTP_condID,' met[0-9]+_lib.$')
  vectorlikes = RNAiName[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                         | str_detect(RNAiName, '^x.REALVECTOR_')]
  batch_vectorlikes = batchID[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                              | str_detect(RNAiName, '^x.REALVECTOR_')]
  # some of the realvectors are technical replicates, we remove them
  rmInd = str_detect(vectorlikes, '^x.REALVECTOR_') & batch_vectorlikes %in% c('met11_lib1','met11_lib2','met11_lib3','')
  vectorlikes = vectorlikes[!rmInd]
  batch_vectorlikes = batch_vectorlikes[!rmInd]
  
  # remove met3_lib4
  vectorlikes = vectorlikes[!(batch_vectorlikes %in% 'met3_lib4')]
  batch_vectorlikes = batch_vectorlikes[!(batch_vectorlikes %in% 'met3_lib4')]
  
  # aggregate the raw tables for cutoff titrating
  mergeTbl = data.frame()
  for (i in 1:length(vectorlikes)){
    tbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_alpha005_RNAi_',vectorlikes[i],'_vs_x.vector_',batch_vectorlikes[i],'.csv',sep = ''))  
    tbl$RNAi = vectorlikes[i]
    tbl$batchID = batch_vectorlikes[i]
    mergeTbl = rbind(mergeTbl, tbl)
  }
  maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
  # N_DE_list = list()
  fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
  fc_seq = seq(1,3,0.1)
  
  N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)), 
                     c(length(fdr_seq), length(fc_seq)))
  for (FDR in fdr_seq){
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$padj < FDR),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = quantile(N_DE,0.9) #
    }
  }
  
  library(pheatmap)
  #dev.off()
  pdf(paste('figures/0_DE_QA/2d_cutoff_titration_71NTP_rawDE_log2FoldChange_',FCtype,'.pdf',sep = ''),height = 10,width = 10) 
  pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
           breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]
  
  # evaluate FC cutoff at a fixed FDR
  # we end up with fdr and eFDR 0.05
  N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
  colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
  for (FC in fc_seq){
    passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$padj < 0.1),]
    passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
    N_DE = table(passTbl$ID)
    zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
    tmp = rep(0, length(zeros))
    names(tmp) = zeros
    N_DE = c(N_DE, tmp)
    N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
  }
  plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
  for (i in 2:ncol(N_DE_fdr01)){
    points(fc_seq, N_DE_fdr01[,i],type = 'o')
  }
  abline(h=5, col = 'red')
  dev.off()
  
  # do this evaluation with real vectors
  realvectors = RNAiName[str_detect(RNAiName, '^x.REALVECTOR_')]
  batch_realvectors = batchID[str_detect(RNAiName, '^x.REALVECTOR_')]
  rmInd = str_detect(realvectors, '^x.REALVECTOR_') & batch_realvectors %in% c('met11_lib1','met11_lib2','met11_lib3','')
  realvectors = realvectors[!rmInd]
  batch_realvectors = batch_realvectors[!rmInd]
  vectorlikes = realvectors
  batch_vectorlikes = batch_realvectors
  
  # aggregate the raw tables for cutoff titrating
  mergeTbl = data.frame()
  for (i in 1:length(vectorlikes)){
    tbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_alpha005_RNAi_',vectorlikes[i],'_vs_x.vector_',batch_vectorlikes[i],'.csv',sep = ''))  
    tbl$RNAi = vectorlikes[i]
    tbl$batchID = batch_vectorlikes[i]
    mergeTbl = rbind(mergeTbl, tbl)
  }
  maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
  # N_DE_list = list()
  fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
  fc_seq = seq(1,3,0.1)
  
  N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)), 
                     c(length(fdr_seq), length(fc_seq)))
  for (FDR in fdr_seq){
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$padj < FDR),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = mean(N_DE) #
    }
  }
  
  library(pheatmap)
  #dev.off()
  pdf(paste('figures/0_DE_QA/2d_cutoff_titration_realVector_techRep_removed_rawDE_log2FoldChange_',FCtype,'.pdf',sep = ''),height = 10,width = 10) 
  pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
           breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]
  
  # evaluate FC cutoff at a fixed FDR
  # we end up with fdr and eFDR 0.05
  N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
  colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
  for (FC in fc_seq){
    passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$padj < 0.1),]
    passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
    N_DE = table(passTbl$ID)
    zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
    tmp = rep(0, length(zeros))
    names(tmp) = zeros
    N_DE = c(N_DE, tmp)
    N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
  }
  plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
  for (i in 2:ncol(N_DE_fdr01)){
    points(fc_seq, N_DE_fdr01[,i],type = 'o')
  }
  abline(h=5, col = 'red')
  dev.off()
  
}
}

################ distribution of #DEG calls in vectorlike without cleaning ###################
rm(list = ls())

# define DEG cutoff
FDRcutoff = 0.01
FCcutoff = 2

# load the vector like information 
library(stringr)
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
allConditions = c()
for (plateID in plateIDs){
  # load the FC matrix from DE
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
}
RNAiName = c()
batchID = c()
for (i in 1:length(allConditions)){
  RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
  batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])
  
}

# load a set of manually curated vectorlike conditions
vectorlikeTbl = read.csv('./../../MetabolicLibrary/input_data/metaData/manually_curated_vectorlike_conditions.csv')
vectorlikeTbl$vectorlike_RNAi = str_remove(vectorlikeTbl$NTP_condID,' met[0-9]+_lib.$')
vectorlikes = RNAiName[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                       | str_detect(RNAiName, '^x.REALVECTOR_')]
batch_vectorlikes = batchID[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                            | str_detect(RNAiName, '^x.REALVECTOR_')]
# some of the realvectors are technical replicates, we remove them
rmInd = str_detect(vectorlikes, '^x.REALVECTOR_') & batch_vectorlikes %in% c('met11_lib1','met11_lib2','met11_lib3','')
vectorlikes = vectorlikes[!rmInd]
batch_vectorlikes = batch_vectorlikes[!rmInd]

# do this evaluation with real vectors
realvectors = RNAiName[str_detect(RNAiName, '^x.REALVECTOR_')]
batch_realvectors = batchID[str_detect(RNAiName, '^x.REALVECTOR_')]
rmInd = str_detect(realvectors, '^x.REALVECTOR_') & batch_realvectors %in% c('met11_lib1','met11_lib2','met11_lib3','')
realvectors = realvectors[!rmInd]
batch_realvectors = batch_realvectors[!rmInd]

# calculate the DE distribution by all spikein controls
N_DE_vectorlike = c()
for (i in 1:length(vectorlikes)){
  tbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_alpha005_RNAi_',vectorlikes[i],'_vs_x.vector_',batch_vectorlikes[i],'.csv',sep = ''))  
  N_DE_vectorlike = c(N_DE_vectorlike, sum(abs(tbl$log2FoldChange_raw) > log2(FCcutoff) & tbl$padj < FDRcutoff,na.rm = T))
} 

#hist(N_DE_vectorlike)
# plot
library(ggplot2)
# Customize the theme for publication quality
my_theme <- theme_bw() +
  theme(
    text = element_text(face = "plain"),  # Remove bold text from all elements
    axis.title = element_text(size = 14),  # Set axis title text size
    axis.text = element_text(size = 12),  # Set axis tick label text size
    legend.title = element_text(size = 12),  # Set legend title text size
    legend.text = element_text(size = 10),  # Set legend item text size
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank(),  # Remove box around the plot
    axis.line = element_line(color = "black")  # Change axis color to black
  )

# Create a histogram using ggplot2 with custom theme
histogram1 <- ggplot(data.frame(N_DE_vectorlike), aes(x = N_DE_vectorlike)) +
  geom_histogram( color = "black", fill = "grey") +
  labs(
    title = 'realvector + NTP (71)',
    x = "number of DEG calls",
    y = "frequency",
  ) +   geom_vline(aes(xintercept = 5), linetype = "dashed", color = "red") + 
  ylim(c(0,13))+
  my_theme

# calculate the DE distribution by real vector spikein controls
N_DE_realVector = c()
for (i in 1:length(realvectors)){
  tbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_alpha005_RNAi_',realvectors[i],'_vs_x.vector_',batch_realvectors[i],'.csv',sep = ''))  
  N_DE_realVector = c(N_DE_realVector, sum(abs(tbl$log2FoldChange_raw) > log2(FCcutoff) & tbl$padj < FDRcutoff,na.rm = T))
} 

# Create a histogram using ggplot2 with custom theme
histogram2 <- ggplot(data.frame(N_DE_realVector), aes(x = N_DE_realVector)) +
  geom_histogram( color = "black", fill = "grey") +
  labs(
    title = 'realvector (4)',
    x = "number of DEG calls",
    y = "frequency",
  ) +   geom_vline(aes(xintercept = 5), linetype = "dashed", color = "red") + 
  my_theme

# overlay together 
df = data.frame(N_DE_vectorlike = c(N_DE_vectorlike, N_DE_realVector),
                type = c(rep('all',length(N_DE_vectorlike)), rep('Spike-in control',length(N_DE_realVector))))
# Create a histogram using ggplot2 with custom theme
histogram3 <- ggplot(df, aes(x = N_DE_vectorlike)) +
  geom_histogram(data = df[df$type == 'all', ], aes(x = N_DE_vectorlike), color = "black", fill = 'grey') +
  geom_histogram(data = df[df$type == 'Spike-in control', ], aes(x = N_DE_vectorlike), color = "black", fill = 'red') +
  labs(
    title = 'realvector + NTP (71)',
    x = "number of DEG calls",
    y = "frequency",
  ) +   geom_vline(aes(xintercept = 5), linetype = "dashed", color = "red") + 
  ylim(c(0,13))+
  my_theme

# Print the histogram
dev.off()
pdf('figures/0_DE_QA/vectorlike_analysis/N_DE_distribution_71NTP.pdf', width = 3, height = 3)
print(histogram1)
print(histogram2)
print(histogram3)
dev.off()
