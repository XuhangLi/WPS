library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
#

# ############# quick identification by PCA with marker genes###########################################################
# ## load and merge all data
# load("./outputs/TPM_preserved_met5_lib6.RData")
# 
# library(edgeR)
# library(dplyr)
# library(tibble)
# input <- readsCount
# coldata =  data.frame(sampleName = colnames(readsCount))
# rownames(coldata) = colnames(input)
# coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
# library(DESeq2)
# dds <- DESeqDataSetFromMatrix(countData = input,
#                               colData = coldata,
#                               design= ~ sampleName)
# 
# # filtering 
# keep <- rowSums(counts(dds)>=10) >= 1
# dds <- dds[keep,]
# 
# #exploratory data analysis (EDA) plot, ie, PCA
# vsd <- vst(dds, blind=TRUE)
# # only pca on feature gene
# set0 = read.csv('./../../../DevelopmentalStage/wrongsampleFeatureGene.csv')
# set0 = set0$x
# selectedData = as.data.frame(assay(vsd))
# selectedData = selectedData[set0,]
# # selectedData = selectedData[,-which(colnames(selectedData) %in% c('x.dib_1_rep1','x.dib_1_rep2','x.dib_1_rep3'))]
# pca <- prcomp(t(selectedData))
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
# percentVar <- round(100 * percentVar)
# pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(selectedData))
# ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
#   geom_point(color = "red") +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   geom_text_repel()
# 
# # confirm suspicious ones
# dataMat = as.data.frame(assay(vsd))
# plot(dataMat$x.adk_1_rep2,dataMat$x.adk_1_rep3)
# #x.vector_well2_rep1
# 
# badSamples = list()
# badSamples[['lib1']] = c('x.ard_1_rep1')
# 
# badSamples[['lib2']] = c('x.pycr_1_rep2')
# 
# badSamples[['lib3']] = c('x.gln_6_rep2')
# 
# badSamples[['lib4']] = c('NOT ANALYZED')
# 
# badSamples[['lib5']] = c('x.afmd_1_rep1','x.F01D4.8_rep2','x.gldc_1_rep2')
# 
# badSamples[['lib6']] = c('x.K01C8.1_rep2')
# 
# save(list = 'badSamples',file = 'outputs/badSample_met3.Rdata')

############# RAW EDA ###########################################################
# just exploratory data analysis to look at data without any treatment
library(edgeR)
library(dplyr)
library(tibble)
libs = c('lib1','lib2','lib3','lib4','lib5','lib6')
exps = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10', 'met11','met13')
skipLib = c() 
for (expID in exps){
  for (libID in libs){
    if(!(paste(expID,libID,sep = '-') %in% skipLib)){
      ## load and merge all data
      load(paste("./outputs/TPM_preserved_",expID,'_',libID,".RData",sep = ''))
      
      # remove the low depth in order not to bias normalization
      readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
      
      input <- readsCount
      coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
      rownames(coldata) = colnames(input)
      coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
      coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
      coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
      
      library(DESeq2)
      dds <- DESeqDataSetFromMatrix(countData = input,
                                    colData = coldata,
                                    design= ~ rep + RNAi)
      
      # filtering 
      keep <- rowSums(counts(dds)>=10) >= 1
      dds <- dds[keep,]
      
      vsd <- vst(dds, blind=TRUE)
      mat <- assay(vsd)
      # mat <- limma::removeBatchEffect(mat, batch = vsd$rep)
      if ('x.vector' %in% unique(as.character(vsd$RNAi))){
        allRNAi = c('x.vector',setdiff(unique(as.character(vsd$RNAi)),'x.vector'))
      }else{
        allRNAi = unique(as.character(vsd$RNAi))
      }
      dataMat = as.data.frame(mat)
      pdf(paste("figures/",expID,'_',libID,'_EDA.pdf',sep = '')) 
      for (i in 1:length(allRNAi)){
        if ('x.vector' %in% allRNAi){
          mat_sub = mat[,dds$RNAi %in% c('x.vector',allRNAi[i])]
        }else{
          mat_sub = mat[,dds$RNAi %in% c(allRNAi[i])]
        }
        
        pca <- prcomp(t(mat_sub))
        percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
        percentVar <- round(100 * percentVar)
        pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat_sub))
        if ('x.vector' %in% allRNAi){
          pcaData$RNAi = dds$RNAi[dds$RNAi %in% c('x.vector',allRNAi[i])]
          pcaData$batch = dds$batchLabel[dds$RNAi %in% c('x.vector',allRNAi[i])]
        }else{
          pcaData$RNAi = dds$RNAi[dds$RNAi %in% c(allRNAi[i])]
          pcaData$batch = dds$batchLabel[dds$RNAi %in% c(allRNAi[i])]
        }
        p = ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
          geom_point(aes(color = RNAi,), size = 3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          scale_color_manual(values=c("grey", "red")) +
          geom_text_repel()
        print(p)
      }
      while (!is.null(dev.list())) {dev.off()}
    }
  }
}

############# gold-standard identification by RNAi-control paired PCA ###########################################################
# this methold may also identify the accidental RNAi failure in one replicate
# first, generate the PCA plots
libs = c('lib1','lib2','lib3','lib4','lib5','lib6')
exps = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11')
skipLib = c() # no vector, we dont analyze 
for (expID in exps){
  for (libID in libs){
    if(!(paste(expID,libID,sep = '-') %in% skipLib)){
      ## load and merge all data
      load(paste("./outputs/TPM_preserved_",expID,'_',libID,".RData",sep = ''))
      
      # remove the low depth in order not to bias normalization
      readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
      
      input <- readsCount
      coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
      rownames(coldata) = colnames(input)
      coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
      coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
      coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
      
      library(DESeq2)
      dds <- DESeqDataSetFromMatrix(countData = input,
                                    colData = coldata,
                                    design= ~ rep + RNAi)
      
      # filtering 
      keep <- rowSums(counts(dds)>=10) >= 1
      dds <- dds[keep,]
      
      vsd <- vst(dds, blind=FALSE)
      mat <- assay(vsd)
      mat <- limma::removeBatchEffect(mat, batch = vsd$rep, design= model.matrix(~ vsd$RNAi))
      
      if ('x.vector' %in% unique(as.character(vsd$RNAi))){
        allRNAi = c('x.vector',setdiff(unique(as.character(vsd$RNAi)),'x.vector'))
      }else{
        allRNAi = unique(as.character(vsd$RNAi))
      }

      dataMat = as.data.frame(mat)
      pdf(paste("figures/",expID,'_',libID,'_badSamplePCA.pdf',sep = '')) 
      for (i in 1:length(allRNAi)){
        if ('x.vector' %in% allRNAi){
          mat_sub = mat[,dds$RNAi %in% c('x.vector',allRNAi[i])]
        }else{
          mat_sub = mat[,dds$RNAi %in% c(allRNAi[i])]
        }
        pca <- prcomp(t(mat_sub))
        percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
        percentVar <- round(100 * percentVar)
        pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat_sub))
        if ('x.vector' %in% allRNAi){
          pcaData$RNAi = dds$RNAi[dds$RNAi %in% c('x.vector',allRNAi[i])]
          pcaData$batch = dds$batchLabel[dds$RNAi %in% c('x.vector',allRNAi[i])]
        }else{
          pcaData$RNAi = dds$RNAi[dds$RNAi %in% c(allRNAi[i])]
          pcaData$batch = dds$batchLabel[dds$RNAi %in% c(allRNAi[i])]
        }
        p = ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
          geom_point(aes(color = RNAi,), size = 3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          scale_color_manual(values=c("grey", "red")) +
          geom_text_repel()+
          ggtitle(paste(paste(expID,libID,i,allRNAi[i],sep = '-'),'\n',
                        paste(names(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]])),collapse = ' '),'\n',
                        paste(as.character(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]]) /1e6),collapse = ' '),sep = '')
                  )
        print(p)
        
        sampleDists <- dist(t(mat_sub))
        sampleDistMatrix <- as.matrix(sampleDists)
        library("RColorBrewer")
        colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
        pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,col=colors)
        
      }
      while (!is.null(dev.list())) {dev.off()}
      pdf(paste("figures/",expID,'_',libID,'_badSampleCorr.pdf',sep = '')) 
      for (z in 1:length(allRNAi)){
        mat_sub = as.data.frame(mat[,dds$RNAi %in% allRNAi[z]])
        vectorsamples = colnames(mat_sub)
        for (j in 1:(length(vectorsamples)-1)){
          for (i in (j+1):length(vectorsamples)){
              rsqr = round(cor(as.matrix(mat_sub[vectorsamples[j]]),as.matrix(mat_sub[vectorsamples[i]]))^2,4)
              plot(as.matrix(mat_sub[vectorsamples[j]]),as.matrix(mat_sub[vectorsamples[i]]),main = paste(allRNAi[z],rsqr),xlab = vectorsamples[j], ylab = vectorsamples[i])
            }
          }
        }
      while (!is.null(dev.list())) {dev.off()}
    }
  }
}

############# inspection and label the quality level library-by-library #########
# our criteria is we allow the natural noisiness and only look for failed sample that is due to technical reason; 
# such sample usually show high distance to any other samples or show strong difference with its other two replicates
# that cannot be explained by biological variations (i.e., other two are close to vectors but this one is very far; or 
# vice versa); we also consider a sample as bad when it show strong de-correlation (i.e. r^2 < 0.95)
# sometimes we see PCA/distance outliers because a sample's depth is low; this usually is caused by lowly expressed genes and 
# does not mean the sample is bad. we still include them to maximize the power (low-depth replicate is still beneficial when # of
# replicate is low); but we exlude the low-depth sample if its variance cannot be explained by lack of depth only 
qLevels <- function(expID, libID, badSet){ #checkSet
  load(paste("./outputs/TPM_preserved_",expID,'_',libID,".RData",sep = ''))
  
  # remove the low depth in order not to bias normalization
  readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
  
  library(edgeR)
  library(dplyr)
  library(tibble)
  input <- readsCount
  coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
  rownames(coldata) = colnames(input)
  coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
  coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
  coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = input,
                                colData = coldata,
                                design= ~ rep + RNAi)
  
  # filtering 
  keep <- rowSums(counts(dds)>=10) >= 1
  dds <- dds[keep,]
  
  # vsd <- vst(dds, blind=TRUE)
  # mat <- assay(vsd)
  # mat <- limma::removeBatchEffect(mat, batch = vsd$rep)
  # allRNAi = setdiff(unique(as.character(vsd$RNAi)),'x.vector')
  # dataMat = as.data.frame(mat)
  dds$Qlevel = 1
  # pdf(paste("figures/",expID,'_',libID,'_badSampleTPM.pdf',sep = '')) 
  # for (i in 1:length(checkSet)){
  #   bs = checkSet[i]
  #   RNAi = str_replace(bs,'_rep.$','_rep')
  #   vectorsamples = colnames(dataMat)[str_detect(colnames(dataMat),'x.vector')]
  #   othersamples = setdiff(colnames(dataMat)[str_detect(colnames(dataMat),RNAi)],bs)
  #   for (j in 1:length(othersamples)){
  #     plot(as.matrix(dataMat[bs]),as.matrix(dataMat[othersamples[j]]),main = paste('self - ',bs,sep = ''),xlab = bs, ylab = othersamples[j])
  #   }
  #   plot(as.matrix(dataMat[othersamples[1]]),as.matrix(dataMat[othersamples[2]]),main = paste('other two - ',othersamples[1],sep = ''),xlab = bs, ylab = othersamples[2])
  #   for (j in 1:length(vectorsamples)){
  #     plot(as.matrix(dataMat[bs]),as.matrix(dataMat[vectorsamples[j]]),main = paste('vector - ',bs,sep = ''),xlab = bs, ylab = vectorsamples[j])
  #   }
  # }
  # vectorsamples = colnames(dataMat)[str_detect(colnames(dataMat),'x.vector')]
  # for (j in 1:(length(vectorsamples)-1)){
  #   for (i in (j+1):length(vectorsamples))
  #   plot(as.matrix(dataMat[vectorsamples[j]]),as.matrix(dataMat[vectorsamples[i]]),main = 'vector vs vector',xlab = vectorsamples[j], ylab = vectorsamples[i])
  # }
  dds$Qlevel[colnames(dds) %in% badSet] = 4
  if (!all(badSet %in% colnames(dds))){
    stop('sample ID error!!!')
  }
  # 
  # # for replot
  # if (length(badSet)>0){
  #   for (i in 1:length(badSet)){
  #     bs = badSet[i]
  #     RNAi = str_replace(bs,'_rep.$','')
  #     mat_sub = mat[,dds$RNAi %in% c('x.vector',RNAi)]
  #     rmID = which(colnames(mat_sub) %in% badSet)
  #     mat_sub = mat_sub[,-rmID]
  #     pca <- prcomp(t(mat_sub))
  #     percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  #     percentVar <- round(100 * percentVar)
  #     pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat_sub))
  #     pcaData$RNAi = dds$RNAi[dds$RNAi %in% c('x.vector',RNAi)][-rmID]
  #     pcaData$batch = dds$batchLabel[dds$RNAi %in% c('x.vector',RNAi)][-rmID]
  #     p <- ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
  #       geom_point(aes(color = RNAi), size = 3) +
  #       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #       ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #       scale_color_manual(values=c("grey", "red")) +
  #       geom_text_repel() + 
  #       ggtitle(paste(paste(expID,libID,i,RNAi,sep = '-')))
  #     print(p)
  #     
  #     sampleDists <- dist(t(mat_sub))
  #     sampleDistMatrix <- as.matrix(sampleDists)
  #     library("RColorBrewer")
  #     colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  #     pheatmap(sampleDistMatrix,
  #              clustering_distance_rows=sampleDists,
  #              clustering_distance_cols=sampleDists,col=colors)
  #   }
  # }
  # save the quality levels
  qLevels =  data.frame(sampleName = colnames(dds), QL = dds$Qlevel)
  rownames(qLevels) = qLevels$sampleName
  qLevels$QL[qLevels$QL == 1] = 'good'
  qLevels$QL[qLevels$QL == 4] = 'bad'
  write.csv(qLevels,paste('outputs/qualityLevels_',expID,'_', libID,'.csv',sep = ''))
  # while (!is.null(dev.list())) {dev.off()}
}
finalPCA <- function(expID, libID, skipLib){
  if(!(paste(expID,libID,sep = '-') %in% skipLib)){
    ## load and merge all data
    load(paste("./outputs/TPM_preserved_",expID,'_',libID,".RData",sep = ''))
    
    # remove the low depth in order not to bias normalization
    readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
    # remove the bad ones
    qLevels = read.csv(paste('outputs/qualityLevels_',expID,'_', libID,'.csv',sep = ''))
    readsCount = readsCount[, colnames(readsCount) %in% qLevels$sampleName[qLevels$QL == 'good']]
    
    input <- readsCount
    coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
    rownames(coldata) = colnames(input)
    coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
    coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
    coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
    
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = input,
                                  colData = coldata,
                                  design= ~ rep + RNAi)
    
    # filtering 
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]
    
    vsd <- vst(dds, blind=TRUE)
    mat <- assay(vsd)
    mat <- limma::removeBatchEffect(mat, batch = vsd$rep)
    
    if ('x.vector' %in% unique(as.character(vsd$RNAi))){
      allRNAi = c('x.vector',setdiff(unique(as.character(vsd$RNAi)),'x.vector'))
    }else{
      allRNAi = unique(as.character(vsd$RNAi))
    }
    dataMat = as.data.frame(mat)
    pdf(paste("figures/",expID,'_',libID,'_cleanSamplePCA.pdf',sep = '')) 
    for (i in 1:length(allRNAi)){
      if ('x.vector' %in% allRNAi){
        mat_sub = mat[,dds$RNAi %in% c('x.vector',allRNAi[i])]
      }else{
        mat_sub = mat[,dds$RNAi %in% c(allRNAi[i])]
      }
      pca <- prcomp(t(mat_sub))
      percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
      percentVar <- round(100 * percentVar)
      pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat_sub))
      if ('x.vector' %in% allRNAi){
        pcaData$RNAi = dds$RNAi[dds$RNAi %in% c('x.vector',allRNAi[i])]
        pcaData$batch = dds$batchLabel[dds$RNAi %in% c('x.vector',allRNAi[i])]
      }else{
        pcaData$RNAi = dds$RNAi[dds$RNAi %in% c(allRNAi[i])]
        pcaData$batch = dds$batchLabel[dds$RNAi %in% c(allRNAi[i])]
      }
      p = ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
        geom_point(aes(color = RNAi,), size = 3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        scale_color_manual(values=c("grey", "red")) +
        geom_text_repel()+
        ggtitle(paste(paste(expID,libID,i,allRNAi[i],sep = '-'),'\n',
                      paste(names(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]])),collapse = ' '),'\n',
                      paste(as.character(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]]) /1e6),collapse = ' '),sep = '')
        )
      print(p)
      
      sampleDists <- dist(t(mat_sub))
      sampleDistMatrix <- as.matrix(sampleDists)
      library("RColorBrewer")
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,col=colors)
    }
    while (!is.null(dev.list())) {dev.off()}
  }
}

skipLib = c() # no vector, we dont analyze 

expID = 'met13' 
libID = 'lib6'
badSet = c()#, 
#checkSet =c('x.gpdh_1_rep1','x.ERROR_let_767_rep1')
qLevels(expID, libID, badSet)
# final PCA plots with bad removed
finalPCA(expID, libID, skipLib)

######################

