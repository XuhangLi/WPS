# ExtraVectorPlates = c('met2','met3','met4','met5')

# # prepare all vectors 
# load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',ExtraVectorPlates[1],'.Rdata',sep = ''))
# batches = levels(batchLabel)
# #pick the deepest control 
# a = sort(colSums(input)[grep('^x.vector',colnames(input))])/1e6
# aveDph = c()
# for (i in 1:length(batches)){
#   aveDph = c(aveDph, mean(a[grep(paste(batches[i],'$',sep = ''),names(a))]))
# }
# maxBatch = which(aveDph == max(aveDph))
# input_extr_ctr = input[,RNAi %in% 'x.vector' & batchLabel %in% batches[maxBatch]]
# batchLabel_extr_ctr = as.character(batchLabel[RNAi %in% 'x.vector' & batchLabel %in% batches[maxBatch]])
# for (i in 2:length(ExtraVectorPlates)){
#   load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',ExtraVectorPlates[i],'.Rdata',sep = ''))
#   batches = levels(batchLabel)
#   #pick the deepest control 
#   a = sort(colSums(input)[grep('^x.vector',colnames(input))])/1e6
#   aveDph = c()
#   for (i in 1:length(batches)){
#     aveDph = c(aveDph, mean(a[grep(paste(batches[i],'$',sep = ''),names(a))]))
#   }
#   maxBatch = which(aveDph == max(aveDph))
#   input_extr_ctr = input_extr_ctr %>% rownames_to_column('gene') %>%
#     inner_join(input[,RNAi %in% 'x.vector' & batchLabel %in% batches[maxBatch]] %>% rownames_to_column('gene'), by = 'gene')
#   rownames(input_extr_ctr) = input_extr_ctr$gene
#   input_extr_ctr = input_extr_ctr[,-1]
#   batchLabel_extr_ctr = c(batchLabel_extr_ctr, as.character(batchLabel[RNAi %in% 'x.vector' & batchLabel %in% batches[maxBatch]]))
# }

library(dplyr)
library(tibble)
library(stringr)

# plateID = 'met1'
# we skip the met3_lib4
for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
  # load the target library set
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  batches = levels(batchLabel)
  library(DESeq2)
  library(limma)
  library(edgeR)
  #mergeTbl = data.frame()
  for (i in 1:length(batches)){ 
    if (!(batches[i] %in% c('met3_lib4','met2_lib5'))){
      # pick the control samples (in the batches belonging to target DE group)
      subsetInd = batchLabel == batches[i]
      # merge selected library with the extra controls
      input_subset = input[,subsetInd] # %>% rownames_to_column('gene') %>%
      # inner_join(input_extr_ctr %>% rownames_to_column('gene'), by = 'gene')
      # rownames(input_subset) = input_subset$gene
      # input_subset = input_subset[,-1]
      batchLabel_subset = colnames(input_subset)
      batchLabel_subset = str_extract(batchLabel_subset,'_rep._met[0-9]+_')
      batchLabel_subset = str_replace(batchLabel_subset,'_met._$','')
      batchLabel_subset = as.factor(batchLabel_subset)
      # batchLabel_subset = as.factor(c(as.character(batchLabel[subsetInd]), batchLabel_extr_ctr))
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      # RNAi_subset[(1+length(RNAi_subset)):(length(batchLabel_extr_ctr)+length(RNAi_subset))] = 'x.vector'
      
      coldata =  data.frame(batchLabel = batchLabel_subset, RNAi = RNAi_subset)
      rownames(coldata) = colnames(input_subset)
      dds <- DESeqDataSetFromMatrix(countData = input_subset,
                                    colData = coldata,
                                    design= ~ batchLabel + RNAi)
      
      # filtering 
      keep <- rowSums(counts(dds)>=10) >= 1
      dds <- dds[keep,]
      
      # library("BiocParallel")
      # register(MulticoreParam(2))
      dds <- DESeq(dds,minReplicatesForReplace=Inf)
      # resultsNames(dds) # lists the coefficients
      # plotDispEsts(dds)
      # assays(dds)[["cooks"]] #cooksCutoff in result function 
      # par(mar=c(8,5,2,2))
      # boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      
      #res <- results(dds, name = 'RNAi_nhr_20_vs_vector',alpha = 0.05)
      # or to shrink log fold changes association with condition:
      # res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
      # Note: We have sped up the apeglm method so it takes roughly about the same amount of time as normal, e.g. ~5 seconds for the pasilla dataset of ~10,000 genes and 7 samples. If fast shrinkage estimation of LFC is needed, but the posterior standard deviation is not needed, setting apeMethod="nbinomC" will produce a ~10x speedup, but the lfcSE column will be returned with NA. A variant of this fast method, apeMethod="nbinomC*" includes random starts.
      
      # write result 
      #list the DE for each individual RNAi 
      for (j in 2:(length(levels(RNAi_subset)))){
        targetGene = levels(RNAi_subset)[j]
        mySample = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
        # supply cooks filter metric (maximum cook for all samples related to tested covariate)
        dds.filt = dds
        
        # update base mean to the basemean between RNA and control
        # we use DEseq2's default setting to optimize power, but we define the base mean by only vectors plus the RNAi in query
        libWiseBaseMean = mcols(dds.filt)$baseMean
        mycounts = counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)]
        myconditions = dds.filt$RNAi[dds.filt$RNAi %in% c('x.vector',targetGene)] 
        myWeights = rep(0,length(myconditions))
        myWeights[myconditions == 'x.vector'] = 0.5 / sum(myconditions == 'x.vector') 
        myWeights[myconditions == targetGene] = 0.5 / sum(myconditions == targetGene)
        
        baseMean = apply(mycounts, 1, result <- function(x){result = weighted.mean(x,myWeights)})
        # baseVar = rowVars(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)])
        mcols(dds.filt)$baseMean = baseMean
        mcols(dds.filt)$baseVar = NA
        
        res <- results(dds.filt,name = str_replace_all(mySample,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
        stat = res$stat
        logFC_raw = res$log2FoldChange
        res <- lfcShrink(dds.filt, coef = str_replace_all(mySample,'-','.'), type="apeglm",apeMethod="nbinomC", res=res)
        res$stat = stat
        res$log2FoldChange_raw = logFC_raw
        
        # attach the counts information
        # filter by mean read counts of each component of contrast(will do in cutoff step)
        # batch effect is causing some trouble here, so we use median instead of mean; for met lib, we want to use mean 
        # ==> change to median alway to ensure robustness
        mVals1 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
        names(mVals1) = rownames(dds.filt)
        mVals2 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'x.vector'])
        names(mVals2) = rownames(dds.filt)
        res = as.data.frame(res)
        res$medianCount_RNAi = mVals1[rownames(res)]
        res$medianCount_ctr = mVals2[rownames(res)]
        res$libWiseBaseMean = libWiseBaseMean # this will be used in the MERGE modeling! (otherwise the pseudocounts will not be comparable across conditions)
        # save
        write.csv(res,paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'.csv',sep = ''))
        # resOrdered <- res[order(res$padj),]
        # outputTbl = as.data.frame(subset(res, padj < 0.05))
        # if (nrow(outputTbl) > 0){
        #   outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
        #   rownames(outputTbl) = 1:nrow(outputTbl)
        #   outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
        # }else{
        #   outputTbl = data.frame(WBID = 'NoHit')
        #   outputTbl$baseMean = NA
        #   outputTbl$log2FoldChange = NA
        #   outputTbl$lfcSE = NA
        #   outputTbl$pvalue = NA
        #   outputTbl$padj = NA
        #   outputTbl$medianCount_RNAi = NA
        #   outputTbl$medianCount_ctr = NA
        #   outputTbl$RNAi = targetGene
        #   outputTbl$libWiseBaseMean = NA
        #   rownames(outputTbl) = 1
        # }
        # outputTbl$batchID = batches[i]
        # mergeTbl = rbind(mergeTbl, outputTbl)
        print(mySample)
      }
    }else if(batches[i] == 'met2_lib5'){
      # pick the control samples (in the batches belonging to target DE group)
      subsetInd = batchLabel == batches[i]
      # merge selected library with the extra controls
      input_subset = input[,subsetInd] # %>% rownames_to_column('gene') %>%
      # inner_join(input_extr_ctr %>% rownames_to_column('gene'), by = 'gene')
      # rownames(input_subset) = input_subset$gene
      # input_subset = input_subset[,-1]
      batchLabel_subset = colnames(input_subset)
      batchLabel_subset = str_extract(batchLabel_subset,'_rep._met[0-9]+_')
      batchLabel_subset = str_replace(batchLabel_subset,'_met._$','')
      batchLabel_subset = as.factor(batchLabel_subset)
      # batchLabel_subset = as.factor(c(as.character(batchLabel[subsetInd]), batchLabel_extr_ctr))
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      # RNAi_subset[(1+length(RNAi_subset)):(length(batchLabel_extr_ctr)+length(RNAi_subset))] = 'x.vector'
      
      coldata =  data.frame(batchLabel = batchLabel_subset, RNAi = RNAi_subset)
      rownames(coldata) = colnames(input_subset)
      
      # add back in the extra replicate 
      loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        newdata = list(input, batchLabel, RNAi)
        return(newdata)
      }
      extra_data <- loadRData('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_extra_met2_lib5.Rdata')
      # merge 
      rep3_input = extra_data[[1]]
      rep3_RNAi = extra_data[[3]]
      keepind = str_detect(colnames(rep3_input),'_rep3_extra_met2_lib5$')
      rep3_input = rep3_input[,keepind]
      rep3_RNAi = rep3_RNAi[keepind]
      
      keepind = !str_detect(colnames(input_subset),'_rep3_met2_lib5$')
      input_subset = input_subset[,keepind]
      coldata = coldata[keepind,]
      
      input_subset <- input_subset %>% rownames_to_column('gene') %>%
        full_join(rep3_input %>% rownames_to_column('gene'), by = 'gene')
      rownames(input_subset) = input_subset$gene
      input_subset = input_subset[,-1]
      input_subset[is.na(input_subset)] = 0
      coldata2 = data.frame(row.names = colnames(rep3_input))
      coldata2$batchLabel = '_rep3'
      coldata2$RNAi = as.character(rep3_RNAi)
      coldata = rbind(coldata, coldata2)
      
      dds <- DESeqDataSetFromMatrix(countData = input_subset,
                                    colData = coldata,
                                    design= ~ batchLabel + RNAi)
      
      # filtering 
      keep <- rowSums(counts(dds)>=10) >= 1
      dds <- dds[keep,]
      
      # library("BiocParallel")
      # register(MulticoreParam(2))
      dds <- DESeq(dds,minReplicatesForReplace=Inf)
      # resultsNames(dds) # lists the coefficients
      # plotDispEsts(dds)
      # assays(dds)[["cooks"]] #cooksCutoff in result function 
      # par(mar=c(8,5,2,2))
      # boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      
      #res <- results(dds, name = 'RNAi_nhr_20_vs_vector',alpha = 0.05)
      # or to shrink log fold changes association with condition:
      # res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
      # Note: We have sped up the apeglm method so it takes roughly about the same amount of time as normal, e.g. ~5 seconds for the pasilla dataset of ~10,000 genes and 7 samples. If fast shrinkage estimation of LFC is needed, but the posterior standard deviation is not needed, setting apeMethod="nbinomC" will produce a ~10x speedup, but the lfcSE column will be returned with NA. A variant of this fast method, apeMethod="nbinomC*" includes random starts.
      
      # write result 
      #list the DE for each individual RNAi 
      for (j in 2:(length(levels(RNAi_subset)))){
        targetGene = levels(RNAi_subset)[j]
        mySample = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
        # supply cooks filter metric (maximum cook for all samples related to tested covariate)
        dds.filt = dds
        
        # update base mean to the basemean between RNA and control
        # we use DEseq2's default setting to optimize power, but we define the base mean by only vectors plus the RNAi in query
        libWiseBaseMean = mcols(dds.filt)$baseMean
        mycounts = counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)]
        myconditions = dds.filt$RNAi[dds.filt$RNAi %in% c('x.vector',targetGene)] 
        myWeights = rep(0,length(myconditions))
        myWeights[myconditions == 'x.vector'] = 0.5 / sum(myconditions == 'x.vector') 
        myWeights[myconditions == targetGene] = 0.5 / sum(myconditions == targetGene)
        
        baseMean = apply(mycounts, 1, result <- function(x){result = weighted.mean(x,myWeights)})
        # baseVar = rowVars(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)])
        mcols(dds.filt)$baseMean = baseMean
        mcols(dds.filt)$baseVar = NA
        
        res <- results(dds.filt,name = str_replace_all(mySample,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
        stat = res$stat
        logFC_raw = res$log2FoldChange
        res <- lfcShrink(dds.filt, coef = str_replace_all(mySample,'-','.'), type="apeglm",apeMethod="nbinomC", res=res)
        res$stat = stat
        res$log2FoldChange_raw = logFC_raw

        # attach the counts information
        # filter by mean read counts of each component of contrast(will do in cutoff step)
        # batch effect is causing some trouble here, so we use median instead of mean; for met lib, we want to use mean 
        # ==> change to median alway to ensure robustness
        mVals1 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
        names(mVals1) = rownames(dds.filt)
        mVals2 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'x.vector'])
        names(mVals2) = rownames(dds.filt)
        res = as.data.frame(res)
        res$medianCount_RNAi = mVals1[rownames(res)]
        res$medianCount_ctr = mVals2[rownames(res)]
        res$libWiseBaseMean = libWiseBaseMean # this will be used in the MERGE modeling! (otherwise the pseudocounts will not be comparable across conditions)
        # save
        write.csv(res,paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'.csv',sep = ''))
        # resOrdered <- res[order(res$padj),]
        # outputTbl = as.data.frame(subset(res, padj < 0.05))
        # if (nrow(outputTbl) > 0){
        #   outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
        #   rownames(outputTbl) = 1:nrow(outputTbl)
        #   outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
        # }else{
        #   outputTbl = data.frame(WBID = 'NoHit')
        #   outputTbl$baseMean = NA
        #   outputTbl$log2FoldChange = NA
        #   outputTbl$lfcSE = NA
        #   outputTbl$pvalue = NA
        #   outputTbl$padj = NA
        #   outputTbl$medianCount_RNAi = NA
        #   outputTbl$medianCount_ctr = NA
        #   outputTbl$RNAi = targetGene
        #   outputTbl$libWiseBaseMean = NA
        #   rownames(outputTbl) = 1
        # }
        # outputTbl$batchID = batches[i]
        # mergeTbl = rbind(mergeTbl, outputTbl)
        print(mySample)
      }
    }
  }
  #mergeTbl = mergeTbl[order(mergeTbl["RNAi"], mergeTbl["log2FoldChange"]),]
  #write.csv(mergeTbl,file = paste('output/DE_master_table_FDR005_',plateID,'.csv',sep = ''))
}

writeLines(capture.output(sessionInfo()), "sessionInfo_DEseq2_one2one.txt")