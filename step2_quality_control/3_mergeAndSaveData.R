library(stringr)

# we merge dataset by plate for future use. 

# after merge all and clean bad, write out two data:
# master data matrix ==> batch-corrected TPM (log2 --> unscaled)
# master raw data ==> all merged readcounts; batch labeled and DE group labeled; one NHR will only appears in one DE group. The same vector rep will not be in the same DE group


########################################################################
library(edgeR)
library(dplyr)
library(tibble)
cleanAndSave <- function(plateID, batchExcluded){
  load(paste("./outputs/TPM_preserved_",plateID,"_lib1.RData",sep = ''))
  batchLabel = rep(paste(plateID,'_lib1',sep = ''),ncol(readsCount))
  input = readsCount
  RNAi = colnames(input)
  RNAi = str_replace(RNAi,'_rep.$','')
  RNAi = str_replace(RNAi,'_well.$','')
  colnames(input) = paste(colnames(input),'_',plateID,'_lib1',sep = '')
  #assign the sample label order to put controls at begining in each batch
  ctrInd = which(str_detect(RNAi,'^x.vector'))
  input = input[,c(ctrInd,setdiff(1:ncol(input),ctrInd))]
  batchLabel = batchLabel[c(ctrInd,setdiff(1:ncol(input),ctrInd))]
  RNAi = RNAi[c(ctrInd,setdiff(1:ncol(input),ctrInd))]
  
  allBatches = paste(plateID,c('lib2','lib3','lib4','lib5','lib6'),sep = '_');
  allBatches = setdiff(allBatches,batchExcluded)
  
  # bad samples
  bad_list_certain = c()
  libExc = str_extract(batchExcluded,'lib.$')
  for(libID in setdiff(c('lib1','lib2','lib3','lib4','lib5','lib6'),libExc)){
    badsamples = read.csv(paste('outputs/qualityLevels_',plateID,'_',libID,'.csv',sep = ''),row.names = 1)
    badsamples = badsamples$sampleName[badsamples$QL == 'bad']
    if(length(badsamples)>0){
      bad_list_certain = c(bad_list_certain,
                           paste(badsamples,plateID,libID,sep = '_'))
    }
  }
  
  for (i in 1:length(allBatches)){
    load(paste("./outputs/TPM_preserved_",allBatches[i],".RData",sep = ''))
    RNAi0 = str_replace(colnames(readsCount),'_rep.$','')
    RNAi0 = str_replace(RNAi0,'_well.$','')
    #assign the sample label order to put controls at begining in each batch
    ctrInd = which(str_detect(RNAi0,'^x.vector'))
    readsCount = readsCount[,c(ctrInd,setdiff(1:ncol(readsCount),ctrInd))]
    RNAi0 = RNAi0[c(ctrInd,setdiff(1:ncol(readsCount),ctrInd))]
    RNAi = c(RNAi, RNAi0)
    
    colnames(readsCount) = paste(colnames(readsCount),'_',allBatches[i],sep = '')
    input <- input %>% rownames_to_column('gene') %>%
      full_join(readsCount %>% rownames_to_column('gene'), by = 'gene')
    rownames(input) = input$gene
    input = input[,-1]
    input[is.na(input)] = 0
    batchLabel = c(batchLabel, rep(allBatches[i],ncol(readsCount)))
  }
  
  # total sample number before any cleaning 
  N_total = ncol(input)
  N_bad = length(bad_list_certain)
  N_low = sum(str_detect(colnames(input),'^x.LOWDEPTH_'))
  N_SHORT = sum(str_detect(colnames(input),'^x.SHORT_'))
  N_VECTORLIKE = sum(str_detect(colnames(input),'^x.VECTORLIKE_'))
  N_RCBVECTOR = sum(str_detect(colnames(input),'^x.RCBVECTOR_'))
  N_MULTIPLE = sum(str_detect(colnames(input),'^x.MULTIPLE_'))
  N_NOSIGNAL = sum(str_detect(colnames(input),'^x.NOSIGNAL_'))
  print(paste('total sample: ',N_total,sep = ''))
  print(paste('total bad sample: ',N_bad,sep = ''))
  print(paste('total low depth sample: ',N_low,sep = ''))
  print(paste('total SHORT RNAi sample: ',N_SHORT,sep = ''))
  print(paste('total VECTORLIKE RNAi sample: ',N_VECTORLIKE,sep = ''))
  print(paste('total RCBVECTOR RNAi sample: ',N_RCBVECTOR,sep = ''))
  print(paste('total MULTIPLE RNAi sample: ',N_MULTIPLE,sep = ''))
  print(paste('total NOSIGNAL RNAi sample: ',N_NOSIGNAL,sep = ''))
  print(paste('total pass rate: ',(N_total - N_bad - N_low)/N_total,sep = ''))
  print(paste('clean pass rate: ',(N_total - N_bad - N_low - N_SHORT - N_VECTORLIKE - N_RCBVECTOR - N_MULTIPLE - N_NOSIGNAL)/N_total,sep = ''))
  # remove all the bad samples and low depth samples
  bad_list_certain = c(bad_list_certain,colnames(input)[str_detect(colnames(input),'^x.LOWDEPTH_')])
  keep = !is.element(colnames(input),bad_list_certain)
  input = input[,keep]
  RNAi = RNAi[keep]
  batchLabel = batchLabel[keep]
  
  if (ncol(input) != (N_total - N_bad - N_low)){
    stop('error occurred!')
  }
  
  # show replication number
  print(sort(table(RNAi)))
  # change metadata to factor
  RNAi = factor(RNAi,levels =c('x.vector',setdiff(unique(RNAi),'x.vector')))
  batchLabel = as.factor(batchLabel)
  
  # save data sets
  # save the merged and cleaned raw readscount labeled with DE group and batches 
  save('input','RNAi','batchLabel',file = paste('outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  
  # write the batch-corrected TPM tables
  repLabel = str_extract(colnames(input),'_rep._met._')
  repLabel = str_replace(repLabel,'_met._$','')
  batchLabel2 = paste(batchLabel,repLabel,sep = '')
  
  coldata =  data.frame(batchLabel = as.factor(batchLabel2), RNAi = RNAi)
  rownames(coldata) = colnames(input)
  
  # met3-lib4 has no vector, so model matrix not full rank; we treat it specially
  if ('met3_lib4_rep1' %in% coldata$batchLabel){
    met3_lib4_input = input[,coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3')]
    met3_lib4_coldata = coldata[coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3'),]
    met3_lib4_RNAi = levels(dropEmptyLevels(met3_lib4_coldata$RNAi))
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = met3_lib4_input,
                                  colData = met3_lib4_coldata,
                                  design= ~ RNAi)
    # filtering 
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]
    # save('dds',file = 'outputs/DEseqInput.Rdata') # we dont save it since we dont do DE on this design matrix 
    TPMmat = counts(dds, normalized = FALSE, replaced = FALSE) 
    TPMmat = sweep(TPMmat, 2, colSums(TPMmat), FUN="/") * 1e6
    mat <- TPMmat + 1
    write.csv(mat,file = paste('outputs/batchCorrectedTPM_met3_lib4.csv',sep = ''))
    # the mean and std table
    meanTbl = data.frame(row.names = rownames(mat))
    stdTbl = data.frame(row.names = rownames(mat))
    RNAiLevels = met3_lib4_RNAi
    for (i in 1:length(RNAiLevels)){
      meanTbl[,RNAiLevels[i]] = rowMeans(mat[,met3_lib4_coldata$RNAi %in% RNAiLevels[i]]) 
      stdTbl[,RNAiLevels[i]] = rowSds(mat[,met3_lib4_coldata$RNAi %in% RNAiLevels[i]]) 
    }
    IDtbl = data.frame(WBID = rownames(mat))
    WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
    IDtbl$GeneNames = WBID$Public.Name[match(IDtbl$WBID,WBID$WormBase.Gene.ID)]
    write.csv(meanTbl,file = paste('outputs/mean_TPM_met3_lib4.csv',sep = ''))
    write.csv(stdTbl,file = paste('outputs/std_TPM_met3_lib4.csv',sep = ''))
    write.csv(IDtbl,file = paste('outputs/geneID_TPM_met3_lib4.csv',sep = ''))
    
    input = input[,!(coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3'))]
    RNAi = dropEmptyLevels(RNAi[!(coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3'))])
    batchLabel = dropEmptyLevels(batchLabel[!(coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3'))])
    coldata = coldata[!(coldata$batchLabel %in% c('met3_lib4_rep1','met3_lib4_rep2','met3_lib4_rep3')),]
    allBatches = setdiff(allBatches,'met3_lib4')
  }
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = input,
                                colData = coldata,
                                design= ~ batchLabel + RNAi)
  # filtering 
  keep <- rowSums(counts(dds)>=10) >= 1
  dds <- dds[keep,]
  # save('dds',file = 'outputs/DEseqInput.Rdata') # we dont save it since we dont do DE on this design matrix 
  TPMmat = counts(dds, normalized = FALSE, replaced = FALSE) 
  TPMmat = sweep(TPMmat, 2, colSums(TPMmat), FUN="/") * 1e6
  mat <- log2(TPMmat + 1)
  mat <- limma::removeBatchEffect(mat, batch = dds$batchLabel, design= model.matrix(~ RNAi))
  mat <- 2^mat
  write.csv(mat,file = paste('outputs/batchCorrectedTPM_',plateID,'.csv',sep = ''))
  # the mean and std table
  meanTbl = data.frame(row.names = rownames(mat))
  stdTbl = data.frame(row.names = rownames(mat))
  RNAiLevels = setdiff(levels(RNAi),'x.vector')
  for (i in 1:length(RNAiLevels)){
    meanTbl[,RNAiLevels[i]] = rowMeans(mat[,RNAi %in% RNAiLevels[i]]) 
    stdTbl[,RNAiLevels[i]] = rowSds(mat[,RNAi %in% RNAiLevels[i]]) 
  }
  # find the deepest library to calculate for vectors
  a = sort(colSums(input)[grep('^x.vector',colnames(input))])/1e6
  aveDph = mean(a[grep('_lib1$',names(a))])
  for (i in 1:length(allBatches)){
    aveDph = c(aveDph, mean(a[grep(paste(allBatches[i],'$',sep = ''),names(a))]))
  }
  maxBatch = which(aveDph == max(aveDph))
  if (maxBatch == 1){
    meanTbl$vector = rowMeans(mat[,RNAi %in% 'x.vector' & batchLabel %in% paste(plateID,'_lib1',sep = '')]) 
    stdTbl$vector = rowSds(mat[,RNAi %in% 'x.vector' & batchLabel %in% paste(plateID,'_lib1',sep = '')]) 
  }else{
    meanTbl$vector = rowMeans(mat[,RNAi %in% 'x.vector' & batchLabel %in% allBatches[(maxBatch-1)]]) 
    stdTbl$vector = rowSds(mat[,RNAi %in% 'x.vector' & batchLabel %in% allBatches[(maxBatch-1)]]) 
  }
  
  IDtbl = data.frame(WBID = rownames(mat))
  WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
  IDtbl$GeneNames = WBID$Public.Name[match(IDtbl$WBID,WBID$WormBase.Gene.ID)]
  write.csv(meanTbl,file = paste('outputs/mean_TPM_',plateID,'.csv',sep = ''))
  write.csv(stdTbl,file = paste('outputs/std_TPM_',plateID,'.csv',sep = ''))
  write.csv(IDtbl,file = paste('outputs/geneID_TPM_',plateID,'.csv',sep = ''))
} 

plateID = 'met13'
#batchExcluded = c('met3_lib4','met9_lib4')
batchExcluded = c()

cleanAndSave(plateID, batchExcluded)

# summary dataset
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13','extra_met2_lib5')
RNAiConditions = c()
batchIDs = c()
condUID = c()
for (plateID in plateIDs){
  load(paste('outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''));
  for (batch in levels(batchLabel)){
    if (plateID != 'extra_met2_lib5'){
      RNAiConditions = c(RNAiConditions, unique(as.character(RNAi[batchLabel %in% batch])))
      batchIDs = c(batchIDs, rep(batch, length(unique(as.character(RNAi[batchLabel %in% batch])))))
    }
    condUID = c(condUID, paste(as.character(RNAi[batchLabel %in% batch]), 
                               rep(batch, length(as.character(RNAi[batchLabel %in% batch])))))
  }
  if (any(unique(levels(RNAi)) != unique(as.character(RNAi)))){stop('check dataset')}
}
# specially, we remove the batch label prefix for extra_met2_lib5
condUID = str_replace(condUID,'extra_met2_lib1$','met2_lib5')
repCounts = sort(table(condUID))
sum(repCounts < 2)
sum(repCounts == 2)
sum(repCounts == 3)
repCounts[repCounts>3]
write.csv(data.frame(RNAi_ID = RNAiConditions, RNAi_name = str_replace(RNAiConditions,'^x.',''), batch_ID = batchIDs),'outputs/RNAi_batch_lookup.csv')

RNAiConditions = table(RNAiConditions)
RNAiConditions[RNAiConditions>1]
summaryTbl = data.frame(RNAi_ID = names(RNAiConditions), independentExp = as.numeric(RNAiConditions))
summaryTbl$geneName = str_replace(summaryTbl$RNAi_ID, '^x.','')
summaryTbl = summaryTbl[!str_detect(summaryTbl$RNAi_ID, 'x.VECTORLIKE_'),]
summaryTbl = summaryTbl[!str_detect(summaryTbl$RNAi_ID, 'x.RCBVECTOR_'),]
summaryTbl = summaryTbl[!str_detect(summaryTbl$RNAi_ID, 'x.NOSIGNAL_'),]
summaryTbl = summaryTbl[!str_detect(summaryTbl$RNAi_ID, 'x.REALVECTOR_'),]
summaryTbl$geneName = str_replace(summaryTbl$geneName,'_L.$','')
summaryTbl$geneName = str_replace(summaryTbl$geneName, 'SHORT_','')
summaryTbl$geneName = str_replace(summaryTbl$geneName, 'MULTIPLE_','')
summaryTbl$geneName = str_replace_all(summaryTbl$geneName,'_','-')
WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
summaryTbl$WBID = WBID$WormBase.Gene.ID[match(summaryTbl$geneName,WBID$Public.Name)]
summaryTbl = summaryTbl[!(summaryTbl$RNAi_ID %in% 'x.vector'),]
metabolicGenes = read.csv('./../input_data/otherTbls/allMetGenes.txt')
iCELgenes = read.csv('./../input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
worcat = read.csv('./../input_data/otherTbls/wormcat_whole_genome_nov-16-2019.csv')
worcat$summary = NA
for (i in 1:nrow(worcat)){
  worcat$summary[i] = paste(worcat$Category.1[i],worcat$Category.2[i],worcat$Category.3[i],worcat$Automated.Description[i],sep = ' - ')
}
summaryTbl$isMetabolic = summaryTbl$WBID %in% metabolicGenes$WBName
summaryTbl$isICEL = summaryTbl$WBID %in% iCELgenes$WormBase.ID
summaryTbl$wormcat = worcat$summary[match(summaryTbl$WBID,worcat$Wormbase.ID)]
write.csv(summaryTbl,'outputs/RNAi_summary_final_dataset.csv')


# check iCEL coverage 
# covered iCEL genes
sum(iCELgenes$WormBase.ID %in% summaryTbl$WBID)

# covered non-dispensible gene
nondispensableGenes = read.csv('./../input_data/otherTbls/non_dispensable_gene_table.csv')
sum(iCELgenes$WormBase.ID[match(nondispensableGenes$nondis,iCELgenes$Gene.Name)] %in% summaryTbl$WBID)

# covered perturbale reactions
nondispensableRxns = read.csv('./../input_data/otherTbls/non_dispensable_gene_table_by_rxn.csv')
sum(!(nondispensableRxns$nondispensable_genes %in% c('','ND;','NA;','TBD;','Unknown;')))
nondispensableRxns$perturbed = F
for (i in 1:nrow(nondispensableRxns)){
  if (!(nondispensableRxns$nondispensable_genes[i] %in% c('','ND;','NA;','TBD;','Unknown;'))){
    genes = strsplit(nondispensableRxns$nondispensable_genes[i] ,';')[[1]]
    genes = iCELgenes$WormBase.ID[match(genes, iCELgenes$Gene.Name)]
    if (any(genes %in% summaryTbl$WBID)){
      nondispensableRxns$perturbed[i] = T
    }
  }
}
sum(nondispensableRxns$perturbed)


# finally, save the publishable data collection with standardized meta data sheet 
# first gather the data 
# summary dataset
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13','extra_met2_lib5')
RNAiCondition_ID = c()
batchID = c()
plateID = c()
sampleUID = c()
dataCollection = list()
for (myplateID in plateIDs){
  load(paste('outputs/cleaned_merged_raw_data_',myplateID,'.Rdata',sep = ''))
  dataCollection[[myplateID]] = input
  sampleUID = c(sampleUID, colnames(input))
  RNAiCondition_ID = c(RNAiCondition_ID, as.character(RNAi))
  batchID = c(batchID, as.character(batchLabel))
  plateID = c(plateID, rep(myplateID, ncol(input)))
}
mataTable = data.frame(sampleUID,RNAiCondition_ID,plateID, batchID)
# load the seeding stage annotation 
seedingStage = read.csv('./../input_data/metaData/special_seeding_stages.csv')
seedingStage$condID = paste('x.',seedingStage$RNAi_name,'_',seedingStage$library,sep = '')
setdiff(seedingStage$condID, paste(RNAiCondition_ID, batchID,sep = '_'))
mataTable$seedingStage = 'L1'
mataTable$seedingStage[paste(RNAiCondition_ID, batchID,sep = '_') %in% seedingStage$condID] = 
  seedingStage$seeding.stage[match(paste(RNAiCondition_ID, batchID,sep = '_')[paste(RNAiCondition_ID, batchID,sep = '_') %in% seedingStage$condID], seedingStage$condID)]

# also save the count data as the plain text matrix
input = dataCollection[[1]]
for (myplateID in plateIDs[2:length(plateIDs)]){
  input <- input %>% rownames_to_column('gene') %>%
    full_join(dataCollection[[myplateID]] %>% rownames_to_column('gene'), by = 'gene')
  rownames(input) = input$gene
  input = input[,-1]
  input[is.na(input)] = 0
}
#all(colnames(input) == sampleUID)
save('dataCollection','mataTable',file = paste('outputs/full_dataset.Rdata',sep = ''))
write.csv(input,file = paste('outputs/readCounts_full_dataset.csv',sep = ''))
write.csv(mataTable,file = paste('outputs/metadata_full_dataset.csv',sep = ''))

# write the averge TPM
readCountsFull = read.csv('outputs/readCounts_full_dataset.csv',row.names = 1)
TPM = sweep(readCountsFull, 2, colSums(readCountsFull), FUN="/") * 1e6
TPM_info = data.frame(WBID = rownames(readCountsFull))
TPM_info$whole_dataset_TPM = rowMeans(TPM)
TPM_info$vector_average_TPM = rowMeans(TPM[,grep('vector|REALVECTOR',colnames(TPM))])
write.csv(TPM_info,file = paste('outputs/WT_TPM.csv',sep = ''))


