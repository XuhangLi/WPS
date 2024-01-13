RNAi_QC <- function(expID, libID, vectorSamples, status, ctr_bam_depth){
  cat('Executing QC pipeline for ',expID, '-',libID,'...\n',sep = '')
  
  # load annotation tables 
  WBID = read.table('./../data/metaData/WBIDtbl.txt',header = T,sep = '\t')
  
  # part I: process the count data 
  
  # ribosomal genes 
  riboGenes = read.table('./../data/metaData/ribosomalGenes.txt',header = T,sep = '\t')
  # sample sheet
  sampleSheet = read.csv(paste('./../data/sampleSheets_',status,'/sampleSheet_',expID,'_',libID,'.txt',sep = ''),header = F)
  
  if (!dir.exists('figures')){
    dir.create('figures')
  }
  if (!dir.exists('outputs')){
    dir.create('outputs')
  }
  
  # in case there is any existing connection
  graphics.off()
  
  # open a pdf for main QC output figures 
  pdf(paste("figures/",expID,'_',libID,'_QC_report.pdf',sep = '')) 
  
  # load read counts
  readsCount = read.table(paste("./../step1_process_raw_data/output/readsCounts/",expID,'-',libID,'_reads_count.txt',sep = ''),header = T,row.names = 1,sep = '\t')
  colnames(readsCount) = sampleSheet$V2[match(colnames(readsCount), sampleSheet$V1)]
  if (sum(is.na(colnames(readsCount))) >0){ # remove empty barcodes 
    cat('empty barcodes contain ',sum(readsCount[,is.na(colnames(readsCount))])/1e6,' mil reads! High number (e.g., > 1 mil) indicates a wrong/swapped barcode!\n',sep = '')
    readsCount = readsCount[,-which(is.na(colnames(readsCount)))]
  }
  
  # remove ribosomal RNA
  readsCount = readsCount[-which(is.element(rownames(readsCount),riboGenes$WormBase.Gene.ID)),]
  
  #plot the depth distribution
  depth = colSums(readsCount)
  depth2  <- depth[order(depth,decreasing=TRUE)]
  depth2 = depth2 / 1000000
  cat(as.character(sum(depth2 < 1)),' samples are low-depth (< 1 mil)\n',sep = '')
  lowSamples = names(depth2[depth2 < 3])
  
  hist(depth2,
       xlab = paste('Final depth distribution across ', length(depth2),' samples',sep = ''),
       ylab = 'Frequency',
       main = 'Depth distribution')
  text(x = 6, y = 10, labels = paste(as.integer(sum(depth2 > 3)/length(depth2)*100),'% samples have depth > 3 mil'), cex = 1.2)
  text(x = 6, y = 7.5, labels = paste(as.integer(sum(depth2 > 5)/length(depth2)*100),'% samples have depth > 5 mil'), cex = 1.2)

  # correct for dsRNA contamination
  # ignore the met1_lib1 suffix in the variable name. Not meaningful!
  ToCorrect_met1_lib = read.table(paste('./../step1_process_raw_data/output/RNAiAnalysis_',status,'/',
                                        expID,'-',libID,'_sorted_corrected_uniqOnly__',expID,'-',libID,'_sorted.bam.csv',sep = ''),header = T,row.names = 1,sep = ',')
  colnames(ToCorrect_met1_lib) = sampleSheet$V2[match(colnames(ToCorrect_met1_lib), sampleSheet$V1)]
  # remove NAs (they are empty barcodes)
  if (sum(is.na(colnames(ToCorrect_met1_lib))) >0){
    ToCorrect_met1_lib = ToCorrect_met1_lib[,-which(is.na(colnames(ToCorrect_met1_lib)))]
  }
  ToCorrect_met1_lib = ToCorrect_met1_lib[,colnames(readsCount)]
  # Filtering: avoid correction if it interfere with the power too much (control reads reduction >50%)
  vectInd = is.element(colnames(ToCorrect_met1_lib),vectorSamples)
  notCorrectedRNAi = character()
  for (i in 1:nrow(ToCorrect_met1_lib)){
    if (any(rownames(ToCorrect_met1_lib)[i] == rownames(readsCount))){
      if (rowMeans(as.matrix(ToCorrect_met1_lib[i,vectInd])) < rowMeans(as.matrix(readsCount[which(rownames(readsCount)==rownames(ToCorrect_met1_lib)[i]),vectInd])) * 0.5){
        cat('\033[31mWarining: dsRNA contamination in gene',rownames(ToCorrect_met1_lib)[i],'cannot be corrected because of loss of power!\033[39m\n')
        notCorrectedRNAi = c(notCorrectedRNAi,rownames(ToCorrect_met1_lib)[i])
      }else{
        readsCount[which(rownames(readsCount)==rownames(ToCorrect_met1_lib)[i]),] = ToCorrect_met1_lib[i,]
      }
    }
  } 
  cat('\033[31mWarining: The read counts of ',as.character(round(100*length(notCorrectedRNAi)/nrow(ToCorrect_met1_lib))),'% RNAi targeted genes cannot be corrected because of loss of power!\033[39m\n',sep = '')
  
  data = readsCount
  # normalize to TPM
  TPM = sweep(data, 2, colSums(data), FUN="/") * 1e6
  # save the processed count data 
  save(TPM,readsCount,notCorrectedRNAi, file=paste("./outputs/processed_read_count_",expID,'_',libID,".RData",sep = ''))
  
  
  # part II: run the RNAi QC analysis 
  
  # check the knockdown efficiency 
  # extract the target gene WBID from the sample names
  sampleNames = colnames(TPM)
  # remove control samples 
  sampleNames = setdiff(sampleNames, vectorSamples)
  sampleNames = str_replace(sampleNames,'^x.','') # remove prefix
  sampleNames = str_replace(sampleNames,'_rep.$','') # remove replicate suffix
  sampleNames = str_replace(sampleNames,'_L.$','') # remove seeding stage suffix (if any)
  sampleNames = str_replace(sampleNames,'_','-') # replace '_' with '-'
  sampleNames = unique(sampleNames)
  # skip incorrect RNAi or low depth samples
  sampleNames = sampleNames[!str_detect(sampleNames,'^LOWDEPTH-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^SHORT-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^VECTORLIKE-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^RCBVECTOR-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^MULTIPLE-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^NOSIGNAL-')]
  sampleNames = sampleNames[!str_detect(sampleNames,'^REALVECTOR-')]
  sampleWBID = WBID$WormBase.Gene.ID[match(sampleNames,WBID$Public.Name)]
  if(sum(!is.na(sampleWBID)) != length(sampleNames)){
    stop(paste('Cannot identify WBID for condition name', sampleNames[is.na(sampleWBID)]))
  }
  
  # calculate the target gene decrease by TPM 
  # 01112024: the codes for this section is very stupid - we assumed up to six replicates for one condition and 
  # try to match each RNAi with the vector sample from the same replicate and calculate the TPM fold change. This 
  # will be updated to a more robust and general algorithm in the future, but seems good for now...
  
  RNAitbl = data.frame(vector1 = numeric(length(sampleWBID)),vector2 = numeric(length(sampleWBID)),vector3 = numeric(length(sampleWBID)),
                       vector4 = numeric(length(sampleWBID)),vector5 = numeric(length(sampleWBID)),vector6 = numeric(length(sampleWBID)),
                       RNAi1 = numeric(length(sampleWBID)),RNAi2 = numeric(length(sampleWBID)),RNAi3 = numeric(length(sampleWBID)),
                       RNAi4 = numeric(length(sampleWBID)),RNAi5 = numeric(length(sampleWBID)),RNAi6 = numeric(length(sampleWBID)))
  rownames(RNAitbl) = sampleNames
  for (i in 1:length(sampleWBID)){
    target = sampleWBID[i]
    tmp = as.matrix(TPM[as.character(target),])
    tmp = as.data.frame(tmp[,order(tmp,decreasing = T)])
    tmp$names = rownames(tmp)
    colnames(tmp) = c('gene','names')
    
    # get the name (some gene may have 'L2' labels)
    name = sampleNames[i]
    name = unique(str_replace(rownames(tmp)[grep(paste('x.',str_replace(name,'-','_'),'_',sep = ''), rownames(tmp))],'_rep.$',''))
    if (length(name) > 1){
      if (length(unique(str_remove(name,'_L2$')))==1){
        # 02222022: fix a bug when both L2 sample and L1 sample appears in the same library
        # we do RNAi QC for them together 
        name = unique(str_remove(name,'_L2$'))
        rownames(tmp)[rownames(tmp)==paste(name,'_L2_rep1',sep = '')] = paste(name,'_rep4',sep = '')
        rownames(tmp)[rownames(tmp)==paste(name,'_L2_rep2',sep = '')] = paste(name,'_rep5',sep = '')
        rownames(tmp)[rownames(tmp)==paste(name,'_L2_rep3',sep = '')] = paste(name,'_rep6',sep = '')
        
      }else{
        stop('Failed to match gene ID with sample name. Check!')
      }
    }
   
    # get control read count for filtering 
    if(any(rownames(readsCount) == target)){
      RNAitbl$ctr_readsCount[i] = median(as.matrix(readsCount[target, colnames(readsCount) %in% vectorSamples]))
    }else{
      RNAitbl$ctr_readsCount[i] = 0
      cat('Notice: in ',expID,'-',libID,', RNAi targeted gene ', target,' is not detected!\n',sep = '')
    }
    
    # calculate FC within each replicate (these codes are stupid..)
    if (any(name == rownames(tmp))){# no replicates (if there is no replicate suffix in the name)
      RNAitbl$RNAi1[i] = tmp[paste('x.',str_replace(name,'-','_'),sep = ''),1]
      RNAitbl$RNAi2[i] = NA
      RNAitbl$RNAi3[i] = NA
      RNAitbl$RNAi4[i] = NA
      RNAitbl$RNAi5[i] = NA
      RNAitbl$RNAi6[i] = NA
      RNAitbl$vector1[i] = mean(tmp[rownames(tmp) %in% vectorSamples,1])
      RNAitbl$vector2[i] = NA
      RNAitbl$vector3[i] = NA
      RNAitbl$vector4[i] = NA
      RNAitbl$vector5[i] = NA
      RNAitbl$vector6[i] = NA
    }else{
      if (any(paste(name,'_rep1',sep = '') == rownames(tmp))){
        RNAitbl$RNAi1[i] = tmp[paste(name,'_rep1',sep = ''),1]
        RNAitbl$vector1[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep1$')],1])
      }else{
        RNAitbl$RNAi1[i] = NA
        RNAitbl$vector1[i] = NA
      }
      if (any(paste(name,'_rep2',sep = '') == rownames(tmp))){
        RNAitbl$RNAi2[i] = tmp[paste(name,'_rep2',sep = ''),1]
        RNAitbl$vector2[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep2$')],1])
      }else{
        RNAitbl$RNAi2[i] = NA
        RNAitbl$vector2[i] = NA
      }    
      if (any(paste(name,'_rep3',sep = '') == rownames(tmp))){
        RNAitbl$RNAi3[i] = tmp[paste(name,'_rep3',sep = ''),1]
        RNAitbl$vector3[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep3$')],1])
      }else{
        RNAitbl$RNAi3[i] = NA
        RNAitbl$vector3[i] = NA
      } 
      if (any(paste(name,'_rep4',sep = '') == rownames(tmp))){
        RNAitbl$RNAi4[i] = tmp[paste(name,'_rep4',sep = ''),1]
        RNAitbl$vector4[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep4$')],1])
      }else{
        RNAitbl$RNAi4[i] = NA
        RNAitbl$vector4[i] = NA
      } 
      if (any(paste(name,'_rep5',sep = '') == rownames(tmp))){
        RNAitbl$RNAi5[i] = tmp[paste(name,'_rep5',sep = ''),1]
        RNAitbl$vector5[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep5$')],1])
      }else{
        RNAitbl$RNAi5[i] = NA
        RNAitbl$vector5[i] = NA
      } 
      if (any(paste(name,'_rep6',sep = '') == rownames(tmp))){
        RNAitbl$RNAi6[i] = tmp[paste(name,'_rep6',sep = ''),1]
        RNAitbl$vector6[i] = mean(tmp[rownames(tmp) %in% vectorSamples[str_detect(vectorSamples,'_rep6$')],1])
      }else{
        RNAitbl$RNAi6[i] = NA
        RNAitbl$vector6[i] = NA
      } 
    }
  }
  
  
  # calcualte the LOG2FC ...
  RNAitbl$logFC1 = log2(RNAitbl$RNAi1 / RNAitbl$vector1)
  RNAitbl$logFC2 = log2(RNAitbl$RNAi2 / RNAitbl$vector2)
  RNAitbl$logFC3 = log2(RNAitbl$RNAi3 / RNAitbl$vector3)
  RNAitbl$logFC4 = log2(RNAitbl$RNAi4 / RNAitbl$vector4)
  RNAitbl$logFC5 = log2(RNAitbl$RNAi5 / RNAitbl$vector5)
  RNAitbl$logFC6 = log2(RNAitbl$RNAi6 / RNAitbl$vector6)
  # the code is too stupid to write comment... this is basically to avoid inf caused by deviding zero
  maxFC = min(c(RNAitbl$logFC1[!is.infinite(RNAitbl$logFC1)],
                RNAitbl$logFC2[!is.infinite(RNAitbl$logFC2)],
                RNAitbl$logFC3[!is.infinite(RNAitbl$logFC3)],
                RNAitbl$logFC4[!is.infinite(RNAitbl$logFC4)],
                RNAitbl$logFC5[!is.infinite(RNAitbl$logFC5)],
                RNAitbl$logFC6[!is.infinite(RNAitbl$logFC6)]),na.rm = T)
  RNAitbl$logFC1[is.infinite(RNAitbl$logFC1) & RNAitbl$logFC1<0] = maxFC
  RNAitbl$logFC2[is.infinite(RNAitbl$logFC2) & RNAitbl$logFC2<0] = maxFC
  RNAitbl$logFC3[is.infinite(RNAitbl$logFC3) & RNAitbl$logFC3<0] = maxFC
  RNAitbl$logFC4[is.infinite(RNAitbl$logFC4) & RNAitbl$logFC4<0] = maxFC
  RNAitbl$logFC5[is.infinite(RNAitbl$logFC5) & RNAitbl$logFC5<0] = maxFC
  RNAitbl$logFC6[is.infinite(RNAitbl$logFC6) & RNAitbl$logFC6<0] = maxFC
  
  # remove any FC when the control expression is less than 10 read count 
  RNAitbl[RNAitbl$ctr_readsCount < 10,] = NA #detection cutoff
  
  # make the dataframe for plotting
  RNAitbl$names = rownames(RNAitbl)
  df = rbind(RNAitbl[,c('logFC1','names')] %>% dplyr::rename(logFC = logFC1),
             RNAitbl[,c('logFC2','names')] %>% dplyr::rename(logFC = logFC2),
             RNAitbl[,c('logFC3','names')] %>% dplyr::rename(logFC = logFC3),
             RNAitbl[,c('logFC4','names')] %>% dplyr::rename(logFC = logFC4),
             RNAitbl[,c('logFC5','names')] %>% dplyr::rename(logFC = logFC5),
             RNAitbl[,c('logFC6','names')] %>% dplyr::rename(logFC = logFC6))
  
  df.summary <- df %>%
    group_by(names) %>%
    summarise(
      sd = sd(logFC, na.rm = TRUE),
      logFC = mean(logFC,na.rm = TRUE)
    )
  df.summary
  
  # label the uncorrected samples (filtered by 50% read reduction threshold)
  df$NotCorrected = is.element(rownames(df),WBID$Public.Name[match(notCorrectedRNAi,WBID$WormBase.Gene.ID)])
  df.summary$NotCorrected = is.element(df.summary$names,WBID$Public.Name[match(notCorrectedRNAi,WBID$WormBase.Gene.ID)])
  RNAitbl$NotCorrected = is.element(rownames(RNAitbl),WBID$Public.Name[match(notCorrectedRNAi,WBID$WormBase.Gene.ID)])
  write.csv(RNAitbl,file = paste('outputs/RNAi_efficiency_',expID,'_',libID,'.csv',sep = ''))
  
  # plot
  p <- ggplot(df, aes(x= reorder(names,logFC), y=logFC)) +
    geom_bar(stat = "identity", data = df.summary,
             fill = 'grey', color = "black") +
    geom_text(data = df.summary,aes(x= reorder(names,logFC), y=0.2, label = ifelse(NotCorrected, "#", "")),
              size = 5, color = 'red') +
    geom_jitter( position = position_jitter(0.2),
                 color = "black") + 
    geom_errorbar(aes(ymin = logFC-sd, ymax = logFC+sd),data = df.summary, width = 0.2) +
    ggtitle('RNAi efficiency') +
    coord_flip()+
    xlab('') + 
    theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8,face="bold"),
          panel.grid.major = element_line(colour = "grey",linetype = "dashed"),
          panel.background = element_blank())+
    geom_hline(yintercept = -1,linetype="dotted", 
               color = "red")
  print(p)
  
  
  # next, plot the detection of dsRNA signals (as FC against background)
  dsRNA = read.csv(paste('./../step1_process_raw_data/output/RNAiAnalysis_',status,'/',expID,'-',libID,'_sorted_corrected_uniqOnly__dsRNA_count.csv',sep = ''),header = F)
  dsRNA_ctr = read.csv(paste('./../step1_process_raw_data/output/RNAiAnalysis_',status,'/',expID,'-',libID,'_sorted_corrected_uniqOnly__dsRNA_count_ctr.csv',sep = ''),header = F)
  # add one pseudocount to the read count to avoid errors in FC calculation 
  dsRNA_ctr$V2 = dsRNA_ctr$V2+1 
  
  # make a variable for sample names without seeding stage info (this is to clean up some cases that 'L2' is labeled before 'rep1')
  # (just to match genes to sample names more easily)
  tmp = str_replace(colnames(readsCount),'_L[0-9]_','_') 
  tmp = str_replace_all(tmp,'-','_')
  # calculate dsRNA CPM
  for (i in 1:nrow(dsRNA)){
    # get the name (some gene may have 'L2' labels)
    name = WBID$Public.Name[match(dsRNA$V1[i],WBID$WormBase.Gene.ID)]
    name = unique(str_replace(tmp[grep(paste('(x.|_)',str_replace(name,'-','_'),'_',sep = ''), tmp)],'_rep.$',''))
    if (any(name == tmp)){# no replicates
      dsRNA$cpm[i] = dsRNA$V2[i]/sum(readsCount[,which(tmp %in% name)])*1e6
      dsRNA_ctr$cpm[i] = dsRNA_ctr$V2[i]/ctr_bam_depth*1e6
    }else{
      rnaiNames = paste(name,'_rep',1:6,sep = '') # for historical reason, we assume no more than 6 replicates 
      dsRNA$cpm[i] = dsRNA$V2[i]/sum(readsCount[,which(tmp %in% rnaiNames)])*1e6
      dsRNA_ctr$cpm[i] = dsRNA_ctr$V2[i]/ctr_bam_depth*1e6
    }
    
    rownames(dsRNA)[i] = WBID$Public.Name[which(WBID$WormBase.Gene.ID==dsRNA$V1[i])]
    rownames(dsRNA_ctr)[i] = WBID$Public.Name[which(WBID$WormBase.Gene.ID==dsRNA$V1[i])]
  }
  # calculate FCs 
  dsRNA$FC = log2(dsRNA$cpm / dsRNA_ctr$cpm)
  dsRNA$names = rownames(dsRNA)
  dsRNA$FC[is.infinite(dsRNA$FC) & dsRNA$FC < 0] = 0 # -inf fold change indicates dsRNA CPM in the RNAi condition is 0, so we just put the FC to 0 that indicate no enrichment of dsRNA
  dsRNA$FC[is.na(dsRNA$FC)] = 0 # nan fold change indicates dsRNA CPM in both control and RNAi condition is 0, so we just put the FC to 0 that indicate no enrichment of dsRNA
  
  # plot
  p<-ggplot(data=dsRNA, aes(x=reorder(names,FC), y=FC)) +
    geom_bar(stat="identity") +
    ggtitle('dsRNA detection') +
    coord_flip() +
    ylab('log2(CPM RNAi / CPM background)') +
    xlab('') + 
    theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8,face="bold"))+
    geom_hline(yintercept = 1,linetype="dotted", 
               color = "red")
  print(p)
  
  # summarize this QC result in text 
  # if average FC decrease in TPM is greater than 2 fold, it is a pass-QC conditions
  pass1 = RNAitbl$names[rowMeans(RNAitbl[,c('logFC1','logFC2','logFC3')],na.rm = TRUE) < -1] 
  pass1 = pass1[!is.na(pass1)] 
  # otherwise we require at least 10-fold enrichment of the dsRNA signals
  pass2 = rownames(dsRNA)[dsRNA$FC > log2(10)]
  cat('INFO: these RNAi conditions show high-confidence in RNAi QC (2-FC KD AND 10-FC detection of dsRNA): ',paste(intersect(pass1,pass2),collapse = ', '),'\n',sep = '')
  # RNAi failed
  cat('\033[31mWARNING: these RNAi conditions FAILED to pass RNAi QC (lack of 2-FC KD or 10-FC detection of dsRNA): ',paste(setdiff(sampleNames,union(pass1,pass2)),collapse = ', '),'\033[39m\n',sep = '')
  cat('\033[31mPlease check their QC figures carefully! \033[39m\n')
  fileConn<-file(paste("figures/",expID,'_',libID,'_summary.txt',sep = ''))
  txt = c()
  for (i in 1:length(depth2)){
    txt = c(txt,paste(names(depth2)[i],depth2[i],sep = '    '))
  }
  writeLines(c(paste('these RNAi conditions show high-confidence in RNAi QC (2-FC KD AND 10-FC detection of dsRNA): ',paste(intersect(pass1,pass2),collapse = ', '),sep = ''),
               paste('these RNAi conditions FAILED to pass RNAi QC (lack of 2-FC KD or 10-FC detection of dsRNA): ',paste(setdiff(sampleNames,union(pass1,pass2)),collapse = ', '),sep = ''),
               paste('depth of these samples is less than 3 mil: ',paste(lowSamples,collapse = ', '),sep = ''),
               'sample    depth',
               txt
               )
             , fileConn)
  close(fileConn)
  
  # finally, detect the cross-contamination of dsRNA by ploting the dsRNA heatmap
  dsRNA_met1 = read.csv(paste('./../step1_process_raw_data/output/RNAiAnalysis_',status,'/',expID,'-',libID,'_sorted_dsRNA__',expID,'-',libID,'_sorted.bam_dsRNAcount.csv',sep = ''),row.names = 1)
  colnames(dsRNA_met1) = sampleSheet$V2[match(colnames(dsRNA_met1), sampleSheet$V1)]
  if (sum(is.na(colnames(dsRNA_met1))) >0){
    dsRNA_met1 = dsRNA_met1[,-which(is.na(colnames(dsRNA_met1)))]
  }
  # skip the QC of bad samples that are already labeled out 
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.LOWDEPTH_')]
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.SHORT_')]
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.VECTORLIKE_')]
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.RCBVECTOR_')]
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.MULTIPLE_')]
  dsRNA_met1 = dsRNA_met1[,!str_detect(colnames(dsRNA_met1),'^x.NOSIGNAL_')]
  
  # calculate CPM
  rownames(dsRNA_met1) = paste(str_replace(WBID$Public.Name[match(rownames(dsRNA_met1), WBID$WormBase.Gene.ID)],'-','_'),sep = '')
  dsRNA_met1_cpm = sweep(dsRNA_met1, 2,colSums(readsCount[,colnames(dsRNA_met1)])/1e6,'/') 
  
  # control bam file cpm
  dsRNA_bg_met1 = read.csv(paste('./../step1_process_raw_data/output/RNAiAnalysis_',status,'/',expID,'-',libID,'_sorted_dsRNA__dsRNA_count_whole_ctrBAM.csv',sep = ''),row.names = 1,header = F)
  dsRNA_bg_met1_cpm = (dsRNA_bg_met1) / ctr_bam_depth * 1e6 
  rownames(dsRNA_bg_met1_cpm) = paste(str_replace(WBID$Public.Name[match(rownames(dsRNA_bg_met1_cpm), WBID$WormBase.Gene.ID)],'-','_'),sep = '')
  
  # calculate FC 
  # To increase the robustness in the visual interpretation of the heatmap, we applied pseudocount at the CPM level here. 
  # this reduces the backgrounds in the heatmap and makes the swapped conditions more easily identifiable
  dsRNA_FC = sweep(dsRNA_met1_cpm, 1,dsRNA_bg_met1_cpm[rownames(dsRNA_met1_cpm),],'/')
  dsRNA_FC_logTPM = sweep(log2(dsRNA_met1_cpm+1), 1,log2(dsRNA_bg_met1_cpm[rownames(dsRNA_met1_cpm),]+1),'-')
  
  # prepare the matrix for plotting
  sorted_mat = dsRNA_FC_logTPM
  # clean up the matrix 
  # first always keep what is in the library and put them at the beginning
  # then, we order the columns to align with the rows and put ctr behind the RNAi 
  colOrder = c(sort(setdiff(colnames(sorted_mat),vectorSamples)),vectorSamples)
  rowOrder = str_replace(str_replace(unique(str_replace(sort(setdiff(colnames(sorted_mat),vectorSamples)),'_rep.$','')),'^x.',''),'_L[0-9]$','')
  keepL = length(rowOrder)
  rowOrder = c(rowOrder,setdiff(rownames(sorted_mat),rowOrder))
  sorted_mat = sorted_mat[rowOrder,colOrder]
  
  # next, we remove genes that are not in the library and not show strong dsRNA signals
  # (only potentially target gene that a swapped condition targets to will be visualized)
  # here we filtered artificially at 2-fold change minimum 
  sorted_mat = rbind(sorted_mat[1:keepL,],sorted_mat[(keepL+1):nrow(sorted_mat),][apply(sorted_mat[(keepL+1):nrow(sorted_mat),],1,max) > log2(2),])
  
  maxVal = max((sorted_mat))
  seq1 = seq(0,maxVal,maxVal/100)
  pheatmap((sorted_mat),
           breaks = seq1,
           color = colorRampPalette((rev(brewer.pal(n = 7, name =
                                                      "RdBu"))))(2*length(seq1)+1)[(length(seq1)):(2*length(seq1)+1)],
           cluster_rows = FALSE,cluster_cols = FALSE,
           fontsize_row = 2,
           main = 'Log2(dsRNA_CPM_RNAi + 1) - log2(dsRNA_CPM_emptyLib + 1)')
  while (!is.null(dev.list())) {dev.off()}
}
nonRNAi_QC <- function(expID, libID, status){
  cat('Formating data for ',expID, '-',libID,'...\n',sep = '')
  
  # process the count data 
  
  # ribosomal genes 
  riboGenes = read.table('./../data/metaData/ribosomalGenes.txt',header = T,sep = '\t')
  # sample sheet
  sampleSheet = read.csv(paste('./../data/sampleSheets_',status,'/sampleSheet_',expID,'_',libID,'.txt',sep = ''),header = F)
  
  if (!dir.exists('figures')){
    dir.create('figures')
  }
  if (!dir.exists('outputs')){
    dir.create('outputs')
  }
  
  # in case there is any existing connection
  graphics.off()
  
  # open a pdf for main QC output figures 
  pdf(paste("figures/NON_RNAi_",expID,'_',libID,'_QC_report.pdf',sep = '')) 
  
  # load read counts
  readsCount = read.table(paste("./../step1_process_raw_data/output/readsCounts/",expID,'-',libID,'_reads_count.txt',sep = ''),header = T,row.names = 1,sep = '\t')
  colnames(readsCount) = sampleSheet$V2[match(colnames(readsCount), sampleSheet$V1)]
  if (sum(is.na(colnames(readsCount))) >0){ # remove empty barcodes 
    cat('empty barcodes contain ',sum(readsCount[,is.na(colnames(readsCount))])/1e6,' mil reads! High number (e.g., > 1 mil) indicates a wrong/swapped barcode!\n',sep = '')
    readsCount = readsCount[,-which(is.na(colnames(readsCount)))]
  }
  
  # remove ribosomal RNA
  readsCount = readsCount[-which(is.element(rownames(readsCount),riboGenes$WormBase.Gene.ID)),]
  
  #plot the depth distribution
  depth = colSums(readsCount)
  depth2  <- depth[order(depth,decreasing=TRUE)]
  depth2 = depth2 / 1000000
  cat(as.character(sum(depth2 < 1)),' samples are low-depth (< 1 mil)\n',sep = '')
  lowSamples = names(depth2[depth2 < 3])
  
  hist(depth2,
       xlab = paste('Final depth distribution across ', length(depth2),' samples',sep = ''),
       ylab = 'Frequency',
       main = 'Depth distribution')
  text(x = 6, y = 10, labels = paste(as.integer(sum(depth2 > 3)/length(depth2)*100),'% samples have depth > 3 mil'), cex = 1.2)
  text(x = 6, y = 7.5, labels = paste(as.integer(sum(depth2 > 5)/length(depth2)*100),'% samples have depth > 5 mil'), cex = 1.2)
  
  data = readsCount
  # normalize to TPM
  TPM = sweep(data, 2, colSums(data), FUN="/") * 1e6
  # save the processed count data 
  save(TPM,readsCount, file=paste("./outputs/NON_RNAi_processed_read_count_",expID,'_',libID,".RData",sep = ''))
  
}
EDA_plots <- function(expID, libID, controlID){
  
  ## load and merge all data
  load(paste("./outputs/processed_read_count_",expID,'_',libID,".RData",sep = ''))
  
  # remove the low depth in order not to bias normalization
  readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
  
  # construct a dds object and use DEseq2 package for EDA
  input <- readsCount
  coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
  rownames(coldata) = colnames(input)
  coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
  coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
  coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
  
  dds <- DESeqDataSetFromMatrix(countData = input,
                                colData = coldata,
                                design= ~ rep + RNAi)
  
  # filtering 
  keep <- rowSums(counts(dds)>=10) >= 1
  dds <- dds[keep,]
  
  # variance stabilization transformation 
  vsd <- vst(dds, blind=FALSE)
  mat <- assay(vsd)
  # remove any potential batch effect associated with replicates 
  mat <- limma::removeBatchEffect(mat, batch = vsd$rep, design= model.matrix(~ vsd$RNAi))
  
  if (controlID %in% unique(as.character(vsd$RNAi))){
    allRNAi = c(controlID,setdiff(unique(as.character(vsd$RNAi)),controlID))
  }else{
    allRNAi = unique(as.character(vsd$RNAi))
  }
  
  dataMat = as.data.frame(mat)
  pdf(paste("figures/",expID,'_',libID,'_initial_EDA.pdf',sep = '')) 
  par(mar = c(0, 0, 0, 0))
  for (i in 1:length(allRNAi)){
    if (controlID %in% allRNAi){
      mat_sub = mat[,dds$RNAi %in% c(controlID, allRNAi[i])]
    }else{
      mat_sub = mat[,dds$RNAi %in% c(allRNAi[i])]
    }
    pca <- prcomp(t(mat_sub))
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    percentVar <- round(100 * percentVar)
    pcaData = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat_sub))
    if (controlID %in% allRNAi){
      pcaData$RNAi = dds$RNAi[dds$RNAi %in% c(controlID,allRNAi[i])]
      pcaData$batch = dds$batchLabel[dds$RNAi %in% c(controlID,allRNAi[i])]
    }else{
      pcaData$RNAi = dds$RNAi[dds$RNAi %in% c(allRNAi[i])]
      pcaData$batch = dds$batchLabel[dds$RNAi %in% c(allRNAi[i])]
    }
    pcaData$RNAi = factor(as.character(pcaData$RNAi), levels = unique(c(allRNAi[i],controlID)))
    p1 = ggplot(pcaData, aes(PC1, PC2,label = name)) + #  or time
      geom_point(aes(color = RNAi,), size = 3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      scale_color_manual(values=c("red","grey")) +
      geom_text_repel()+
      ggtitle(paste(paste(expID,libID,i,allRNAi[i],sep = '-'),'\n',
                    paste(names(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]])),collapse = ' '),'\n',
                    paste(as.character(colSums(counts(dds)[,dds$RNAi %in% allRNAi[i]]) /1e6),collapse = ' '),sep = '')
      )
    print(p1)

    sampleDists <- dist(t(mat_sub))
    sampleDistMatrix <- as.matrix(sampleDists)
    library("RColorBrewer")
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,col=colors)
    
    mat_sub = as.data.frame(mat[,dds$RNAi %in% allRNAi[i]])
    vectorsamples = colnames(mat_sub)
    p2 = ggpairs(mat_sub,
                lower = list(continuous = wrap("points", alpha = 0.5,    size=0.1))) 
    # rasterize image to reduce file size 
    width_in_inches <- 7
    height_in_inches <- 7
    dpi <- 300
    png(filename = "temp_plot.png", width = width_in_inches * dpi, height = height_in_inches * dpi, res = dpi)
    print(p2)
    dev.off()
    plot_raster <- png::readPNG("temp_plot.png")
    plot.new()
    rasterImage(plot_raster, 0, 0, 1, 1)
    unlink("temp_plot.png")
    
  }
  while (!is.null(dev.list())) {dev.off()}
}



qLevels <- function(expID, libID, badSet){
  load(paste("./outputs/processed_read_count_",expID,'_',libID,".RData",sep = ''))
  
  # remove the low depth in order not to bias normalization
  readsCount = readsCount[, !str_detect(colnames(readsCount),'^x.LOWDEPTH_')]
  
  # create a dds object for easy operation
  input <- readsCount
  coldata =  data.frame(sampleName = as.factor(colnames(readsCount)))
  rownames(coldata) = colnames(input)
  coldata$RNAi = str_replace(coldata$sampleName,'_rep.$','')
  coldata$RNAi = as.factor(str_replace(coldata$RNAi,'_well.$',''))
  coldata$rep = as.factor(str_extract(coldata$sampleName,'_rep.$'))
  
  invisible(suppressMessages(dds <- DESeqDataSetFromMatrix(countData = input,
                                colData = coldata,
                                design= ~ rep + RNAi)))
  
  # filtering 
  keep <- rowSums(counts(dds)>=10) >= 1
  dds <- dds[keep,]
  
  # ignore the number conversions - just for historical reason - it is basically bad and good 
  dds$Qlevel = 1
  dds$Qlevel[colnames(dds) %in% badSet] = 4
  if (!all(badSet %in% colnames(dds))){
    stop('The input sample ID cannot be matched to data. Check for typos!')
  }
 
  # save the quality levels
  qLevels =  data.frame(sampleName = colnames(dds), QL = dds$Qlevel)
  rownames(qLevels) = qLevels$sampleName
  qLevels$QL[qLevels$QL == 1] = 'good'
  qLevels$QL[qLevels$QL == 4] = 'bad'
  write.csv(qLevels,paste('outputs/qualityLevels_',expID,'_', libID,'.csv',sep = ''))
}

