import pysam
import argparse
import pandas as pd
import collections


def parseHeader(header):
    # extract the sequences that should be the UMI and BC:
    tmp1 = header.split(" ")
    tmp2 = tmp1[0].split(':')
    #use hard index to improve performance
    return tmp2[-2]

def isUniqueMapping(read):
    for tag in read.tags:
        if tag[0] == 'NH':#handle multiple mapping
            NHtag = int(tag[1])
            break
    if NHtag == 1:
        return True
    else:
        return False

ap = argparse.ArgumentParser()
ap.add_argument('--inbam_RNAi', help='Mapped reads in bam format, target bam for target RNAi  ')
ap.add_argument('--inbam_control', help='Mapped reads in bam format, control bam to get background dsRNA')
ap.add_argument('--inbam_ToCount', help='a list of the bam file you want to count for dsRNA')
ap.add_argument('--RNAiBC', help='A table file for BCs of RNAi samples')
ap.add_argument('--controlBC', help='BCs of all samples intended to be control when RNAi sample is excluded')
ap.add_argument('--ToCountBC', help='BCs of all samples intended to count in ToCorrect bam ')
ap.add_argument('--G2Ttable', help='the annotation table used in ESAT')
ap.add_argument('--outFile', help='Output file to save corrected counts')


args = ap.parse_args()
controlBCFile=args.controlBC
RNAiBCFile=args.RNAiBC
ToCountBCFile=args.ToCountBC
G2TtableFile=args.G2Ttable
inbam_ToCorrect=args.inbam_ToCount
outFile=args.outFile

# setup the bam file connection
print('Loading and indexing bam files...')
bam_RNAi = pysam.AlignmentFile(args.inbam_RNAi, "rb")
bam_ctr = pysam.AlignmentFile(args.inbam_control, "rb")
# iterate through the ToCount bam 
bam_ToCount_files = set(pd.read_csv(inbam_ToCorrect,header=None,sep= ',')[0])
bam_ToCount_list = {}
# bam_ToCount_name_indexed = {} # we giving up handling multiple mapping, now only look at unique mapping
for bamfile in bam_ToCount_files:
    bam_ToCount_list[bamfile] = pysam.AlignmentFile(bamfile, "rb")
    # to do optimal counting for multiple mapping reads, we build name index for the ToCount bam file
    # bam_ToCount_name_indexed[bamfile] = pysam.IndexedReads(bam_ToCount_list[bamfile])
    # bam_ToCount_name_indexed[bamfile].build()

# load the RNAi-barcode table
RNAiBCtbl = pd.read_csv(RNAiBCFile,index_col=0,header=None,sep= ',').T
controlBCs_all = set(pd.read_csv(controlBCFile,header=None,sep= ',')[0])
ToCountBCs = set(pd.read_csv(ToCountBCFile,header=None,sep= ',')[0])

# organize the RNAi targets by chromosome, and correct chr by chr
targetRNAiGenes = RNAiBCtbl.columns
targetRNAiGenes = set(targetRNAiGenes) - set(['vector'])
G2Ttable = pd.read_csv(G2TtableFile, sep='\t')
geneByChr = dict()
for gene in targetRNAiGenes:
    targetchr = G2Ttable.loc[G2Ttable['name2']==gene]['chrom'].iloc[0]
    if geneByChr.has_key(targetchr):
        geneByChr[targetchr].append(gene)
    else:
        geneByChr[targetchr] = [gene]


# each ToCount bam file will have a table written out
outputTbl = dict()
for bamfile in bam_ToCount_list:
    outputTbl[bamfile] = pd.DataFrame(columns = ToCountBCs)
dsRNAcount_ctr = dict()

# iterate each gene
for targetchr in geneByChr:
    #1. find genes that needs correction
    #       gene names, and gene position 
    # only look at the target chromosome 
    G2Ttable_sub = G2Ttable.loc[G2Ttable['chrom']==targetchr].reset_index(drop=True)
    # gene 2 table index dict
    ind2gene = G2Ttable_sub.to_dict()['name2']
    gene2ind = dict()
    for key, value in ind2gene.items(): 
       if value in gene2ind: 
           gene2ind[value].append(key) 
       else: 
           gene2ind[value]=[key] 
           
    # construct the mataRNA of each gene   
    gene2info = dict()
    gene2strd = dict()
    for gene,ind in gene2ind.items():
        exonStarts = G2Ttable_sub.iloc[ind,]['exonStarts'] #convert to 0-based for pysam
        exonEnds = G2Ttable_sub.iloc[ind,]['exonEnds']#convert to 0-based for pysam, so the length = end-start 
        startPos = []
        endPos = []
        for i in range(0,len(exonStarts)):
            isoStart = exonStarts.iloc[i,]
            isoEnd = exonEnds.iloc[i,]
            startPos+=[int(j) for j in isoStart.split(',')[:-1]]
            endPos+=[int(j) for j in isoEnd.split(',')[:-1]]
        startPos = [x - 1 for x in startPos]#convert to 0-based for pysam
        endPos = [x - 1 for x in endPos]#convert to 0-based for pysam
        IsExon = collections.OrderedDict(zip(range(min(startPos),max(endPos)+1), [0]*(max(endPos)+1-min(startPos)))) #use a pos:status dict
        for i in range(0,len(startPos)):
            for j in range(startPos[i],endPos[i]+1):
                IsExon[j]=1
        # we remove the introns to get the metaRNA
        metaRNA = []
        for pos,status in IsExon.items():
            if status:
                metaRNA.append(pos)
        gene2info[gene]=metaRNA # nested Dict
        gene2strd[gene]= G2Ttable_sub.iloc[ind[0],]['strand']
    
    # correct each gene 
    for RNAiWBID in geneByChr[targetchr]:
        print('Correcting for %s dsRNA...' % RNAiWBID)
        if args.inbam_RNAi == args.inbam_control: #if the same library is used for control (not recommend)
            controlBCs = controlBCs_all - set(RNAiBCtbl[RNAiWBID])
        else:
            controlBCs = controlBCs_all
  
         # map the density     
        metaRNA = set(gene2info[RNAiWBID])
                      
        # we first count the reverse strand of the target RNAi gene, in the control lirbary, to compare with the dsRNA count        
        dsRNA_ctr_count = 0
        if gene2strd[RNAiWBID]=='-':#gene is reverse strand   
            for read in bam_ctr.fetch(targetchr, min(metaRNA), max(metaRNA)):#we look at 5' end pile!
                if not read.is_reverse:
                    if read.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                        myPosition=read.reference_start - read.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                    else:
                        myPosition=read.reference_start
                    if myPosition in metaRNA:
                        if isUniqueMapping(read): #only look at unique map
                            if parseHeader(read.query_name) in controlBCs:
                                dsRNA_ctr_count+=1     
        else:
            for read in bam_ctr.fetch(targetchr, min(metaRNA), max(metaRNA)):#we look at 5' end pile!
                if read.is_reverse:
                    #5' end is the largest coordinate (original read is reverse complimented)
                    if read.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                        myPosition=read.reference_end + read.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                    else:
                        myPosition=read.reference_end
                    if myPosition in metaRNA:
                        if isUniqueMapping(read):
                            if parseHeader(read.query_name) in controlBCs:
                                dsRNA_ctr_count+=1   
        
        dsRNAcount_ctr[RNAiWBID] = dsRNA_ctr_count

        #count each condition 
        for bamfile in bam_ToCount_list:
            outputTbl[bamfile].loc[RNAiWBID] = [0] * len(ToCountBCs) #initial with all zero count
            for read in bam_ToCount_list[bamfile].fetch(targetchr,min(metaRNA), max(metaRNA)):
                if (gene2strd[RNAiWBID]=='+') == read.is_reverse:
                    # get 5' end position
                    if gene2strd[RNAiWBID]=='+':
                        if read.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                            myPosition=read.reference_end + read.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                        else:
                            myPosition=read.reference_end
                    else:
                        if read.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                            myPosition=read.reference_start - read.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                        else:
                            myPosition=read.reference_start
                     # count or not           
                    if myPosition in metaRNA: # we need to count this read to its barcode
                        myBC = parseHeader(read.query_name)
                        if myBC in ToCountBCs:
                            if isUniqueMapping(read):
                                outputTbl[bamfile].loc[RNAiWBID,myBC] = outputTbl[bamfile].loc[RNAiWBID,myBC] + 1
                                                                        
                            
# write output  
for bamfile in outputTbl:             
    outputTbl[bamfile].to_csv(outFile+'_dsRNA_count_matrix.csv',index=True,header=True)

dsRNAcount_ctr=pd.DataFrame.from_dict(dsRNAcount_ctr, orient='index')
dsRNAcount_ctr.to_csv(outFile+'_ctr_dsRNA_count_expanded.csv',index=True,header=False)

bam_RNAi.close()
bam_ctr.close()
for bamfile in bam_ToCount_files:
    bam_ToCount_list[bamfile].close()
