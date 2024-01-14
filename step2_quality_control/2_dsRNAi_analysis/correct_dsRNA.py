import pysam
import argparse
import pandas as pd
import collections
import numpy as np

# constant parameters
bin_size_ctr = 100 #how many bases are tolerated in interval calling for control calculation
bin_size_dsRNA = 100 #how many bases are tolerated in interval calling 

minIntensity = 100# we assume detection limit as 100 reads per gene (assume 2*50)
reads_kept_perc_signal = 0.5 #we require the main bins to keep over 50% reads for target gene to count. Noise level is assumed to be less than 10%
reads_kept_perc_dsRNA = 0.9 #we require to filter out 90% of dsRNA contamination
dsRNA_absCutoff_per_bin = 5 #a contaminated bin must have more than 5 dsRNA reads (because we found if the reverse str is less than 5, the mRNA str usually has no contamination)

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
ap.add_argument('--inbam_RNAi', help='Mapped reads in bam format, target bam for correction ')
ap.add_argument('--inbam_control', help='Mapped reads in bam format, control bam to get reference intervals')
ap.add_argument('--inbam_ToCorrect', help='a list of the bam file you want to recount with the clean regions')
ap.add_argument('--RNAiBC', help='A table file for BCs of RNAi samples')
ap.add_argument('--controlBC', help='BCs of all samples intended to be control when RNAi sample is excluded')
ap.add_argument('--ToCorrectBC', help='BCs of all samples intended to correct in ToCorrect bam')
ap.add_argument('--G2Ttable', help='the annotation table used in ESAT')
ap.add_argument('--ExtLen', help='the end extension length used in ESAT')
ap.add_argument('--outFile', help='Output file to save corrected counts')

# we will correct for the target gene as well as the flanking genes whose 3' extensive collapsed into target gene's gene body 
# once a countable region is redefined, this gene in all barcodes in the ToCount bams will be recounted
#1. find genes that needs correction
#       gene names, and gene position (same meta transcript as in ESAT)
#2. iterate through every gene
#       1. find the control reads positions
#       2. calculate counts for each bin
#       3. remove the low density bins 
#       4. find and remove the RNAi contaminated bins by the same procedure 
#       6. count the reads for all barcodes in all ToCount files

# In V2 version, we use 5' end position instead of density to estimate the countable region (the RT2 start site)
# also, we use regular bins instead of adaptive interval to handle control and dsRNA region calculation consistently and simontanously 
# the  3' extension is set in the same way as ESAT
# we first added but then removed the multiple mapping "proper" method as ESAT; because haveing it causes too much ~200G mem usage and slow speed(we need to account for multiple mapping in both countable region, contamination region calculation and the couting step!)

# in V3 version, we further optimized the running speed
# in V3.2, we exclude introns in the density analysis, since they may influence the proper noise cutoff estimation!(i.e WBGene00010275)

# Remaining problem: if the two RNAi target genes are next to each other in the genome and they are in the same RNAi library (bam), their dsRNA will cross-contaminate each other. 
# we didnt check for this because it is rare. When the output is used, need to check if one gene is corrected twice.

args = ap.parse_args()
controlBCFile=args.controlBC
RNAiBCFile=args.RNAiBC
ToCorrectBCFile=args.ToCorrectBC
G2TtableFile=args.G2Ttable
ExtLen_ori=int(args.ExtLen)
inbam_ToCorrect=args.inbam_ToCorrect
outFile=args.outFile

# setup the bam file connection
print('Loading and indexing bam files, this may take long time...')
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
ToCorrectBCs = set(pd.read_csv(ToCorrectBCFile,header=None,sep= ',')[0])

# organize the RNAi targets by chromosome, and correct chr by chr
targetRNAiGenes = RNAiBCtbl.columns
targetRNAiGenes = set(targetRNAiGenes) - set(['vector'])
G2Ttable = pd.read_csv(G2TtableFile,sep= '\t')
geneByChr = dict()
for gene in targetRNAiGenes:
    targetchr = G2Ttable.loc[G2Ttable['name2']==gene]['chrom'].iloc[0]
    if geneByChr.has_key(targetchr):
        geneByChr[targetchr].append(gene)
    else:
        geneByChr[targetchr] = [gene]

# initialize outputs
dsRNAcount = dict()
dsRNAcount_ctr = dict()
# each ToCount bam file will have a table written out
outputTbl = dict()
for bamfile in bam_ToCount_list:
    outputTbl[bamfile] = pd.DataFrame(columns = ToCorrectBCs)

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
    
    # calculate the exact extension length of each gene
    loc2status = dict()# not key ==> no gene; '+' ,'-','b'
    for gene, metaRNA in gene2info.items(): 
        strd = gene2strd[gene]
        for loc in metaRNA:
            if loc in loc2status: 
                loc2status[loc] = 'b'
            else:
                loc2status[loc] = strd
    # get extension length
    gene2extL = dict()
    for gene in gene2info:
        if gene2strd[gene] == '+':
            ext_loc = range(max(gene2info[gene])+1,max(gene2info[gene])+ExtLen_ori+1)
        else:
            ext_loc = range(min(gene2info[gene])-1,min(gene2info[gene])-ExtLen_ori-1,-1)
        for loc in ext_loc:
            if loc2status.has_key(loc):
                if loc2status[loc] == gene2strd[gene] or loc2status[loc] == 'b':
                    gene2extL[gene] = ext_loc.index(loc)
                    break
        if not gene2extL.has_key(gene):
            gene2extL[gene] = 1000
    
    # correct each gene 
    for RNAiWBID in geneByChr[targetchr]:
        print('Correcting for %s dsRNA...' % RNAiWBID)
        if args.inbam_RNAi == args.inbam_control: #if the same library is used for control (not recommend)
            controlBCs = controlBCs_all - set(RNAiBCtbl[RNAiWBID])
        else:
            controlBCs = controlBCs_all
        
        # find the genes to be corrected
        gene2correct = [RNAiWBID]
        checkRegion = set(gene2info[RNAiWBID])
        for gene,metaRNA in gene2info.items():
            for loc in metaRNA:
                if loc in checkRegion:
                    gene2correct.append(gene)
                    break
            # apply the 3' extension
            if gene2strd[gene] == '+':
                ext_loc = range(max(metaRNA)+1,max(metaRNA)+gene2extL[gene]+1)
            else:
                ext_loc = range(min(metaRNA)-gene2extL[gene],min(metaRNA))
            for loc in ext_loc:
                if loc in checkRegion:
                    gene2correct.append(gene)
                    break
        gene2correct = set(gene2correct)
        
        #2. iterate through every gene
        #       1. find the control reads positions
        #       2. calculate counts for each bin
        #       3. remove the low density bins 
        #       4. remove the RNAi contaminated bins
        #       6. count the reads for RNAi barcode 
        #       * before doing this, finding the dsRNA contaminated region by the same procedure
        
        # we first define the dsRNA contsminated region 
        # RNAi cotaminated region is defined by the revserse strand mapping; and is defined as a large, continous area, due to uncertainty of RT initial site
         # inital a density dict variable(dsRNA has no extension)
        dsRNAdensity = collections.OrderedDict() #the complementary strand!
        for key in gene2info[RNAiWBID]:
            dsRNAdensity[key] = 0
         # map the density     
        metaRNA = set(gene2info[RNAiWBID])
        if gene2strd[RNAiWBID]=='-':#gene is reverse strand   
            for read in bam_RNAi.fetch(targetchr, min(metaRNA), max(metaRNA)):#we look at 5' end pile!
                if not read.is_reverse:
                    if read.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                        myPosition=read.reference_start - read.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                    else:
                        myPosition=read.reference_start
                    if myPosition in metaRNA:
                        if isUniqueMapping(read): #only look at unique map
                            if parseHeader(read.query_name) in set(RNAiBCtbl[RNAiWBID]):
                                dsRNAdensity[myPosition]=dsRNAdensity[myPosition]+1     
        else:
            for read in bam_RNAi.fetch(targetchr, min(metaRNA), max(metaRNA)):#we look at 5' end pile!
                if read.is_reverse:
                    #5' end is the largest coordinate (original read is reverse complimented)
                    if read.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                        myPosition=read.reference_end + read.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                    else:
                        myPosition=read.reference_end
                    if myPosition in metaRNA:
                        if isUniqueMapping(read):
                            if parseHeader(read.query_name) in set(RNAiBCtbl[RNAiWBID]):
                                dsRNAdensity[myPosition]=dsRNAdensity[myPosition]+1                         
        #  make new bins only covering metaRNA
        binStart = range(0,len(dsRNAdensity.keys()),bin_size_dsRNA); 
        binEnd = range(bin_size_dsRNA-1,len(dsRNAdensity.keys()),bin_size_dsRNA);
        if len(binEnd) < len(binStart):
            binEnd.append(len(dsRNAdensity.keys())-1)
        # calculate counts for each bin
        intensityCount = [0] * len(binStart)
        positions = dsRNAdensity.keys()
        for i in range(0,len(binStart)):
            for ind in range(binStart[i],binEnd[i]+1):
                intensityCount[i]+=dsRNAdensity[positions[ind]]
        # remove the low density bins 
        if sum(intensityCount) != 0:
            intensityNorm = [float(x) / float(sum(intensityCount)) for x in intensityCount]
        else:
            intensityNorm = [0] * len(intensityCount)
            lowDensityCutoff = 1
            print('Notice: %s RNAi did not find dsRNA contamination!:'% RNAiWBID)
            
        dsRNAcount[RNAiWBID] = sum(intensityCount)
        # find the proper cutoff
        for cutoff in np.arange(1,0,-0.005):
            if sum([x for x in intensityNorm if x > cutoff]) > reads_kept_perc_dsRNA:
                lowDensityCutoff = cutoff
                break  
        # cutoff low density and label the un-countable region
        # initial a count indicator variable 
        NotToCount = collections.OrderedDict()
        for key in gene2info[RNAiWBID]:
            NotToCount[key] = 0
        # label regions
        for i in range(0,len(binStart)):
            if intensityNorm[i] > lowDensityCutoff and intensityCount[i] > dsRNA_absCutoff_per_bin:
                for ind in range(binStart[i],binEnd[i]+1):
                    NotToCount[positions[ind]] = 1
        
        
        # we also additionally count the reverse strand of the target RNAi gene, in the control lirbary, to compare with the dsRNA count        
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


        # then, we loop through every gene to correct
        print('Recounting...')
        for gene in gene2correct:
            # 1. find the control reads positions
            # go through the genome region of this gene, and label the density of the reads
            metaRNA = set(gene2info[gene])
            # inital a density dict variable
            if gene2strd[gene]=='+':
                refDensity = collections.OrderedDict()
                for key in gene2info[gene]:
                    refDensity[key] = 0
                PileStart = min(metaRNA)
                PileEnd = max(metaRNA)+gene2extL[gene] #set extension 
                for loc in range(max(metaRNA)+1,max(metaRNA)+gene2extL[gene]+1):#set extension 
                    refDensity[loc] = 0
            else:
                refDensity = collections.OrderedDict()
                PileStart = min(metaRNA)-gene2extL[gene]
                PileEnd = max(metaRNA)
                for loc in range(min(metaRNA)-gene2extL[gene],min(metaRNA)):#set extension 
                    refDensity[loc] = 0
                for key in gene2info[gene]:
                    refDensity[key] = 0
                    
            if gene2strd[gene]=='-':#reverse strand   
                for read in bam_ctr.fetch(targetchr, PileStart, PileEnd):#we look at 5' end pile!
                    if read.is_reverse:
                        #5' end is the largest coordinate (original read is reverse complimented)
                        if read.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                            myPosition=read.reference_end + read.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                        else:
                            myPosition=read.reference_end
                        if refDensity.has_key(myPosition):
                            if isUniqueMapping(read):
                                if parseHeader(read.query_name) in controlBCs:
                                    refDensity[myPosition]=refDensity[myPosition]+1
            else:#forward strand   
                for read in bam_ctr.fetch(targetchr, PileStart, PileEnd):#we look at 5' end pile!
                    if not read.is_reverse:
                        #5' end is the largest coordinate (original read is reverse complimented)
                        if read.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                            myPosition=read.reference_start - read.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                        else:
                            myPosition=read.reference_start
                        if refDensity.has_key(myPosition):
                            if isUniqueMapping(read):
                                if parseHeader(read.query_name) in controlBCs:
                                    refDensity[myPosition]=refDensity[myPosition]+1
                               
            # 2. calculate counts for each bin
            # first make bins
            binStart = range(0,len(refDensity.keys()),bin_size_ctr); 
            binEnd = range(bin_size_ctr-1,len(refDensity.keys()),bin_size_ctr);
            if len(binEnd) < len(binStart):
                binEnd.append(len(refDensity.keys())-1)
            # count each bin's total reads
            intensityCount = [0] * len(binStart)
            positions = refDensity.keys()
            for i in range(0,len(binStart)):
                for ind in range(binStart[i],binEnd[i]+1):
                    intensityCount[i]+=refDensity[positions[ind]]
            # 3. remove the low density bins 
            if sum(intensityCount) != 0:
                intensityNorm = [float(x) / float(sum(intensityCount)) for x in intensityCount]
            else:
                intensityNorm = [0] * len(intensityCount)
            # find the proper cutoff
            for cutoff in np.arange(1,0,-0.005):
                if sum([x for x in intensityNorm if x > cutoff]) > reads_kept_perc_signal:
                    lowDensityCutoff = cutoff
                    break  
            # QC: check if the reads are skewed towards 3'end
            if sum(intensityCount) > minIntensity: # we assume delection limit as 100 reads per gene
                maxPos = 0.5*(binStart[intensityNorm.index(max(intensityNorm))]+binEnd[intensityNorm.index(max(intensityNorm))])
                if (gene2strd[gene] == '+' and (maxPos < 0.5*(len(metaRNA)))) or \
                     (gene2strd[gene] == '-' and (maxPos > 0.5*(len(metaRNA))+gene2extL[gene])):
                         print('Warning: in RNAi sample %s, %s skews towards 5-end!'% (RNAiWBID,gene))
            else:
                print('Warning: in RNAi sample %s, %s is below detection limit, the correction may be invalid!'% (RNAiWBID,gene))
            # remove low density bins and label the countable region
            # initial a count indicator variable 
            if gene2strd[gene]=='+':
                ToCount = collections.OrderedDict()
                for key in gene2info[gene]:
                    ToCount[key] = 0
                for loc in range(max(metaRNA)+1,max(metaRNA)+gene2extL[gene]+1):#set extension 
                    ToCount[loc] = 0
            else:
                ToCount = collections.OrderedDict()
                for loc in range(min(metaRNA)-gene2extL[gene],min(metaRNA)):#set extension 
                    ToCount[loc] = 0
                for key in gene2info[gene]:
                    ToCount[key] = 0
                    
            for i in range(0,len(binStart)):
                if intensityNorm[i] > lowDensityCutoff:
                    for ind in range(binStart[i],binEnd[i]+1):
                            ToCount[positions[ind]] = 1
                   
            # 4. remove the RNAi contaminated region               
            #  update the ToCount dict
            for pos in ToCount:
                if NotToCount.has_key(pos):
                    if NotToCount[pos] == 1:
                        ToCount[pos] = 0
            # a quick QC
            countablePos = 0
            for pos in ToCount:
                countablePos+=ToCount[pos]
            if countablePos == 0 and sum(intensityCount) > 0:
                print('Warning: in RNAi sample %s, %s has no countable region! Except zero counts!'% (RNAiWBID,gene))
            else:
                print('info: the countable region of %s is %d.'% (gene ,countablePos))
    
            #  6. count the reads for all barcodes in all ToCount files
            # bam_ToCount_list = {} bam_ToCount_name_indexed = {}
            for bamfile in bam_ToCount_list:
                if gene not in outputTbl[bamfile].index:
                    outputTbl[bamfile].loc[gene] = [0] * len(ToCorrectBCs) #initial with all zero count
                else:
                    raise Exception("%s gene is corrected twice! Not allowed!" % gene)
                for read in bam_ToCount_list[bamfile].fetch(targetchr,start =  PileStart, stop = PileEnd):
                    if (gene2strd[gene]=='-') == read.is_reverse:
                        # get 5' end position
                        if gene2strd[gene]=='-':
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
                        if ToCount.has_key(myPosition):
                            if ToCount[myPosition] == 1: # we need to count this read to its barcode
                                myBC = parseHeader(read.query_name)
                                if myBC in ToCorrectBCs:
                                    #for tag in read.tags:
                                    #   if tag[0] == 'NH':#handle multiple mapping
                                    #       NHtag = int(tag[1])          
                                    #if NHtag == 1:
                                    if isUniqueMapping(read):
                                        outputTbl[bamfile].loc[gene,myBC] = outputTbl[bamfile].loc[gene,myBC] + 1
                                    #else:# check if optimal multiple mapping applicable
                                    #    N_hit_gene = 0
                                    #   for read2 in bam_ToCount_name_indexed[bamfile].find(read.query_name):
                                    #       # again, get 5' end position
                                    #      if read2.is_reverse:
                                    #           if read2.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                                    #               myPosition2=read2.reference_end + read2.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                                    #           else:
                                    #               myPosition2=read2.reference_end
                                    #       else:
                                    #           if read2.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                                    #               myPosition2=read2.reference_start - read2.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                                    #           else:
                                    #               myPosition2=read2.reference_start
                                    #       # check if read2 is in a gene exon
                                    #       N_hit_gene+= loc2status.has_key(myPosition2)
                                    #       if N_hit_gene >= 2:
                                    #           break
                                    #   if N_hit_gene == 1: # qualified
                                    #       outputTbl[bamfile].loc[gene,myBC] = outputTbl[bamfile].loc[gene,myBC] + 1
                                                
                            
# write output  
if len(outputTbl) == 1:
    bamfile, table = next(iter(outputTbl.items()))
    table.to_csv(outFile+'_recount_matrix.csv',index=True,header=True)
else:
    for bamfile in outputTbl:
	print('Multiple Bam file is being corrected. Saving the recount table in the bam file location...')             
        outputTbl[bamfile].to_csv(bamfile+'_recount_matrix.csv',index=True,header=True)

dsRNAcount=pd.DataFrame.from_dict(dsRNAcount, orient='index')
dsRNAcount.to_csv(outFile+'_target_dsRNA_count.csv',index=True,header=False)

dsRNAcount_ctr=pd.DataFrame.from_dict(dsRNAcount_ctr, orient='index')
dsRNAcount_ctr.to_csv(outFile+'_ctr_dsRNA_count.csv',index=True,header=False)

bam_RNAi.close()
bam_ctr.close()
for bamfile in bam_ToCount_files:
    bam_ToCount_list[bamfile].close()
