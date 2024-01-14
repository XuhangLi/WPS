import pysam
import argparse


#this code is the optimized version of collapseDup_fast, where memory usage is minimzed 
#the v2 version:
#counting is done by aligning the 5' end of each reads, instead of leftmost coordinate in the aligned chromosome ( which is reference_start tag)
#also the PCR duplicates are written out for quality controls


def parseHeader(header):
    # extract the sequences that should be the BC:
    tmp1 = header.split(" ")
    tmp2 = tmp1[0].split(':')
    #use hard index to improve performance
    return tmp2[-2] #+ tmp2[-1] #only extract the Barcode

def parseBarcode(text,bcDict):
    # extract the sequences that should be the BC:
    text = text.replace("\n","")
    text = text.replace("\r","")
    tmp1 = text.split(":")
    tmp2 = tmp1[1].split(',')
    for key in tmp2:
        bcDict[key] = tmp1[0]
    return bcDict #+ tmp2[-1] #only extract the Barcode


ap = argparse.ArgumentParser()
ap.add_argument('--inbam', help='Mapped reads in bam format')
ap.add_argument('--outbam', help='Output file suffix to save reads')
ap.add_argument('--barcodeTable', help='barcode table for spliting bam file')

#ap.add_argument("--u",help="a binary flag used to indicate that only uniquely mapped reads will be considered. By default uniquely mapped reads are defined as reads with MAPQ=60",action="store_true")
#ap.add_argument('--end',action='store_true',help='End position of read will not be considered (do NOT use False here! function not supported!')


args = ap.parse_args()
bam=args.inbam
out=args.outbam
barcodeT = args.barcodeTable


barcodeTable = open(barcodeT,"r")
#parse barcode table
bcDict = dict()
for lines in barcodeTable:
    bcDict=parseBarcode(lines,bcDict)


bam_header = pysam.Samfile(bam, 'rb').header
#create file handle for each barcode
for name in set(bcDict.values()):
    exec(name + '=pysam.Samfile("' + out + '_' + name +'.bam", "wb", header=bam_header)')


samfile = pysam.Samfile(bam, "rb" )
skippedReads = 0;
processedReads = 0;
for read in samfile.fetch():
    #write to corresponding file
    try:
        exec(bcDict[parseHeader(read.query_name)]+'.write(read)')
    except:
        skippedReads+=1
    if processedReads%1000000==1:
        print ('processed '+str(int(processedReads/1000000))+' million reads;')
    processedReads+=1
    

for name in set(bcDict.values()):
    exec(name + '.close()')
print ("DONE!")
print ("processed "+str(processedReads)+" reads!")
print ("skiped "+str(skippedReads)+" reads!")


