#!/bin/bash

# check and install dependencies 
check_and_install() {
    package=$1
    if python -c "import $package" > /dev/null 2>&1; then
        echo "$package is already installed."
    else
        echo "$package is not installed. Installing..."
        pip install $package
    fi
}

# Check and install pandas
check_and_install pandas

# Check and install pysam
check_and_install pysam

# Check and install numpy
check_and_install numpy

# set parameters 

outFileName="met7-lib2_sorted" # prefix of the output files. For standard WPS, please use the input bam file name
inputBamFile="./../../step1_process_raw_data/output/star/met7-lib2_sorted.bam" # path to the input bam file for the RNAi library to be analyzed
RNAiBCsheet="./../../data/RNAiBCsheets_before_cleaning/RNAiBCsheet_met7_lib2.csv" # corresponding RNAiBCsheet produced in the previous step.
ctrBamFile="bam_ctr/Dev_SS3_trim_read_sorted.bam" # no-RNAi control bam file. If no custom file available, you can use the same control file here.
ctrBCsheet="bam_ctr/bcSet_full.txt" # all barcodes in the control bam file
BCtoCorrect="bcSet_full.txt" # barcodes in the input bam file that we want to correct on. We recommend including all possible barcodes. 
G2Ttable="./../../data/metaData/WS279_G2Ttable.selected.txt" # the Gene-to-transcript table generated in the previous raw data processing step

echo ${inputBamFile} > BamToCorrect.txt 
awk -F, '$2!="X"' "$RNAiBCsheet" > "trimmedRNAiBCsheet.csv"


python correct_dsRNA.py --inbam_RNAi ${inputBamFile} --inbam_control ${ctrBamFile} --inbam_ToCorrect BamToCorrect.txt --RNAiBC trimmedRNAiBCsheet.csv --controlBC ${ctrBCsheet} --ToCorrectBC ${BCtoCorrect} --G2Ttable ${G2Ttable} --ExtLen 1000 --outFile ${outFileName} >> ${outFileName}_dsRNAcorrect.log

python count_dsRNA.py --inbam_RNAi ${inputBamFile} --inbam_control ${ctrBamFile} --inbam_ToCount BamToCorrect.txt --RNAiBC ${RNAiBCsheet} --controlBC ${ctrBCsheet} --ToCountBC ${BCtoCorrect} --G2Ttable ${G2Ttable} --outFile ${outFileName} >> ${outFileName}_dsRNAcount.log


rm BamToCorrect.txt
rm trimmedRNAiBCsheet.csv
