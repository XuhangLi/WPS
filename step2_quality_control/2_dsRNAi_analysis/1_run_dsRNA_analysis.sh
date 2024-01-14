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

outFileName="met7-lib2_sorted"
inputBamFile="./../../step1_process_raw_data/output/star/met7-lib2_sorted.bam"
RNAiBCsheet="./../../data/RNAiBCsheets_before_cleaning/RNAiBCsheet_met7_lib2.csv"
ctrBamFile="bam_ctr/Dev_SS3_trim_read_sorted.bam"
ctrBCsheet="bam_ctr/bcSet_full.txt"
BCtoCorrect="bcSet_full.txt"
G2Ttable="./../../data/metaData/WS279_G2Ttable.selected.txt"

echo ${inputBamFile} > BamToCorrect.txt 
awk -F, '$2!="X"' "$RNAiBCsheet" > "trimmedRNAiBCsheet.csv"


python correct_dsRNA.py --inbam_RNAi ${inputBamFile} --inbam_control ${ctrBamFile} --inbam_ToCorrect BamToCorrect.txt --RNAiBC trimmedRNAiBCsheet.csv --controlBC ${ctrBCsheet} --ToCorrectBC ${BCtoCorrect} --G2Ttable ${G2Ttable} --ExtLen 1000 --outFile ${outFileName} >> ${outFileName}_dsRNAcorrect.log

python count_dsRNA.py --inbam_RNAi ${inputBamFile} --inbam_control ${ctrBamFile} --inbam_ToCount BamToCorrect.txt --RNAiBC ${RNAiBCsheet} --controlBC ${ctrBCsheet} --ToCountBC ${BCtoCorrect} --G2Ttable ${G2Ttable} --outFile ${outFileName} >> ${outFileName}_dsRNAcount.log


rm BamToCorrect.txt
rm trimmedRNAiBCsheet.csv
