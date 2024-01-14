#!/bin/bash

inputBamFile="./../../step1_process_raw_data/output/star/met7-lib2_sorted.bam"
bcTable="./../../data/RNAiBCsheets_before_cleaning/bcTable_met7_lib2.txt"

python splitBam.py --inbam ${inputBamFile} \
--outbam met7 \
--barcodeTable ${bcTable}
