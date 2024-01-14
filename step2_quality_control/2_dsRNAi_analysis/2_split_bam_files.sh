#!/bin/bash

inputBamFile="./../../step1_process_raw_data/output/star/met7-lib2_sorted.bam" # bam file to be demultiplexed
bcTable="./../../data/RNAiBCsheets_before_cleaning/bcTable_met7_lib2.txt" # barcode table, generated in the previous step

python splitBam.py --inbam ${inputBamFile} \
--outbam met7 \
--barcodeTable ${bcTable}
