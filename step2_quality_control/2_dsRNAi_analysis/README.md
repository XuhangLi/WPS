# Run dsRNA analysis
This is a walkthrough to execute dsRNA analysis using an in-house Python script. The inputs needed for this analysis are generated in the previous steps of this pipeline. Please note, if you change the location of any files produced by the previous steps, you may encounter errors or need to modify the input paths manually.
### programs
There are three python scripts for this analysis, including:
* **correct_dsRNA.py**: this script corrects the dsRNA contaminations on gene expression quantification. It generates a recounted read count table for contaminated genes in a given library. 
* **count_dsRNA.py**: this script quantifies the dsRNA read counts for both target genes in a libary and other user-supplied genes that can be potential swapped clones.
* **splitBam.py**: this script is a helper function for the convenience of the users. It demultiplexed a bam file for the entire library into individual bam files for each RNAi condition. This is useful when inspecting individual condition and its dsRNA signals is desired. 

### usage
We use met7-lib2 as an example to illustrate how to use the scripts. 

For a seamless analysis, we provided two wrapper bash scripts to run the python programs easily, including:
* _1_run_dsRNA_analysis.sh_: this script will execute correct_dsRNA.py and count_dsRNA.py with inputs defined in the script. The outputs will be save in the same folder as _1_run_dsRNA_analysis.sh_.
* _2_split_bam_files.sh_: this script will execute splitBam.py with inputs defined in the script. The outputs will be save in the same folder as _2_split_bam_files.sh_. Splitting bam files is optional for WPS pipeline.

__NOTICE__: To perform a test run of met7-lib2 dsRNA analysis, you need to download the bam files for both met7-lib2 and a control library that does not target to any gene (but from animals fed on vector control RNAi). These files were excluded from GitHub deposition due to their large size (see [gitignore](./../../.gitignore)). Please download from [here](wormflux link) and save the files in corresponding directory (see instructions in the downloaded files).

To run analysis:

```
sh 1_run_dsRNA_analysis.sh
```
and
```
sh 2_split_bam_files.sh
```
The parameters are explained by the comments in these bash scripts. 

### output
**1_run_dsRNA_analysis.sh** will generate five csv files as follows:
* met7-lib2_sorted_recount_matrix.csv: a matrix for recounted read counts for all dsRNA-contaminated genes in the input library.
* met7-lib2_sorted_target_dsRNA_count.csv: the dsRNA read counts of RNAi targeted genes in corresponding conditions (three replicates added up) in the input library.
* met7-lib2_sorted_ctr_dsRNA_count.csv: the dsRNA read counts of RNAi targeted genes in the no-RNAi control library. This is the  library-wise total count. 
* met7-lib2_sorted_dsRNA_count_matrix.csv: a matrix for the dsRNA read counts for all genes specified in the input RNAiBCtable.
* met7-lib2_sorted_ctr_dsRNA_count_expanded.csv: the dsRNA read counts of all specified genes in the input RNAiBCtable, but using library-wise total count in the no-RNAi control library.

In real data analysis, one should move these five files to "/outputs/RNAiAnalysis_before_cleaning/" (or "/outputs/RNAiAnalysis_after_cleaning/" if it is a rerun after updating the bad clones) such that they will be identified automatically in the next-step RNAi quality control. 
