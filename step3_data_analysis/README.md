# WPS data analysis

## Introduction 
The last step of WPS data processing pipeline is running the WPS DE analysis that is specially design to leverage the power of massively parallel profiled conditions. This step also serve as a showcase on how to use the data produced and organized in the previous steps for downstream analysis. 

To learn more about **WPS DE framework**, see the **wpsDE** package [here](https://github.com/XuhangLi/wpsDE).

## Procedures
**wpsDE** is highly integrated with WPS data. Simply follow the [__1_WPS_DE_analysis walkthrough__](1_WPS_DE_analysis.Rmd).

## Cautions
Distinct from conventional DE analysis that analyzes on a per-condition-basis, **wpsDE** performs whole-dataset level statistical modeling to maximize the statistical power. This means you should combine all data collected for a WPS project (as long as they are experimented in the same setup). In the walkthrough, we showcased a setting that performed WPS of one RNAi plate (~96 conditions) in triplicates. If your experiment involves multiple plates, they should be combined before inputting to **wpsDE**. 

## Computational performance 
We acknowledge that **wpsDE** takes significant computing time on real data. For instance, the showcase dataset that includes ~100 conditions in triplicates (~300 samples) takes ~4 hours to complete on single core of a regular desktop computer. This computation time scales linearly with the number of libraries and number of genes. We have identified the bottleneck and are working on encoding parallelization functions to increase the performance, which we expect a linear decrease with more cores involved.  

