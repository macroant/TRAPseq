# A TRAPseq protocol and analysis pipeline to identify mutant-cell specific transcriptome
Targeted purification of polysomal mRNA (TRAP-Seq) is a recently developed transcriptome profiling technique that maps translating mRNAs under specific tagged-cells. Herein, we describe our TRAPseq protocol for transcriptome profiling. Our protocol outlines all necessary steps for transcriptome profiling including the bioinformatic pipeline steps required to process, analyze, and identify differentially expressed genes (DEGs). <BR>
see, https://www.google.com/ for more TRAPseq experiment details.

<img src="https://" height="400" width="500">
## 

To replicate results presented in the manuscript above, please run the script after installation of all required libraries.


```bash
$ bash hello_work.sh para_1 para_2
# for example
$ bash hello_work.sh para_3 para_4
```

##Install required libraries
We are using "R version 3.4.1" here.
For R, 'DESeq2' library is required, which can be installed via bioconductor https://bioconductor.org/packages/release/bioc/html/DESeq2.html
```R
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("DESeq2")
```
## Useful References
