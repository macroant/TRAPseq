# A TRAPseq protocol and analysis pipeline to identify mutant-cell specific transcriptome
Targeted purification of polysomal mRNA (TRAP-Seq) is a recently developed transcriptome profiling technique that maps translating mRNAs under specific tagged-cells. Herein, we describe our TRAPseq protocol for transcriptome profiling. Our protocol outlines all necessary steps for transcriptome profiling including the bioinformatic pipeline steps required to process, analyze, and identify differentially expressed genes (DEGs). <BR>
see, https://www.google.com/ for more TRAPseq experiment details.

<img src="https://github.com/macroant/TRAPseq/blob/main/doc/overview.png" height="400" width="700">
## 

To replicate results presented in the manuscript above, please run the script after installation of all required libraries.


```bash
$ bash hello_work.sh para_1 para_2
# for example
$ bash hello_work.sh para_3 para_4
```

```bash
$ STAR --runThreadN threads_num --runMode genomeGenerate --genomeDir /path-to-index \
		 --genomeFastaFiles /path-to-fasta \
		 --sjdbGTFfile /path-to-GTF-file \
		 --sjdbOverhang 99
```

##Install required libraries
We are using "R version 3.4.1" here.
For R, 'DESeq2' library is required, which can be installed via bioconductor https://bioconductor.org/packages/release/bioc/html/DESeq2.html
```R
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("DESeq2")
```

## Workflow
The bioinformatics pipeline used for TRAPseq analysis consists of six main steps. Each of these steps performs a required function in processing the raw sequencing files (Fastq files) to obtain high-quality DEG list while the optional enrichment step performs the biological interpretation of DEGs once they are called. The six bioinformatics analysis steps are detailed below.<BR>
    
## <a name="workflow"></a> Step 1 Trimming & Alignment
Reads quality trimming and cleaned sequencing reads alignment to reference genomes. 
Firstly, [fastp](https://github.com/OpenGene/fastp) is run to remove low quality bases from raw sequencing reads while also removing any potential adapter sequences. <BR>
```bash
$ fastp -i R1.fastq.gz -I R2.fastq.gz \
    -o R1.fastq.gz.out.fq.gz -O R2.fastq.gz.out.fq.gz \
    -j test.json -h test.html -w threads_num
```
Once the reads are cleaned, alignment to the reference genomes is performed using [STAR](https://github.com/alexdobin/STAR). In this study, reads are aligned to the mouse genome by STAR 2.7.9. Following alignment, unaligned reads are removed, and alignment statistics are calculated.<BR>
```bash
$ STAR --genomeDir /path-to-STAR-index --readFilesIn R1.fastq.gz.out.fq.gz R2.fastq.gz.out.fq.gz \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
    --runThreadN threads_num --outFileNamePrefix output_test
```

## <a name="workflow"></a> Step 2 Raw Reads Quantification & Normalization
Assign aligned reads to genomic features (genes) and summarize raw reads count for each genomic feature and each sample. Normalize raw read counts to remove technical biases, e.g. sequencing depth and gene length, and make normalized gene expression values directly comparable within and across samples.
The gene abundance quantification tool [featureCount](https://subread.sourceforge.net/featureCounts.html) is used to calculate raw read count for each gene and each sample. We used the gene expression detection threshold where genes were selected if there were ≥6 reads in at least 20% of samples. <BR>
The relative abundance of a transcript is calculated by Fragments per Kilobase per Million Mapped Fragments (FPKM) and Transcripts per million (TPM).<BR>

## <a name="workflow"></a> Step 3 Detection of Outlier Samples
Evaluate and remove outlier samples can significantly improve the reliability of DEGs and downstream functional analysis. We used [robustPCA](https://cran.r-project.org/web/packages/rrcov/index.html) for outlier detection as it was found to result in stable performance and  the default cutoff value allows correct identification of outliers without yielding false positives in an evaluation study.<BR>

## <a name="workflow"></a> Step 4 Cell Type Deconvolution
TRAPseq enriches specific cells than other cells, cell type deconvolution analysis can help reveal the cell type composition of TRAPseq transcriptome profiles. We estimated the cellular compositions of the TRAPseq samples through deconvolution analysis. We obtained gene expression signatures for major tubules segments from a mouse micro-dissected tubules dataset (GSE150338, PMID: 33769951), and used as input into [CIBERSORTx](https://cibersortx.stanford.edu) to estimate the cellular composition of our TRAPseq data. <BR>

## <a name="workflow"></a> Step 5 Detection of Differential Expressed Genes & Filtering
This step detects genes with statistically significant differences or changes in read counts or expression levels between two experimental conditions. The statistical significance test is conducted using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) with negative binomial distribution and Benjamini-Hochberg adjusted p-values (FDR). To screen for more meaningful DEGs, we only consider genes that are protein-coding genes and have the same expression change direction in comparison of SKO/WT and SKO/DKO. <BR>

## <a name="workflow"></a> Step 6 Enrichment Analysis (optional analysis step)
Biological interpretations of selected DEGs are performed using various bioinformatics tools, including [DAVID](https://david.ncifcrf.gov), [Metascape](https://metascape.org/gp/index.html#/main/step1), and commercal tool [IPA](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/) or [MetaCore](https://clarivate.com/products/biopharma/discovery-clinical-regulatory/early-research-intelligence-solutions/).

## Useful References
