# GeneLab bioinformatics processing protocol for 454 and IonTorrent amplicon sequencing data - in review

> **The document [`GL-DPPD-XXXX.md`](GL-DPPD-XXXX.md) holds an overview and example code of how GeneLab processes 454 and IonTorrent amplicon datasets. Exact processing code for specific datasets that have been released is available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and is also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

**Date:** November 28, 2020  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee (GeneLab Science Team)

**Approved by:**  

# Details for reviewing
* The primary GeneLab protocol document for review, holding main steps/programs and example code, is ['GL-DPPD-XXXX.md'](GL-DPPD-XXXX.md).

* [GLDS-72](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/) was done as an example. 
  * The ['example-output' directory](example-output) in this repo holds everything except the read files due to size. 
  * The processing is implemented as a Snakemake workflow. The primary snakefile is [here](example-output/processing_info/Snakefile).
  * The `README.txt` in the ['example-output' directory](example-output) describes the output contents and can be found [here](example-output/README.txt), and is also pasted below.
  
<br>
Thanks! 🙂

---

Contents of ['example-output' directory](example-output):

```
##################################################################################
## This directory holds processed data for NASA GLDS-72                         ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/                 ##
##                                                                              ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
## Based on GL-DPPD-XXXX for 454 or IonTorrent amplicon data                    ##
## Annotated code is in "processing_info.tar"                                   ##
##################################################################################

Summary of contents:

    - README.txt                        - this file

    - FastQC_Outputs/                   - multiQC summary reports of FastQC runs

    - Trimmed_Sequence_Data/            - primer-trimmed sequences and trimmed-read count table
        - *.fastq.gz                    - primer-trimmed sequence fastq files
        - trimmed-read-counts.tsv       - starting and trimmed-read counts per sample
        - cutadapt.log                  - log file from cutadapt

    - Filtered_Sequence_Data/           - quality-filtered sequence fastq files and filtered-read count table
        - *.fastq.gz                    - quality-filtered reads
        - filtered-read-counts.tsv      - starting and filtered-read counts per sample
        - bbduk.log                     - log file from bbduk

    - Final_Outputs/                    - primary output files
        - OTUs.fasta                    - fasta file of generated OTUs
        - counts.tsv                    - count table of OTUs across samples
        - taxonomy.tsv                  - taxonomy of OTUs
        - taxonomy-and-counts.biom      - biom-formatted output of taxonomy and counts
        - taxonomy-and-counts.tsv       - combined table of counts and taxonomy across samples
        - read-count-tracking.tsv       - read counts at each processing step

    - processing_info.tar               - a compressed directory holding info related to processing
        - unique-sample-IDs.txt         - single-column file of unique sample identifiers
        - Snakefile                     - snakemake workflow file
        - full-R-processing.R           - all R processing code used
        - primers.fa                    - primer sequences
        - snakemake-run.log             - log file from workflow execution
        - R-processing.log              - log file from R processing
        - vsearch.log                   - log file from vsearch
```
