# Bioinformatics pipeline for 454 and IonTorrent amplicon sequencing data  

> **This page holds an overview and some example code of how GeneLab processes 454 and IonTorrent amplicon datasets. Exact processing code for specific datasets that have been released is available in this repository [here](GLDS_Processing_Scripts) and is also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** October 30, 2020  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee (GeneLab Analysis Team)

**Approved by:**  

---

# Table of contents  

- [Software used](#software-used)
- [Reference databases used](#reference-databases-used)
- [General processing overview with example code](#general-processing-overview-with-example-code)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [Compile Raw Data QC](#compile-raw-data-qc)
  - [2. Trim Primers](#2-trim-primers)
  - [3. Quality filtering](#3-quality-filtering)
  - [4. Filtered Data QC](#4-filtered-data-qc)
    - [Compile Filtered Data QC](#compile-filtered-data-qc)
  - [5. Generating OTUs and counts per sample](#5-generating-otus-and-counts-per-sample)
    - [Dereplicate individual samples](#dereplicate-individual-samples)
    - [Generate OTUs](#generate-otus)
  - [6. Generating taxonomy and additional outputs](#6-generating-taxonomy-and-additional-outputs)
    - [Assigning taxonomy](#assigning-taxonomy)
    - [Generating and writing outputs](#generating-and-writing-outputs)

---

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|-------------:|
|FastQC|`fastqc -v`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`multiqc -v`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`cutadapt --version`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|bbduk|`bbduk.sh --version`|[https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)|
|vsearch|`vsearch --version`|[https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)|
|R|`R --version` (at command line) | [https://www.r-project.org/](https://www.r-project.org/)|
|DECIPHER|`packageVersion("DECIPHER")`|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|`packageVersion("biomformat")`|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|

>**\*** Exact versions are available along with the [processing code](glds_processing_scripts) for each specific dataset.

# Reference databases used

|Program used| Database| Relevant Links|
|:-----|:-----:|--------:|
|DECIPHER| SILVA SSU r138 | [http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData](http://www2.decipher.codes/Classification/TrainingSets/)|
|DECIPHER| UNITE v2020 | [http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2020_February2020.RData](http://www2.decipher.codes/Classification/TrainingSets/)|

---

# General processing overview with example code

> Exact processing code for specific datasets is available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory of this repository, as well as being provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

---

## 1. Raw Data QC

```
fastqc -o raw_fastqc_output *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input files:**

* fastq, compressed or uncompressed

**Output files:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)


<br>  

### Compile Raw Data QC

```
multiqc -o raw_multiqc_output raw_fastqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input file types:**

* fastqc.zip (FastQC output data)

**Output file types:**

* multiqc_report.html (multiqc output html summary)
* multiqc_data (directory containing multiqc output data)

<br>  

---

## 2. Trim Primers

The location and orientation of primers in the data is important to understand in deciding how to do this step. `cutadapt` has many options for primer identification and removal. They are described in detail on their documentation page here: [https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)  

The following example code shows how it was done for [GLDS-72](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/), which was 454 sequencing of the 16S gene using these primers:  
* forward: 5’-AGAGTTTGATCCTGGCTCAG-3’  
* reverse: 5’- GCTGCCTCCCGTAGGAGT-3’  


```
cutadapt -g AGAGTTTGATCCTGGCTCAG -a GCTGCCTCCCGTAGGAGT \
         -o sample-1_trimmed.fastq.gz sample-1_raw.fastq.gz
```

**Parameter Definitions:**

*	`-g` – specifies the forward primer 

*	`-a` – specifies the reverse primer

*	`-o` – specifies output primer-trimmed reads

*	`sample-1_raw.fastq.gz` – positional argument specifying the input reads


**Input files:**

* fastq, compressed or uncompressed (original reads)

**Output files:**

* fastq, compressed or uncompressed (trimmed reads)

<br>

---

## 3. Quality filtering

The following is an example from a [GLDS-72](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/) sample that used 454 sequencing with the following 16S primers:  
* forward: 5’-AGAGTTTGATCCTGGCTCAG-3’
* reverse: 5’- GCTGCCTCCCGTAGGAGT-3’

```
bbduk.sh in=sample-1_trimmed.fastq.gz out=sample-1_filtered.fastq.gz \
         qtrim=r trimq=10 mlf=0.5 minavgquality=20 minlength=50
```

**Parameter Definitions:**

*	`in=` – the input file

*	`out=` – the output file

*	`qtrim=` – specifies the direction to trim low-quality bases from

*	`trimq=` – the minimum quality to keep while trimming low-quality bases

*	`mlf=` – allowed minimum length after trimming is half the initial read size (read filtered out if shorter)

*	`minavgquality=` – minimum average quality allowed after trimming (read filtered out if below)

*	`minlength=` – allowed minimum length (backup to `mlf` setting if starting read is shorter than 100 bps)

**Input file types:**

* fastq, compressed or uncompressed (primer-trimmed reads)

**Output file types:**

* fastq, compressed or uncompressed (filtered reads)

<br>

---

## 4. Filtered Data QC
```
fastqc -o filtered_fastqc_output *filtered.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*filtered.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input files:**

* fastq, compressed or uncompressed (filtered reads)

**Output files:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)

<br>

### Compile Filtered Data QC
```
multiqc -o filtered_multiqc_output  filtered_fastqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`filtered_fastqc_output` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input files:**

* fastqc.zip (FastQC output data)

**Output files:**

* multiqc_report.html (multiqc output html summary)
* multiqc_data (directory containing multiqc output data)

<br>

---

## 5. Generating OTUs and counts per sample

### Dereplicate individual samples
```
vsearch --derep_fulllength sample-1_filtered.fastq.gz --strand both --output sample-1_derep.fasta --sizeout --relabel "sample=sample-1;seq_" 
```

**Parameter Definitions:**

*	`--derep_fulllength` – the vsearch command
*	`sample-1_filtered.fastq.gz` – input file, provided as a positional argument
*  `--strand both` - specifies to check both orientations
*  `--output` - designate the name of the output fasta file
*  `--sizeout` - incorporates abundance information in the sequence header
*  `--relabel` - relabel the headers of the sequences starting with this prefix

**Input files:**

* sample-1_filtered.fastq.gz

**Output files:**

* sample-1_derep.fasta

<br>

### Generate OTUs

#### Combining all individual sample dereplicated sequences for further processing
```bash
cat *_derep.fasta > all-samples-seqs.fasta
rm *_derep.fasta
```

#### Dereplicating combined sequences
```bash
vsearch --derep_fulllength all-samples-seqs.fasta --strand both --output all-samples_derep.fasta --sizein --sizeout
```

**Parameter Definitions:**

*	`--derep_fulllength` – the vsearch command
*	`all-samples-seqs.fasta` – input file, provided as a positional argument
*  `--strand both` - specifies to check both orientations
*  `--output` - designate the name of the output fasta file
*  `--sizein` - take into account abundance information in header
*  `--sizeout` - incorporates abundance information in the sequence header


**Input files:**

* all-samples-seqs.fasta

**Output files:**

* all-samples_derep.fasta


#### Clustering to get representative sequences
```bash
vsearch --cluster_size all-samples_derep.fasta --id 0.97 --strand both --sizein --sizeout --relabel "OTU_" --centroids rep-seqs.fasta
```

**Parameter Definitions:**

*	`--cluster_size ` – the vsearch command
*	`all-samples_derep.fasta` – input file, provided as a positional argument
*  `--id` - specifies the percent identity to cluster at
*  `--strand both` - specifies to consider both orientations
*  `--sizein` - take into account abundance information in header
*  `--sizeout` - incorporates abundance information in the sequence header
*  `--relabel` - relabel the headers of the sequences starting with this prefix
*  `--centroids` - designate the name of the output fasta file holding representative sequences


**Input files:**

* all-samples_derep.fasta

**Output files:**

* rep-seqs.fasta


#### Removing singletons
```bash
vsearch --sortbysize rep-seqs.fasta --minsize 2 --output rep-seqs-no-singletons.fasta
```

**Parameter Definitions:**

*	`--sortbysize ` – the vsearch command
*	`rep-seqs.fasta` – input file, provided as a positional argument
*  `--minsize` - minimum cluster size to be retained
*  `--output` - designate the name of the output fasta file holding filtered representative sequences


**Input files:**

* rep-seqs.fasta

**Output files:**

* rep-seqs-no-singletons.fasta

#### Chimera check and removal

```bash
vsearch --uchime_denovo rep-seqs-no-singletons.fasta --sizein --nonchimeras OTUs.fasta --relabel "OTU_"
```

**Parameter Definitions:**

*	`--uchime_denovo ` – the vsearch command
*	`rep-seqs-no-singletons.fasta` – input file, provided as a positional argument
*  `--sizein` - take into account abundance information in header
*  `--nonchimeras` - designate the name of the output fasta file holding filtered representative sequences
*  `--relabel` - relabel the headers of the sequences starting with this prefix

**Input files:**

* rep-seqs-no-singletons.fasta

**Output files:**

* OTUs.fasta

<br>

### Map reads to OTUs
```bash
vsearch --usearch_global all-samples-seqs.fasta -db OTUs.fasta --sizein --id 0.97 --otutabout - | sed 's/^#OTU ID/OTU_ID/' > counts.tsv
```

**Parameter Definitions:**

*	`--usearch_global ` – the vsearch command
*	`all-samples-seqs.fasta` – input file, provided as a positional argument (note this is our combined individual-dereplicated samples' fasta of reads) 
*  `--db` - specifies the reference OTUs to map to
*  `--id` - specifies the minimum percent identity to enable mapping a read to an OTU
*  `--otutabout` - designates the output of the count table (here going to stout in order to be modified in the subsequent `sed` command
*  `| sed ... > counts.tsv` - renaming the first column header and writing to `counts.tsv`


**Input files:**

* all-samples-seqs.fasta
* OTUs.fasta

**Output files:**

* counts.tsv

<br>

---

## 6. Generating taxonomy and additional outputs

> The following is performed within `R`.


### Assigning taxonomy

Reading in OTU sequences:
```R
dna <- readDNAStringSet("OTUs.fasta")
```

Downloading the reference R taxonomy object:
```R
download.file( url=“http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData”, destfile=“SILVA_SSU_r138_2019.RData”)
```

**Parameter Definitions:**  

*	`download.file()` – the function we are calling, with the following parameters set within it

*	`url=` – specifying the url address of the file to download

*	`destfile=` – specifying the path/name of the file after downloading

<br>

Loading taxonomy object:
```R
load(“SILVA_SSU_r138_2019.RData”)
```

Classifying sequences:
```R
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand=“both”, processors=NULL)
```

**Parameter Definitions:**  

*	`tax_info <-` – specifies the variable that will store the results within in our R environment

*	`IdTaxa()` – the DECIPHER function we are calling, with the following parameters set within it

*	`test=` – specifying the “dna” object created above holding the sequences we want to classify

*	`trainingSet=` – specifying the reference database we downloaded and loaded above

*	`strand=“both”` – specifying to check taxonomy assignment in both orientations

*	`processors=NULL` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

<br>

### Generating and writing outputs
Creating table of taxonomy and setting any that are unclassified as "NA", and writing out:

```R
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

tax_tab <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL
otu_ids <- names(tax_info)
tax_tab <- data.frame("OTU_ID"=otu_ids, tax_tab, check.names=FALSE)

write.table(tax_tab, "taxonomy.tsv", sep = "\t", quote=F, row.names=FALSE)
```

Reading in counts table and generating additional outputs:

```R
otu_tab <- read.table("counts.tsv", sep="\t", header=TRUE, check.names=FALSE)

    # generating and writing out biom file format
biom_object <- make_biom(data=otu_tab, observation_metadata=tax_tab)
write_biom(biom_object, "../Final_Outputs/taxonomy-and-counts.biom")

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, otu_tab)
write.table(tax_and_count_tab, "../Final_Outputs/taxonomy-and-counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)

    # making and writing out final count summary table
cutadapt_tab <- read.table("../Trimmed_Sequence_Data/trimmed-read-counts.tsv", sep="\t", header=TRUE)
bbduk_tab <- read.table("../Filtered_Sequence_Data/filtered-read-counts.tsv", sep="\t", header=TRUE)[,c(1,3)]
otu_tab <- read.table("../Final_Outputs/counts.tsv", sep="\t", header=TRUE, check.names=FALSE, row.names=1)
mapped_sums <- colSums(otu_tab)
mapped_tab <- data.frame(sample=names(mapped_sums), mapped_to_OTUs=mapped_sums, row.names=NULL)

t1 <- merge(cutadapt_tab, bbduk_tab)
count_summary_tab <- merge(t1, mapped_tab)

write.table(count_summary_tab, "../Final_Outputs/read-count-tracking.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

---
---
