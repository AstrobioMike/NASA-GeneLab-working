# GeneLab removal of human reads from metagenomics datasets

It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in [GeneLab's data repository](https://genelab-data.ndc.nasa.gov/genelab/projects). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. 

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|------:|
|kraken2|`kraken2 -v`|[https://github.com/DerrickWood/kraken2/wiki](https://github.com/DerrickWood/kraken2/wiki)|

## 1. Build kraken2 database

```bash
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
kraken2-build --download-taxonomy --db kraken2-human-db/
kraken2-build --build --db kraken2-human-db/ --threads 30
kraken2-build --clean --db kraken2-human-db/
```

**Parameter Definitions:**

* `--download-library` - specifies the references to download (here just the human reference genome)
* `--db` - specifies the directory we are putting the database in
* `--threads` - specifies the number of threads to use
* `--no-masking` - prevents [masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences) of low-complexity sequences
* `--download-taxonomy` - downloads taxonomic mapping information
* `--build` - specifies to construct kraken2-formatted database
* `--clean` - specifies to remove unnecessarily intermediate files

**Input data:**

* None

**Output data:**

* kraken2 database files (hash.k2d, opts.k2d, and taxo.k2d)

---

## 2. Filter out human-classified reads

**Example if paired-end reads**

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names --paired \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1_R#.fastq sample-1-R1.fq.gz sample-1-R2.fq.gz
        
# renaming and gzipping output files
mv sample-1_R_1.fastq sample-1-R1-human-reads-removed.fastq && gzip sample-1-R1-human-reads-removed.fastq
mv sample-1_R_2.fastq sample-1-R2-human-reads-removed.fastq && gzip sample-1-R2-human-reads-removed.fastq
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--paired` - specifies input reads are paired-end
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* `--unclassified-out` - name of output files of reads that do were not classified (the `#` symbol gets replaced with "_1" and "_2" in the output file names)
* last two positional arguments are the input read files

**Input data:**

* sample-1-R1.fq.gz (gzipped forward-reads fastq file)
* sample-1-R2.fq.gz (gzipped reverse-reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
* sample-1-R1-human-reads-removed.fastq.gz (human-read removed, gzipped forward-reads fastq file)
* sample-1-R2-human-reads-removed.fastq.gz (human-read removed, gzipped reverse-reads fastq file)

**Example if single-end reads**

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1-human-reads-removed.fastq sample-1.fq.gz

# gzipping output file
gzip sample-1-human-reads-removed.fastq
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* `--unclassified-out` - name of output files of reads that do were not classified (the `#` symbol gets replaced with "_1" and "_2" in the output file names)
* last positional argument is the input read file

**Input data:**

* sample-1.fq.gz (gzipped reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
* sample-1-human-reads-removed.fastq.gz (human-read removed, gzipped reads fastq file)

---
