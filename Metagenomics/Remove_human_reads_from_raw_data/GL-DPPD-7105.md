# GeneLab removal of human reads from metagenomics datasets

It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in [GeneLab's data repository](https://genelab-data.ndc.nasa.gov/genelab/projects). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. 

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|------:|
|kraken2|`kraken2 -v`|[https://github.com/DerrickWood/kraken2/wiki](https://github.com/DerrickWood/kraken2/wiki)|

## 1. Build kraken2 database

### 1a. Download human reference

```
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
```

**Parameter Definitions:**

* 
* 

**Input data:**

* 

**Output data:**

* 
* 
