# GeneLab removal of human reads from metagenomic datasets

It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in [GeneLab's data repository](https://genelab-data.ndc.nasa.gov/genelab/projects). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. 

The Snakemake workflow we use to do this is included above, system-specific variables should be set in the "config.yaml" file. The database we use was built with kraken2 v2.1.1 as detailed below, and by default will be downloaded to run with the above workflow (it's ~4.3 GB uncompressed). 

See the "config.yaml" file for more details on usage. An quick example can be run as-is with the files included in this repository.

---

## Kraken2 human database build

> The following was performed on 29-Nov-2020 with kraken v2.1.1.

**Downloading human reference (takes ~2 minutes as run here):**

```bash
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
```

**Downloading NCBI taxonomy info needed (takes ~10 minutes):**

```bash
kraken2-build --download-taxonomy --db kraken2-human-db/
```

**Building database (takes ~20 minutes as run here):**

```bash
kraken2-build --build --db kraken2-human-db/ --threads 30
```

**Removing intermediate files:**

```bash
kraken2-build --clean --db kraken2-human-db/
```

---

## Download database as built on 29-Nov-2020
The reference database 3GB compressed and ~4.3GB uncompressed. It can be downloaded and unpacked with the following:

```bash
curl -L -o kraken2-human-db.tar.gz https://ndownloader.figshare.com/files/25627058

tar -xzvf kraken2-human-db.tar.gz
```
