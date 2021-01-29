# Kilpisjarvi-MAGs

Bioinformatics workflow used in the article:

> Pessi I. S.,  Viitamäki S., Eronen-Rasimus E., Delmont T. O., Luoto M., Hultman J. 2020. Truncated denitrifiers dominate the denitrification pathway in tundra soil metagenomes. BioRxiv, doi: [10.1101/2020.12.21.419267](https://doi.org/10.1101/2020.12.21.419267).

## Contacts

**Igor S Pessi**  
Postdoctoral Researcher  
[E-mail](mailto:igor.pessi@helsinki.fi)

**Jenni Hultman**  
Principal Investigator  
[E-mail](mailto:jenni.hultman@helsinki.fi)

## Table of contents

1. [Pre-processing of raw data](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/01-pre-processing.md)
2. [Read-based analyses](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/02-read-based.md)
3. [Metagenome assembling](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/03-assembling.md)
4. [Binning of MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/04-MAG-binning.md)
5. [Re-assembly of hybrid MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/05-hybrid-assembling.md)
6. [Working with MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/06-working-with-MAGs.md)
7. [Denitrifier MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/07-denitrifier-MAGs.md)

## Before starting

### You will need to have these softwares installed and in your path

* fastQC v0.11.9: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* multiQC v1.8: https://multiqc.info/
* Cutadapt v1.16: https://cutadapt.readthedocs.io/
* pycoQC v2.5.0.21: https://github.com/a-slide/pycoQC/
* Porechop v0.2.4: https://github.com/rrwick/Porechop/
* seqtk v1.3: https://github.com/lh3/seqtk/
* METAXA v2.2: https://microbiology.se/software/metaxa2/
* DIAMOND v0.9.14: https://github.com/bbuchfink/diamond/
* KEGG database release 86: https://www.genome.jp/kegg/
  * The KEGG database is available for download for paying subscribers only :(
  * Alternatively, you can use their free-of-charge online tool BlastKOALA: https://www.kegg.jp/blastkoala/
* MEGAHIT v1.1.1.2: https://github.com/voutcn/megahit/
* Flye v2.7.1: https://github.com/fenderglass/Flye/
* bowtie v2.3.5: http://bowtie-bio.sourceforge.net/bowtie2/
* SAMtools v1.9: http://www.htslib.org/
* pilon v1.23: https://github.com/broadinstitute/pilon/
* metaQUAST v5.0.2: http://bioinf.spbau.ru/metaquast/
* anvi’o v6.2: https://merenlab.org/software/anvio/
* Unicycler v0.4.8: https://github.com/rrwick/Unicycler/
* GTDB-Tk v1.3.0 and GTDB release 05-RS95: https://gtdb.ecogenomic.org/
* KOfam database release 11/2020: ftp://ftp.genome.jp/pub/db/kofam/
* R: https://www.r-project.org/
* R packages
  * base v3.6.2
  * tidyverse v1.3.0 (https://cran.r-project.org/web/packages/tidyverse)
  * keggR v0.9.1 (https://github.com/igorspp/keggR)
  * multcomp v1.4.15 (https://cran.r-project.org/web/packages/multcomp)
  * vegan v2.5.6 (https://cran.r-project.org/web/packages/vegan)
  * ape v5.3 (https://cran.r-project.org/web/packages/ape)
  * pheatmap v1.0.12 (https://cran.r-project.org/web/packages/pheatmap)
  * VennDiagram v1.6.20 (https://cran.r-project.org/web/packages/VennDiagram)

We have used an Atos BullSequana X400 system running the Red Hat Enterprise Linux Server 7.7 (Maipo).  
You should be able to run the analysis with any UNIX-based OS.

### Define number of threads to use

You should change below to the number of cores available in your system:

```bash
NTHREADS=40
```

### Get sample metadata

We have prepared a metadata file that will help us downloading the FASTQ files and running some of the scripts.  
You can download it [here](https://raw.githubusercontent.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/main/sample_metadata.tsv), or using the command line:

```bash
wget https://raw.githubusercontent.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/main/sample_metadata.tsv
```
