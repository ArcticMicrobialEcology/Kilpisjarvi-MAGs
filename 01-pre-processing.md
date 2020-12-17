# Pre-processing of raw metagenomic data

## Before starting

### You will need to have these softwares installed and in your path

* fastQC v0.11.9: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* multiQC v1.8: https://multiqc.info/
* Cutadapt v1.16: https://cutadapt.readthedocs.io/
* pycoQC v2.5.0.21: https://github.com/a-slide/pycoQC/
* Porechop v0.2.4: https://github.com/rrwick/Porechop/
* seqtk v1.3: https://github.com/lh3/seqtk/
* METAXA v2.2: https://microbiology.se/software/metaxa2/
* DIAMOND v0.9.25.126: https://github.com/bbuchfink/diamond/
* KEGG database release 86: https://www.genome.jp/kegg/
  * The KEGG database is available for download for paying subscribers only :(
  * Alternatively, you can use their free-of-charge online tool BlastKOALA: https://www.kegg.jp/blastkoala/

We have used an Atos BullSequana X400 system running the Red Hat Enterprise Linux Server 7.7 (Maipo).  
You should be able to run the analysis with any UNIX-based OS.

### Define number of threads to use

You should change below to the number of cores available in your system:

```bash
NTHREADS=40
```

### Get sample metadata

We have prepared two metadata files that will be used to download the FASTQ files and analyse the data in R.  
You can download them [here](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/sample_metadata_illumina.txt) and [here](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/sample_metadata_nanopore.txt), or using the command line:

```bash
wget https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/sample_metadata_illumina.txt
wget https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/sample_metadata_nanopore.txt
```


## Illumina data

### Download raw data from ENA

```bash
mkdir RAW_ILLUMINA

while read PREFIX R1_FTP_PATH R2_FTP_PATH; do
  curl $R1_FTP_PATH > RAW_ILLUMINA/$PREFIX.R1.fastq.gz
  curl $R2_FTP_PATH > RAW_ILLUMINA/$PREFIX.R2.fastq.gz
done < <(cut -f 6,7,8 sample_metadata_illumina.txt | sed '1d')
```

### Check raw data with fastQC and multiQC

```bash
mkdir RAW_ILLUMINA/FASTQC

fastqc RAW_ILLUMINA/*.fastq.gz \
       --outdir RAW_ILLUMINA/FASTQC \
       --threads $NTHREADS

multiqc RAW_ILLUMINA/FASTQC \
        --outdir RAW_ILLUMINA/MULTIQC \
        --interactive
```

### Trim adaptors and do quality filtering with Cutadapt

```bash
mkdir TRIMMED_ILLUMINA

PREFIX=`cut -f 6 sample_metadata_illumina.txt | sed '1d'`

for FILE in $PREFIX; do
  cutadapt RAW_ILLUMINA/$FILE.R1.fastq.gz \
           RAW_ILLUMINA/$FILE.R2.fastq.gz \
           -o TRIMMED_ILLUMINA/$FILE.R1.fastq \
           -p TRIMMED_ILLUMINA/$FILE.R2.fastq \
           -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
           -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
           -m 50 \
           -j $NTHREADS \
           --nextseq-trim 20
done
```

### Check trimmed data with fastQC and multiQC

```bash
mkdir TRIMMED_ILLUMINA/FASTQC

fastqc TRIMMED_ILLUMINA/*.fastq \
       --outdir TRIMMED_ILLUMINA/FASTQC \
       --threads $NTHREADS

multiqc TRIMMED_ILLUMINA/FASTQC \
        --outdir TRIMMED_ILLUMINA/MULTIQC \
        --interactive
```

### Pool NextSeq and NovaSeq data for each sample

```bash
mkdir POOLED_ILLUMINA

while read SAMPLE PREFIX; do
  cat TRIMMED_ILLUMINA/$PREFIX.R1.fastq >> POOLED_ILLUMINA/$SAMPLE.R1.fastq
  cat TRIMMED_ILLUMINA/$PREFIX.R2.fastq >> POOLED_ILLUMINA/$SAMPLE.R2.fastq
done < <(cut -f 1,6 sample_metadata_illumina.txt | sed '1d')
```


## Nanopore data


### Download raw (basecalled) data from ENA

> **NOTE:**  
> The Nanopore data deposited at ENA has been already basecalled with GPU guppy v4.0.11 and checked with pycoQC v2.5.0.21.  
> These were the commands that were used (DO NOT RUN):
>
>```bash
> for SAMPLE in m11216 m12208; do
>   # Run guppy
>   guppy_basecaller -i RAW_NANOPORE/$SAMPLE/*/fast5_skip \
>                    -s BASECALLED_NANOPORE/$SAMPLE \
>                    -c ont-guppy/data/dna_r9.4.1_450bps_hac.cfg \
>                    --device auto \
>                    --qscore_filtering
>
>   # Run pycoQC
>   pycoQC --summary_file BASECALLED_NANOPORE/$SAMPLE/sequencing_summary.txt \
>          --html_outfile BASECALLED_NANOPORE/PYCOQC/$SAMPLE.html
>
>   # Merge passed reads
>   cat BASECALLED_NANOPORE/$SAMPLE/pass/*.fastq | gzip > BASECALLED_NANOPORE/$SAMPLE.fastq.gz
> done
>```


```bash
mkdir BASECALLED_NANOPORE

while read PREFIX FTP_PATH; do
  curl $R1_FTP_PATH > BASECALLED_NANOPORE/$PREFIX.R1.fastq.gz
  curl $R2_FTP_PATH > BASECALLED_NANOPORE/$PREFIX.R2.fastq.gz
done < <(cut -f 6,7,8 sample_metadata_nanopore.txt | sed '1d')
```

### Trim adapters with Porechop

```bash
mkdir TRIMMED_NANOPORE

PREFIX=`cut -f 6 sample_metadata_nanopore.txt | sed '1d'`

for FILE in $PREFIX; do
  porechop -i BASECALLED_NANOPORE/$FILE.fastq.gz \
           -o TRIMMED_NANOPORE/$FILE.fastq \
           --discard_middle \
           --threads $NTHREADS
done
```

### Check trimmed data with fastQC and multiQC

```bash
mkdir TRIMMED_NANOPORE/FASTQC

fastqc TRIMMED_NANOPORE/*.fastq \
       --outdir TRIMMED_NANOPORE/FASTQC \
       --threads $NTHREADS

multiqc TRIMMED_NANOPORE/FASTQC \
        --outdir TRIMMED_NANOPORE/MULTIQC \
        --interactive
```

## Next step

Continue to the [read-based analyses](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/02-read-based.md).
