# Pre-processing of raw metagenomic data

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
  curl $FTP_PATH > BASECALLED_NANOPORE/$PREFIX.fastq.gz
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
