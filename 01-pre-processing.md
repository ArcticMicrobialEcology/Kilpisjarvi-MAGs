# Pre-processing of raw metagenomic data

We have a total of 69 samples.  
Most of them were sequenced only with Illumina NextSeq/NovaSeq; two samples were also sequenced with Nanopore MinION.

## Illumina data

The 69 samples were sequenced using a mix of Illumina NextSeq and NovaSeq, with some samples being sequenced with both.  
That is why we have a total of 88 paired FASTQ files.  
After quality control, we will merge the files from the same sample.

### Download raw data from ENA

```bash
mkdir RAW_ILLUMINA

while read SAMPLE RUN R1_FTP_PATH R2_FTP_PATH; do
  curl $R1_FTP_PATH > RAW_ILLUMINA/$SAMPLE.$RUN.R1.fastq.gz
  curl $R2_FTP_PATH > RAW_ILLUMINA/$SAMPLE.$RUN.R2.fastq.gz
done < <(cut -f 1-4 sample_metadata.tsv | sed '1d')
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

while read SAMPLE RUN; do
  cutadapt RAW_ILLUMINA/$SAMPLE.$RUN.R1.fastq.gz \
           RAW_ILLUMINA/$SAMPLE.$RUN.R2.fastq.gz \
           -o TRIMMED_ILLUMINA/$SAMPLE.$RUN.R1.fastq \
           -p TRIMMED_ILLUMINA/$SAMPLE.$RUN.R2.fastq \
           -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
           -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
           -m 50 \
           -j $NTHREADS \
           --nextseq-trim 20
done < <(cut -f 1-2 sample_metadata.tsv | sed '1d')
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

### Pool NextSeq and NovaSeq runs

```bash
mkdir POOLED_ILLUMINA

while read SAMPLE RUN; do
  cat TRIMMED_ILLUMINA/$SAMPLE.$RUN.R1.fastq >> POOLED_ILLUMINA/$SAMPLE.R1.fastq
  cat TRIMMED_ILLUMINA/$SAMPLE.$RUN.R2.fastq >> POOLED_ILLUMINA/$SAMPLE.R2.fastq
done < <(cut -f 1-2 sample_metadata.tsv | sed '1d')
```

### Get taxonomic profiles

Because there are large differences in library size, we are going to resample the dataset to 2,000,000 reads per sample.

```bash
# Resample FASTQ files with seqtk
mkdir RESAMPLED_ILLUMINA

SAMPLES=`cut -f 1 sample_metadata.tsv | sed '1d' | uniq`

for SAMPLE in $SAMPLES; do
  seqtk sample -s100 POOLED_ILLUMINA/$SAMPLE.R1.fastq 2000000 > RESAMPLED_ILLUMINA/$SAMPLE.R1.fastq
  seqtk sample -s100 POOLED_ILLUMINA/$SAMPLE.R2.fastq 2000000 > RESAMPLED_ILLUMINA/$SAMPLE.R2.fastq
done

# Find 16S rRNA gene sequences with METAXA
mkdir METAXA

for SAMPLE in $SAMPLES; do
  metaxa2 -1 RESAMPLED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 RESAMPLED_ILLUMINA/$SAMPLE.R2.fastq \
          -o METAXA/$SAMPLE \
          --align none \
          --graphical F \
          --cpu $NTHREADS \
          --plus
done

# Download SILVA database formatted for mothur
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz
tar zxf silva.nr_v138_1.tgz

# Classify sequences against SILVA with mothur
for SAMPLE in $SAMPLES; do
  mothur "#classify.seqs(fasta=METAXA/$SAMPLE.bacteria.fasta, template=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax)"
  mothur "#classify.seqs(fasta=METAXA/$SAMPLE.archaea.fasta,  template=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax)"
done
```

## Nanopore data

### Download basecalled data from ENA

> **NOTE:**  
> The Nanopore data deposited at ENA has been already basecalled with GPU guppy v4.0.11 and checked with pycoQC v2.5.0.21.  
> These were the commands that were used (DO NOT RUN):
>
>```bash
> for SAMPLE in m11216 m12208; do
>   # Run guppy
>   guppy_basecaller -i RAW_NANOPORE/$SAMPLE/ \
>                    -s BASECALLED_NANOPORE/$SAMPLE \
>                    -c ont-guppy/data/dna_r9.4.1_450bps_hac.cfg \
>                    --device auto \
>                    --qscore_filtering
>
>   # Run pycoQC
>   pycoQC --summary_file BASECALLED_NANOPORE/$SAMPLE/sequencing_summary.txt \
>          --html_outfile BASECALLED_NANOPORE/$SAMPLE/pycoQC_report.html
>
>   # Merge passed reads
>   cat BASECALLED_NANOPORE/$SAMPLE/pass/*.fastq | gzip > BASECALLED_NANOPORE/$SAMPLE.fastq.gz
> done
>```

```bash
mkdir BASECALLED_NANOPORE

curl $FTP_PATH > BASECALLED_NANOPORE/m11216.NANOPORE.fastq.gz
curl $FTP_PATH > BASECALLED_NANOPORE/m12208.NANOPORE.fastq.gz
```

### Trim adapters with Porechop

```bash
mkdir TRIMMED_NANOPORE

for SAMPLE in m11216 m12208; do
  porechop -i BASECALLED_NANOPORE/$SAMPLE.fastq.gz \
           -o TRIMMED_NANOPORE/$SAMPLE.fastq \
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

Continue to [metagenome assembling](02-assembling.md).
