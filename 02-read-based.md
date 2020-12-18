# Read-based analyses

### Resample FASTQ files with seqtk

Because there are large differences in library size, we are going to resample the dataset to 2,000,000 reads per sample.

```bash
mkdir RESAMPLED_ILLUMINA

SAMPLES=`cut -f 1 sample_metadata_illumina.txt | sed '1d' | sort | uniq`

for SAMPLE in $SAMPLES; do
  seqtk sample -s100 POOLED_ILLUMINA/$SAMPLE.R1.fastq 2000000 > RESAMPLED_ILLUMINA/$SAMPLE.R1.fastq
  seqtk sample -s100 POOLED_ILLUMINA/$SAMPLE.R2.fastq 2000000 > RESAMPLED_ILLUMINA/$SAMPLE.R2.fastq
done
```

### Annotate reads with METAXA

```bash
mkdir METAXA

# Run METAXA
for SAMPLE in $SAMPLES; do
  metaxa2 -1 RESAMPLED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 RESAMPLED_ILLUMINA/$SAMPLE.R2.fastq \
          -o METAXA/$SAMPLE \
          --align none \
          --graphical F \
          --cpu $NTHREADS \
          --plus

  metaxa2_ttt -i METAXA/$SAMPLE.taxonomy.txt \
              -o METAXA/$SAMPLE
done

# Summarise results
for NUM in seq 1 8; do
  metaxa2_dc METAXA/*level_"$NUM".txt \
             -o METAXA/summary_level_"$NUM".txt
done
```

### Annotate reads against the KEGG database with DIAMOND

```bash
mkdir KEGG

# Convert reads to FASTA, rename headers and concatenate R1+R2
for SAMPLE in $SAMPLES; do
  seqtk seq -A RESAMPLED_ILLUMINA/$SAMPLE.R1.fastq |
  awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R1", "READ", ++i; next}{print}' >> RESAMPLED_ILLUMINA/$SAMPLE.fasta

  seqtk seq -A RESAMPLED_ILLUMINA/$SAMPLE.R2.fastq |
  awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R2", "READ", ++i; next}{print}' >> RESAMPLED_ILLUMINA/$SAMPLE.fasta
done

# Format KEGG db for use with DIAMOND
KEGG_DB_DIR=$HOME/kegg/genes # Change here to the location of the PROKARYOTES.pep.gz file in your system

diamond makedb -in $KEGG_DB_DIR/PROKARYOTES.pep.gz \
               -db $KEGG_DB_DIR/PROKARYOTES

# Run DIAMOND
for SAMPLE in $SAMPLES; do
  diamond blastx --query RESAMPLED_ILLUMINA/$SAMPLE.fasta \
                 --out KEGG/$SAMPLE.txt \
                 --db $KEGG_DB_DIR/PROKARYOTES \
                 --outfmt 6 \
                 --evalue 0.00001 \
                 --id 60 \
                 --max-target-seqs 1 \
                 --max-hsps 1 \
                 --threads $NTHREADS
done
```


## Next step

Continue to [metagenome assembling](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/03-assembling.md) or to the [read-based R analyses](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/02-read-based.R).
