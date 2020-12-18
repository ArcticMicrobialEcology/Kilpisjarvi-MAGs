# Working with MAGs

### Dereplicate MAGs with fastANI

We will dereplicate the Illumina and Nanopore MAGs to create a set of non-redundant MAGs.  
Then we will manually go through the hybrid MAGs and replace the respective Nanopore/Illumina MAGs if they are of better quality.

```bash
printf '%s\t%s\n' name path > MAGS/fastANI_paths.txt

# Get Illumina and Nanopore MAGs
for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
  for MAG in `ls BINNING/$ASSEMBLY/FINAL_MAGS | sed 's/.fa//g'`; do
    FILE=`ls BINNING/$ASSEMBLY/FINAL_MAGS/$MAG.fa`

    printf '%s\t%s\n' $MAG $FILE
  done
done >> fastANI_paths.txt

# Run fastANI
anvi-dereplicate-genomes --fasta-text-file MAGS/fastANI_paths.txt \
                         --output-dir MAGS/FASTANI \
                         --program fastANI \
                         --num-threads $NTHREADS \
                         --skip-fasta-report \
                         --similarity-threshold 0.90

# Get dereplicated MAGs
mkdir NR_MAGS

for MAG in `sed '1d' MAGS/FASTANI/CLUSTER_REPORT.txt | cut -f 3 | sort`; do
  cp BINNING/*/FINAL_MAGS/$MAG.fa NR_MAGS
done
```

### Do phylogenomic analysis with GTDB-Tk

```bash
gtdbtk classify_wf --genome_dir MAGS/NR_MAGS \
                   --out_dir MAGS/GTDB \
                   --extension fa \
                   --cpus $NTHREADS \
                   --pplacer_cpus $NTHREADS
```

### Import non-redundant MAGs to anvi'o

```bash
# Concatenate MAGs
cat MAGS/NR_MAGS/*.fa > MAGS/NR_MAGs.fa

# Build a contigs database
anvi-gen-contigs-database --contigs-fasta MAGS/NR_MAGs.fa \
                          --output-db-path MAGS/CONTIGS.db \
                          --project-name NR_MAGS

# Find single-copy genes with HMMER
anvi-run-hmms --contigs-db MAGS/CONTIGS.db \
              --num-threads $NTHREADS

# Get taxonomy for single copy genes
anvi-run-scg-taxonomy --contigs-db MAGS/CONTIGS.db \
                      --num-threads $NTHREADS

# Map reads with bowtie
mkdir MAGS/MAPPING

bowtie2-build MAGS/NR_MAGs.fa \
              MAGS/MAPPING/NR_MAGs

SAMPLES=`cut -f 1 sample_metadata_illumina.txt | sed '1d' | sort | uniq`

for SAMPLE in $SAMPLES; do
  bowtie2 -1 TRIMMED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 TRIMMED_ILLUMINA/$SAMPLE.R2.fastq \
          -S MAGS/MAPPING/$SAMPLE.sam \
          -x MAGS/MAPPING/NR_MAGs \
          --threads $NTHREADS \
          --no-unal

  samtools view -F 4 -bS MAGS/MAPPING/$SAMPLE.sam | samtools sort > MAGS/MAPPING/$SAMPLE.bam
  samtools index MAGS/MAPPING/$SAMPLE.bam
done

# Build profile databases
for SAMPLE in $SAMPLES; do
  anvi-profile --input-file MAGS/MAPPING/$SAMPLE.bam \
               --output-dir MAGS/PROFILES/$SAMPLE \
               --contigs-db MAGS/CONTIGS.db \
               --num-threads $NTHREADS
done

# Merge profiles
anvi-merge MAGS/PROFILES/*/PROFILE.db \
           --output-dir MAGS/MERGED_PROFILES \
           --contigs-db MAGS/CONTIGS.db

# Create map of splits to MAGs
for SPLIT in `sqlite3 MAGS/CONTIGS.db 'SELECT split FROM splits_basic_info'`; do
  MAG=`echo $SPLIT | awk -F '_' -v OFS='_' '{print $1, $2, $3}'`
  printf '%s\t%s\n' $SPLIT $MAG
done > MAGS/NR_MAGs_splits.txt

# Import collection
anvi-import-collection MAGS/NR_MAGs_splits.txt \
                       --contigs-db MAGS/CONTIGS.db \
                       --pan-or-profile-db MAGS/MERGED_PROFILES/PROFILE.db \
                       --collection-name NR_MAGS
```

### Summarize MAGs

```bash
anvi-summarize --contigs-db MAGS/CONTIGS.db \
               --pan-or-profile-db MAGS/MERGED_PROFILES/PROFILE.db \
               --output-dir MAGS/SUMMARY \
               --collection-name NR_MAGS
```

### Annotate MAGs against COG with DIAMOND

```bash
anvi-run-ncbi-cogs --contigs-db MAGS/CONTIGS.db \
                   --num-threads $NTHREADS
```

### Annotate MAGs against KEGG with DIAMOND

```bash
# Get amino acid sequences for gene calls
anvi-get-sequences-for-gene-calls --contigs-db MAGS/CONTIGS.db \
                                  --output-file MAGS/gene_calls.faa \
                                  --get-aa-sequences

# Run DIAMOND
KEGG_DB_DIR=$HOME/kegg/genes # Change here to the location of the PROKARYOTES.pep.gz file in your system

diamond blastp --query MAGS/gene_calls.faa \
               --out MAGS/gene_calls_KEGG.txt \
               --db $KEGG_DB_DIR/PROKARYOTES \
               --outfmt 6 \
               --max-target-seqs 1 \
               --max-hsps 1 \
               --threads $NTHREADS
```

### Annotate MAGs against KOfam with HMMER

```bash
mkdir MAGS/KOFAM

KOFAM_DB_DIR=$HOME/KOfam # Change here to where you have downloaded and decompressed the KOfam database

# Run hmmer
while read knum threshold score_type profile_type F_measure nseq nseq_used alen mlen eff_nseq re_pos definition; do
  if [[ $score_type == "full" ]]; then
    hmmsearch -T $threshold --cpu $NTHREADS --tblout MAGS/KOFAM/$knum.hmm.txt $KOFAM_DB_DIR/profiles/$knum.hmm  MAGS/gene_calls.faa
  else
    hmmsearch --domT $threshold --cpu $NTHREADS --domtblout MAGS/KOFAM/$knum.hmm.txt $KOFAM_DB_DIR/profiles/$knum.hmm MAGS/gene_calls.faa
  fi
done < <(sed '1d' $KOFAM_DB_DIR/ko_list)

# Concatenate results
while read knum threshold score_type profile_type F_measure nseq nseq_used alen mlen eff_nseq re_pos definition; do
  if [[ $score_type == "full" ]]; then
    sed '/^#/d' MAGS/KOFAM/$knum.hmm.txt |
    tr -s ' ' '\t' |
    cut -f 1,3,5,6
  else
    sed '/^#/d' MAGS/KOFAM/$knum.hmm.txt |
    tr -s ' ' '\t' |
    cut -f 1,4,7,8
  fi
done < <(sed '1d' $KOFAM_DB_DIR/ko_list) > MAGS/gene_calls_KOFAM.txt
```
