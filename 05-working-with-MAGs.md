# Working with MAGs

### Import MAGs to anvi'o

```bash
mkdir MAGs

# Concatenate MAGs
cat BINNING/FINAL_MAGs/*.fa > MAGs/FINAL_MAGs.fa

# Build a contigs database
anvi-gen-contigs-database --contigs-fasta MAGs/FINAL_MAGs.fa \
                          --output-db-path MAGs/CONTIGS.db \
                          --project-name FINAL_MAGs \
                          --num-threads $NTHREADS

# Find single-copy genes with HMMER
anvi-run-hmms --contigs-db MAGs/CONTIGS.db \
              --num-threads $NTHREADS

# Get taxonomy for single copy genes
anvi-run-scg-taxonomy --contigs-db MAGs/CONTIGS.db \
                      --num-threads $NTHREADS
```

### Create a database of redundant MAGs

```bash
# Build an empty profile database
anvi-profile --contigs-db MAGs/CONTIGS.db \
             --output-dir MAGs/PROFILE \
             --sample-name EMPTY \
             --blank-profile \
             --skip-hierarchical-clustering

# Create map of splits to MAGs
for SPLIT in `sqlite3 MAGs/CONTIGS.db 'SELECT split FROM splits_basic_info'`; do
  MAG=`echo $SPLIT | sed -E 's/_[0-9]+_split_[0-9]+//'`
  printf '%s\t%s\n' $SPLIT $MAG
done > MAGs/FINAL_MAGs_splits.txt

# Import collection
anvi-import-collection MAGs/FINAL_MAGs_splits.txt \
                       --contigs-db MAGs/CONTIGS.db \
                       --pan-or-profile-db MAGs/PROFILE/PROFILE.db \
                       --collection-name FINAL_MAGs
```

### Summarize redundant MAGs

```bash
anvi-summarize --contigs-db MAGs/CONTIGS.db \
               --pan-or-profile-db MAGs/PROFILE/PROFILE.db \
               --output-dir MAGs/SUMMARY \
               --collection-name FINAL_MAGs
```

### Do phylogenomic analysis with GTDB-Tk

```bash
gtdbtk classify_wf --genome_dir BINNING/FINAL_MAGS \
                   --out_dir MAGs/GTDB \
                   --extension fa \
                   --cpus $NTHREADS \
                   --pplacer_cpus $NTHREADS
```

### Annotate MAGs (general annotation)

```bash
# Annotate against KEGG KOfams
anvi-run-kegg-kofams --contigs-db MAGs/CONTIGS.db \
                     --num-threads $NTHREADS

# Estimate metabolism based on KEGG KOfams
anvi-estimate-metabolism --contigs-db MAGs/CONTIGS.db \
                         --profile-db MAGs/PROFILE/PROFILE.db \
                         --collection-name FINAL_MAGs

# Export amino acid sequences
anvi-get-sequences-for-gene-calls --contigs-db MAGs/CONTIGS.db \
                                  --output-file MAGs/gene_calls.faa \
                                  --get-aa-sequences 

# Annotate MAGs against dbCAN
run_dbcan.py MAGs/gene_calls.faa protein \
             --dia_cpu $NTHREADS \
             --hmm_cpu $NTHREADS \
             --hotpep_cpu $NTHREADS \
             --out_dir MAGs/DBCAN \
             --db_dir dbcan-db
```

### Annotate MAGs (denitrification genes)

```bash
mkdir DENITRIFIERS

# Annotate genes with KofamScan
KOFAM_DB_DIR=$HOME/KOFAM_DB # Change here to the location of the KOfam database in your system

for KO in K00368 K15864 K04561 K00376; do
  exec_annotation MAGs/gene_calls.faa \
                  --profile $KOFAM_DB_DIR/profiles/$KO.hmm \
                  --ko-list $KOFAM_DB_DIR/ko_list \
                  --format detail-tsv \
                  --cpu $NTHREADS | grep ^* > DENITRIFIERS/KOFAM_SCAN/$KO.txt
done

# Align with MAFFT
for KO in K00368 K15864 K04561 K00376; do
  cut -f 2 DENITRIFIERS/KOFAM_SCAN/$KO.txt | sort | uniq > DENITRIFIERS/$KO/genes.txt
  
  seqtk subseq MAGs/gene_calls.faa DENITRIFIERS/$KO/genes.txt > DENITRIFIERS/$KO/genes.faa

  mafft --auto \
        --reorder \
        --thread $NTHREADS \
        DENITRIFIERS/$KO/genes.faa > DENITRIFIERS/$KO/genes.aln.faa
done
```

### Dereplicate MAGs with fastANI

We will dereplicate the Illumina and Nanopore MAGs to create a set of non-redundant MAGs.  

```bash
# Make list of Illumina and Nanopore MAGs
printf '%s\t%s\t%s\t%s\t%s\n' name bin_id collection_id profile_db_path contigs_db_path > MAGs/internal_genomes.txt

for MAG in `sqlite3 MAGs/PROFILE/PROFILE.db 'SELECT bin_name FROM collections_bins_info'`; do
  printf '%s\t%s\t%s\t%s\t%s\n' $MAG $MAG FINAL_MAGs $PWD/MAGs/PROFILE/PROFILE.db $PWD/MAGs/CONTIGS.db
done >> MAGs/internal_genomes.txt

# Run fastANI
anvi-dereplicate-genomes --internal-genomes MAGs/internal_genomes.txt \
                         --output-dir MAGs/NR_MAGs \
                         --program fastANI \
                         --similarity-threshold 0.90 \
                         --representative-method Qscore \
                         --num-threads $NTHREADS
```

### Get the relative abundance of non-redundant MAGs with coverM

```bash
coverm genome --coupled POOLED_ILLUMINA/*.fastq \
              --genome-fasta-files MAGs/NR_MAGs/GENOMES/*.fa \
              --output-file MAGs/NR_MAGs/coverM_relabund.txt \
              --output-format sparse \
              --min-covered-fraction 0 \
              --threads $NTHREADS
```
