# MAG binning

Here we will bin MAGs for each of the four (co-)assemblies separately.

```bash
mkdir BINNING
```

### Define assembly and create list of sample names

```bash
ASSEMBLY=UPLAND_CO   # OR
ASSEMBLY=FEN_CO      # OR
ASSEMBLY=M11216_NANO # OR
ASSEMBLY=M12208_NANO

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'M11216_NANO' ]]; then
  SAMPLES=`awk -F '\t' '{if ($5 == "upland") {print $1}}' sample_metadata.tsv | uniq`
fi

if [[ $ASSEMBLY == 'FEN_CO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  SAMPLES=`awk -F '\t' '{if ($5 == "fen") {print $1}}' sample_metadata.tsv | uniq`
fi
```

### Rename contigs and select those >2,500 bp

```bash
mkdir BINNING/$ASSEMBLY

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  FASTA_FILE=ASSEMBLIES/$ASSEMBLY/final.contigs.fa
fi

if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  FASTA_FILE=ASSEMBLIES/$ASSEMBLY/pilon.fasta
fi

anvi-script-reformat-fasta $FASTA_FILE \
                           --output-file BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                           --report-file BINNING/$ASSEMBLY/CONTIGS_reformat.txt \
                           --prefix $ASSEMBLY \
                           --min-len 2500 \
                           --simplify-names
```

### Build a contigs database

```bash
anvi-gen-contigs-database --contigs-fasta BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                          --output-db-path BINNING/$ASSEMBLY/CONTIGS.db \
                          --project-name $ASSEMBLY \
                          --num-threads $NTHREADS
```

### Find single-copy genes with HMMER

```bash
anvi-run-hmms --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
              --num-threads $NTHREADS
```

### Get taxonomy for single copy genes against GTDB

```bash
anvi-run-scg-taxonomy --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                      --num-threads $NTHREADS
```

### Map reads with bowtie

```bash
mkdir BINNING/$ASSEMBLY/MAPPING

bowtie2-build BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
              BINNING/$ASSEMBLY/MAPPING/contigs

for SAMPLE in $SAMPLES; do
  bowtie2 -1 POOLED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 POOLED_ILLUMINA/$SAMPLE.R2.fastq \
          -S BINNING/$ASSEMBLY/MAPPING/$SAMPLE.sam \
          -x BINNING/$ASSEMBLY/MAPPING/contigs \
          --threads $NTHREADS \
          --no-unal

  samtools view -F 4 -bS BINNING/$ASSEMBLY/MAPPING/$SAMPLE.sam | samtools sort > BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam
  samtools index BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam
done
```

### Build profile databases

```bash
mkdir BINNING/$ASSEMBLY/PROFILES

for SAMPLE in $SAMPLES; do
  anvi-profile --input-file BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam \
               --output-dir BINNING/$ASSEMBLY/PROFILES/$SAMPLE \
               --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
               --num-threads $NTHREADS
done
```

### Merge profiles

```bash
if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  anvi-merge BINNING/$ASSEMBLY/PROFILES/*/PROFILE.db \
             --output-dir BINNING/$ASSEMBLY/MERGED_PROFILES \
             --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
             --skip-hierarchical-clustering
fi

if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  anvi-merge BINNING/$ASSEMBLY/PROFILES/*/PROFILE.db \
             --output-dir BINNING/$ASSEMBLY/MERGED_PROFILES \
             --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
             --enforce-hierarchical-clustering
fi
```

### Bin MAGs

The co-assemblies are really huge and thus impossible to load in the interactive interface for binning.  
We will then first use CONCOCT to pre-cluster the contigs into 100 clusters, which will then be suitable for the interactive interface.

```bash
# Illumina co-assemblies
if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  # Pre-cluster with CONCOCT
  anvi-cluster-contigs --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                       --profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
                       --collection-name CONCOCT \
                       --num-threads $NTHREADS \
                       --driver concoct \
                       --clusters 100 \
                       --just-do-it

  # Create individual contigs and profile databases
  mkdir BINNING/$ASSEMBLY/CONCOCT_SPLIT

  anvi-split --pan-or-profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
             --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
             --output-dir BINNING/$ASSEMBLY/CONCOCT_SPLIT \
             --collection-name CONCOCT \
             --skip-variability-tables

  CLUSTERS=`ls BINNING/$ASSEMBLY/CONCOCT_SPLIT`

  # Bin MAGs
  for CLUSTER in $CLUSTERS; do
    anvi-interactive --profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                     --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db
  done

  # Call MAGs
  for CLUSTER in $CLUSTERS; do
    anvi-rename-bins --profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                     --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                     --report-file BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/renamed_bins.txt \
                     --collection-to-read DEFAULT \
                     --collection-to-write FINAL \
                     --min-completion-for-MAG 50 \
                     --max-redundancy-for-MAG 10 \
                     --prefix "$ASSEMBLY"_"$CLUSTER" \
                     --call-MAGs
  done

  # Refine MAGs
  for CLUSTER in $CLUSTERS; do
    for MAG in `sqlite3 BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
      anvi-refine --profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                  --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                  --collection-name FINAL \
                  --bin-id $MAG
    done
  done

  # Summarise MAGs
  for CLUSTER in $CLUSTERS; do
    anvi-summarize --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                   --pan-or-profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                   --output-dir BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER.SUMMARY \
                   --collection-name FINAL
  done
fi

# Nanopore individual assemblies
if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  # Bin MAGs
  anvi-interactive --profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
                   --contigs-db BINNING/$ASSEMBLY/CONTIGS.db

  # Call MAGs
  anvi-rename-bins --profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
                   --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                   --report-file BINNING/$ASSEMBLY/renamed_bins.txt \
                   --collection-to-read DEFAULT \
                   --collection-to-write FINAL \
                   --min-completion-for-MAG 50 \
                   --max-redundancy-for-MAG 10 \
                   --prefix $ASSEMBLY \
                   --call-MAGs

  # Refine MAGs
  for MAG in `sqlite3 BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
    anvi-refine --profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
                --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                --collection-name FINAL \
                --bin-id $MAG
  done

  # Summarise MAGs
  anvi-summarize --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                 --pan-or-profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
                 --output-dir BINNING/$ASSEMBLY/SUMMARY \
                 --collection-name FINAL
fi
```

### Rename MAGs

```bash
mkdir BINNING/FINAL_MAGs

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'M11216_NANO' ]]; then
  PREFIX=KUL
  COUNT=`ls BINNING/FINAL_MAGs | grep KUL | wc -l`
fi

if [[ $ASSEMBLY == 'FEN_CO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  PREFIX=KWL
  COUNT=`ls BINNING/FINAL_MAGs | grep KWL | wc -l`
fi

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  for CLUSTER in $CLUSTERS; do
    for MAG in `sqlite3 BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
      COUNT=`expr $COUNT + 1`
      NEWMAG=`printf '%s_%04d\n' $PREFIX $COUNT`

      printf '%s\t%s\n' $MAG $NEWMAG >> BINNING/FINAL_MAGs/$ASSEMBLY.renamed_MAGs.txt

      anvi-script-reformat-fasta BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER.SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                                 --output-file BINNING/FINAL_MAGs/$NEWMAG.fa \
                                 --prefix $NEWMAG \
                                 --simplify-names
    done
  done
fi

if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  for MAG in `sqlite3 BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
    COUNT=`expr $COUNT + 1`
    NEWMAG=`printf '%s_%04d\n' $PREFIX $COUNT`

    printf '%s\t%s\n' $MAG $NEWMAG >> BINNING/FINAL_MAGs/$ASSEMBLY.renamed_MAGs.txt

    anvi-script-reformat-fasta BINNING/$ASSEMBLY/SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               --output-file BINNING/FINAL_MAGs/$NEWMAG.fa \
                               --prefix $NEWMAG \
                               --simplify-names
  done
fi
```

## Next step

Continue to [gene-centric analyses](04-gene.centric.md).
