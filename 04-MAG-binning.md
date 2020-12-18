# MAG binning

Here we will bin MAGs for each of the four (co-)assemblies separately.
This is done by assigning the variables $ASSEMBLY and $SAMPLES.

### Define assembly and create list of sample names

```bash
ASSEMBLY=UPLAND_CO   # OR
ASSEMBLY=FEN_CO      # OR
ASSEMBLY=M11216_NANO # OR
ASSEMBLY=M12208_NANO

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'M11216_NANO' ]]; then
  SAMPLES=`cat UPLAND_SAMPLES.txt`
fi

if [[ $ASSEMBLY == 'FEN_CO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  SAMPLES=`cat FEN_SAMPLES.txt`
fi
```

### Rename contigs and select those >2,500 bp

```bash
mkdir -p BINNING/$ASSEMBLY

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  anvi-script-reformat-fasta ASSEMBLIES/$ASSEMBLY/final.contigs.fa \
                             --output-file BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                             --report-file BINNING/$ASSEMBLY/CONTIGS_reformat.txt \
                             --prefix $ASSEMBLY \
                             --min-len 2500 \
                             --simplify-names
fi

if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  anvi-script-reformat-fasta ASSEMBLIES/$ASSEMBLY/pilon.fasta \
                             --output-file BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                             --report-file BINNING/$ASSEMBLY/CONTIGS_reformat.txt \
                             --prefix $ASSEMBLY \
                             --min-len 2500 \
                             --simplify-names
fi
```

### Prepare files for anvi'o

```bash
# Build a contigs database
anvi-gen-contigs-database --contigs-fasta BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                          --output-db-path BINNING/$ASSEMBLY/CONTIGS.db \
                          --project-name $ASSEMBLY

# Find single-copy genes with HMMER
anvi-run-hmms --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
              --num-threads $NTHREADS

# Get taxonomy for single copy genes
anvi-run-scg-taxonomy --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                      --num-threads $NTHREADS

# Map reads with bowtie
mkdir BINNING/$ASSEMBLY/MAPPING

bowtie2-build BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
              BINNING/$ASSEMBLY/MAPPING/contigs

for SAMPLE in $SAMPLES; do
  bowtie2 -1 TRIMMED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 TRIMMED_ILLUMINA/$SAMPLE.R2.fastq \
          -S BINNING/$ASSEMBLY/MAPPING/$SAMPLE.sam \
          -x BINNING/$ASSEMBLY/MAPPING/contigs \
          --threads $NTHREADS \
          --no-unal

  samtools view -F 4 -bS BINNING/$ASSEMBLY/MAPPING/$SAMPLE.sam | samtools sort > BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam
  samtools index BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam
done

# Build profile databases
mkdir BINNING/$ASSEMBLY/PROFILES

for SAMPLE in $SAMPLES; do
  anvi-profile --input-file BINNING/$ASSEMBLY/MAPPING/$SAMPLE.bam \
               --output-dir BINNING/$ASSEMBLY/PROFILES/$SAMPLE \
               --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
               --num-threads $NTHREADS
done

# Merge profiles
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

### Bin MAGs with anvi'o

```bash
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
  anvi-split --pan-or-profile-db BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db \
             --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
             --collection-name CONCOCT \
             --output-dir BINNING/$ASSEMBLY/CONCOCT_SPLIT \
             --skip-variability-tables

  # Bin MAGs
  for CLUSTER in `seq 1 100 | awk -v OFS='_' '{print "Bin", $0}'`; do
    anvi-interactive --profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                     --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db
  done

  # Call MAGs
  for CLUSTER in `seq 1 100 | awk -v OFS='_' '{print "Bin", $0}'`; do
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
  for CLUSTER in `seq 1 100 | awk -v OFS='_' '{print "Bin", $0}'`; do
    for MAG in `sqlite3 BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
      anvi-refine --profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                  --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                  --collection-name FINAL \
                  --bin-id $MAG
    done
  done

  # Summarise MAGs
  for CLUSTER in `seq 1 100 | awk -v OFS='_' '{print "Bin", $0}'`; do
    anvi-summarize --contigs-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                   --pan-or-profile-db BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                   --output-dir BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER.SUMMARY \
                   --collection-name FINAL
  done
fi

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
mkdir BINNING/$ASSEMBLY/FINAL_MAGS

if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
  if [[ $ASSEMBLY == 'UPLAND_CO' ]]; then
    PREFIX=KUL
  fi

  if [[ $ASSEMBLY == 'FEN_CO' ]]; then
    PREFIX=KWL
  fi

  COUNT=0

  for CLUSTER in `seq 1 100 | awk -v OFS='_' '{print "Bin", $0}'`; do
    for MAG in `sqlite3 BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
      COUNT=`expr $COUNT + 1`
      NEWMAG=`printf '%s_%04d\n' $PREFIX $COUNT`

      printf '%s\t%s\n' $MAG $NEWMAG >> BINNING/$ASSEMBLY/renamed_MAGs.txt

      anvi-script-reformat-fasta BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER.SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                                 --output-file BINNING/$ASSEMBLY/FINAL_MAGS/$NEWMAG.fa \
                                 --prefix $NEWMAG \
                                 --simplify-names
    done
  done
fi

if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
  if [[ $ASSEMBLY == 'M11216_NANO' ]]; then
    PREFIX=KUL
    COUNT=`ls BINNING/UPLAND_CO/FINAL_MAGS/ | egrep -o [0-9]{4} | sort -r | head -n 1`
  fi

  if [[ $ASSEMBLY == 'M12208_NANO' ]]; then
    PREFIX=KWL
    COUNT=`ls BINNING/FEN_CO/FINAL_MAGS/ | egrep -o [0-9]{4} | sort -r | head -n 1`
  fi

  for MAG in `sqlite3 BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | grep MAG`; do
    COUNT=`expr $COUNT + 1`
    NEWMAG=`printf '%s_%04d\n' $PREFIX $COUNT`

    printf '%s\t%s\n' $MAG $NEWMAG >> BINNING/$ASSEMBLY/renamed_MAGs.txt

    anvi-script-reformat-fasta BINNING/$ASSEMBLY/SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               --output-file BINNING/$ASSEMBLY/FINAL_MAGS/$NEWMAG.fa \
                               --prefix $NEWMAG \
                               --simplify-names
  done
fi
```

## Next step

Continue to the [re-assembly of hybrid MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/05-hybrid-assembling.md).
