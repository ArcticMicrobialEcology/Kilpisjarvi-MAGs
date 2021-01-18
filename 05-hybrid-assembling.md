# Re-assembly of hybrid MAGs

### Match Illumina and Nanopore MAGs with fastANI

```bash
mkdir HYBRID_MAGS

printf '%s\t%s\n' name path > HYBRID_MAGS/fastANI_paths.txt

# Get Illumina MAGs
for ASSEMBLY in UPLAND_CO FEN_CO; do
  for MAG in `cut -f2 BINNING/$ASSEMBLY/renamed_MAGs.txt`; do
    FILE=`ls BINNING/FINAL_MAGS/$MAG.fa`

    printf '%s\t%s\n' $MAG $FILE
  done
done >> HYBRID_MAGS/fastANI_paths.txt

# Get Nanopore bins
for ASSEMBLY in M11216_NANO M12208_NANO; do
  for BIN in `sqlite3 BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db 'SELECT bin_name FROM collections_bins_info WHERE collection_name LIKE "FINAL"' | sort`; do
    FILE=`ls BINNING/$ASSEMBLY/SUMMARY/bin_by_bin/$BIN/$BIN-contigs.fa`

    printf '%s\t%s\n' $BIN $FILE
  done
done >> HYBRID_MAGS/fastANI_paths.txt

# Compute ANI
anvi-dereplicate-genomes --fasta-text-file HYBRID_MAGS/fastANI_paths.txt \
                         --output-dir HYBRID_MAGS/FASTANI \
                         --program fastANI \
                         --num-threads $NTHREADS \
                         --similarity-threshold 0.90 \
                         --skip-fasta-report
```

### Re-assemble pairs of Illumina/Nanopore MAGs

Here we will go manually through the CLUSTER_REPORT.txt file generated in the previous step and find the matching pairs of Illumina and Nanopore MAGs.  
We will then prepare a file like this (without the headers):

Illumina MAG | Illumina path                                                                                                              | Nanopore path                                                                                 |
------------ | -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
KUL_0013     | BINNING/UPLAND_CO/CONCOCT_SPLIT/Bin_46.SUMMARY/bin_by_bin/UPLAND_CO_Bin_46_MAG_00002/UPLAND_CO_Bin_46_MAG_00002-contigs.fa | BINNING/M11216_NANO/SUMMARY/bin_by_bin/M11216_NANO_MAG_00005/M11216_NANO_MAG_00005-contigs.fa |
...          | ...                                                                                                                        | ...                                                                                           |

And let's save it as "HYBRID_MAGS/reassembly_list.txt".

```bash
while read MAG ILLUMINA_PATH NANOPORE_PATH; do
  mkdir HYBRID_MAGS/$MAG

  if [[ $MAG == KUL* ]]; then
    I_ASSEMBLY=UPLAND_CO
    N_ASSEMBLY=M11216_NANO
    SAMPLE=m11216
  fi

  if [[ $MAG == KWL* ]]; then    
    I_ASSEMBLY=FEN_CO
    N_ASSEMBLY=M12208_NANO
    SAMPLE=m12208
  fi

  # Make list of contigs from the Illumina MAG
  grep '>' $ILLUMINA_PATH | sed 's/>//' > HYBRID_MAGS/$MAG/illumina_contigs.txt

  # Make list of contigs from the Nanopore MAG
  grep '>' $NANOPORE_PATH | sed 's/>//' > HYBRID_MAGS/$MAG/nanopore_contigs.txt

  # Make list of Illumina reads mapping to Illumina MAG
  samtools view BINNING/$I_ASSEMBLY/MAPPING/$SAMPLE.bam | grep -F -f HYBRID_MAGS/$MAG/illumina_contigs.txt | cut -f 1 > HYBRID_MAGS/$MAG/reads_illumina_MAG.txt

  # Make list of Illumina reads mapping to the Nanopore MAG
  samtools view BINNING/$N_ASSEMBLY/MAPPING/$SAMPLE.bam | grep -F -f HYBRID_MAGS/$MAG/nanopore_contigs.txt | cut -f 1 > HYBRID_MAGS/$MAG/reads_nanopore_MAG.txt

  # Get mapped Illumina reads
  cat HYBRID_MAGS/$MAG/reads_nanopore_MAG.txt HYBRID_MAGS/$MAG/reads_illumina_MAG.txt | sort | uniq > HYBRID_MAGS/$MAG/illumina_reads.txt

  seqtk subseq TRIMMED_ILLUMINA/$SAMPLE.R1.fastq HYBRID_MAGS/$MAG/illumina_reads.txt > HYBRID_MAGS/$MAG/illumina_reads_R1.fastq
  seqtk subseq TRIMMED_ILLUMINA/$SAMPLE.R2.fastq HYBRID_MAGS/$MAG/illumina_reads.txt > HYBRID_MAGS/$MAG/illumina_reads_R2.fastq

  # Assemble with Unicycler
  unicycler --short1 HYBRID_MAGS/$MAG/illumina_reads_R1.fastq \
            --short2 HYBRID_MAGS/$MAG/illumina_reads_R2.fastq \
            --long $NANOPORE_PATH \
            --out HYBRID_MAGS/$MAG/UNICYCLER \
            --threads $NTHREADS \
            --mode bold

  # Refine bins with anvi'o
  anvi-script-reformat-fasta HYBRID_MAGS/$MAG/UNICYCLER/assembly.fasta \
                             --output-file HYBRID_MAGS/$MAG/assembly_clean.fasta \
                             --prefix $MAG \
                             --simplify-names

  anvi-gen-contigs-database --contigs-fasta HYBRID_MAGS/$MAG/assembly_clean.fasta \
                            --output-db-path HYBRID_MAGS/$MAG/CONTIGS.db \
                            --project-name $MAG

  anvi-run-hmms --contigs-db HYBRID_MAGS/$MAG/CONTIGS.db \
                --num-threads $NTHREADS

  bowtie2-build HYBRID_MAGS/$MAG/assembly_clean.fasta \
                HYBRID_MAGS/$MAG/contigs

  bowtie2 -1 POOLED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 POOLED_ILLUMINA/$SAMPLE.R2.fastq \
          -S HYBRID_MAGS/$MAG/$SAMPLE.sam \
          -x HYBRID_MAGS/$MAG/contigs \
          --threads $NTHREADS \
          --no-unal

  samtools view -F 4 -bS HYBRID_MAGS/$MAG/$SAMPLE.sam | samtools sort > HYBRID_MAGS/$MAG/$SAMPLE.bam
  samtools index HYBRID_MAGS/$MAG/$SAMPLE.bam

  anvi-profile --input-file HYBRID_MAGS/$MAG/$SAMPLE.bam \
               --output-dir HYBRID_MAGS/$MAG/PROFILE \
               --contigs-db HYBRID_MAGS/$MAG/CONTIGS.db \
               --min-contig-length 0 \
               --cluster-contigs \
               --num_threads $NTHREADS

  anvi-interactive --profile-db HYBRID_MAGS/$MAG/PROFILE/PROFILE.db \
                   --contigs-db HYBRID_MAGS/$MAG/CONTIGS.db

  anvi-summarize --contigs-db HYBRID_MAGS/$MAG/CONTIGS.db \
                 --pan-or-profile-db HYBRID_MAGS/$MAG/PROFILE/PROFILE.db \
                 --output-dir HYBRID_MAGS/$MAG/SUMMARY \
                 --collection-name DEFAULT
done < HYBRID_MAGS/reassembly_list.txt
```

## Next step

Continue to [working with MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/06-working-with-MAGs.md).
