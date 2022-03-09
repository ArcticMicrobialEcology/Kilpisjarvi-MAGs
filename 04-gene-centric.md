# Gene-centric analyses

```bash
mkdir GENE_CENTRIC
```

### Map reads to contigs with CoverM

```bash
mkdir GENE_CENTRIC/COVERM

for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
  if [[ $ASSEMBLY == 'UPLAND_CO' ]]; then
    SAMPLES=`awk -F '\t' '{if ($4 == "upland") {print $1}}' sample_metadata.tsv | uniq`
  fi

  if [[ $ASSEMBLY == 'FEN_CO' ]]; then
    SAMPLES=`awk -F '\t' '{if ($4 == "fen") {print $1}}' sample_metadata.tsv | uniq`
  fi

  R1=`printf '%s\n' $SAMPLES | awk -v ORS=" " '{print "POOLED_ILLUMINA/" $0 ".R1.fastq"}' | sed 's/$/\n/'`
  R2=`printf '%s\n' $SAMPLES | awk -v ORS=" " '{print "POOLED_ILLUMINA/" $0 ".R2.fastq"}' | sed 's/$/\n/'`

  coverm contig -1 $R1 \
                -2 $R2 \
                --reference BINNING/$ASSEMBLY/CONTIGS_2500nt.fa \
                --output-file GENE_CENTRIC/COVERM/$ASSEMBLY.txt \
                --methods rpkm \
                --min-covered-fraction 0 \
                --threads $NTHREADS
done
```

### Export amino acid sequences

```bash
mkdir GENE_CENTRIC/GENE_CALLS

for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
  anvi-get-sequences-for-gene-calls --contigs-db BINNING/$ASSEMBLY/CONTIGS.db \
                                    --output-file GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.faa \
                                    --get-aa-sequences
done
```

### Get map of gene calls to splits and splits to bins

```bash
for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
  # Gene calls to splits
  sqlite3 BINNING/$ASSEMBLY/CONTIGS.db 'SELECT * FROM genes_in_splits' | tr "|" "\t" | awk -F "\t" -v OFS="\t" '{print $2, $3}' | sort -k 2,2 >  GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.gene_calls_to_splits.txt

  # Splits to bins
  if [[ $ASSEMBLY == 'UPLAND_CO' || $ASSEMBLY == 'FEN_CO' ]]; then
    for CLUSTER in `ls BINNING/$ASSEMBLY/CONCOCT_SPLIT | grep -v SUMMARY`; do
      sqlite3 BINNING/$ASSEMBLY/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT * FROM collections_of_splits' | grep 'FINAL' | tr "|" "\t" | awk -F "\t" -v OFS="\t" -v CLUSTER=$CLUSTER '{print CLUSTER, $3, $4}'
    done | sort -k 1,1 > GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.splits_to_bins.txt
  fi

  if [[ $ASSEMBLY == 'M11216_NANO' || $ASSEMBLY == 'M12208_NANO' ]]; then
    sqlite3 BINNING/$ASSEMBLY/MERGED_PROFILES/PROFILE.db 'SELECT * FROM collections_of_splits' | grep 'FINAL' | tr "|" "\t" | awk -F "\t" -v OFS="\t" '{print $3, $4}' | sort -k 1,1 > GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.splits_to_bins.txt
  fi
done
```

### Annotate genes with KofamScan

```bash
# Download KOfam database
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
gunzip ko_list.gz

wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
tar zxf profiles.tar.gz

# Run KofamScan
mkdir GENE_CENTRIC/KOFAM_SCAN

KO_LIST=`cut -f 1 ko_list | sed '1d'`

for KO in $KO_LIST do;
  mkdir GENE_CENTRIC/KOFAM_SCAN/$KO

  for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
    exec_annotation GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.faa \
                    --profile profiles/$KO.hmm \
                    --ko-list ko_list \
                    --format detail-tsv \
                    --cpu $NTHREADS | grep ^* > GENE_CENTRIC/KOFAM_SCAN/$KO/$ASSEMBLY.txt
  done 
done
```

### Analyses of denitrification genes

Here we will analyse further the denitrification genes:

- nirK: K00368
- nirS: K15864
- norB: K04561
- nosZ: K00376

```bash
# Get hits
for KO in K00368 K15864 K04561 K00376; do
  for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
    cut -f 2 GENE_CENTRIC/KOFAM_SCAN/$KO/$ASSEMBLY.txt | sort | uniq > GENE_CENTRIC/KOFAM_SCAN/$KO/$ASSEMBLY.gene_calls.txt
  done

  for ASSEMBLY in UPLAND_CO FEN_CO M11216_NANO M12208_NANO; do
    seqtk subseq GENE_CENTRIC/GENE_CALLS/$ASSEMBLY.faa GENE_CENTRIC/KOFAM_SCAN/$KO/$ASSEMBLY.gene_calls.txt | sed -e "s/>/>$ASSEMBLY./"
  done > GENE_CENTRIC/KOFAM_SCAN/$KO/assemblies.faa
done

# Align with MAFFT
for KO in K00368 K15864 K04561 K00376; do
  mafft --auto \
        --reorder \
        --thread $NTHREADS \
        GENE_CENTRIC/KOFAM_SCAN/$KO/assemblies.faa > GENE_CENTRIC/KOFAM_SCAN/$KO/assemblies.aln.faa
done
```

We will then check the alignments manually for the presence of conserved residues at key positions.  
Sequences that do not contain the correct amino acid will be removed.

## Next step

Continue to [working with MAGs](05-working-with-MAGs.md).