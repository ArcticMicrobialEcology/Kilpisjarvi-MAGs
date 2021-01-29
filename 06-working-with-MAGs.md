# Working with MAGs

### Dereplicate MAGs with fastANI

We will dereplicate the Illumina and Nanopore MAGs to create a set of non-redundant MAGs.  
Then we will manually go through the hybrid MAGs and replace the respective Nanopore/Illumina MAGs if they are of better quality.

```bash
printf '%s\t%s\n' name path > MAGS/fastANI_paths.txt

# Get Illumina and Nanopore MAGs
for MAG in `ls BINNING/FINAL_MAGS | sed 's/.fa//g'`; do
  FILE=`ls BINNING/FINAL_MAGS/$MAG.fa`

  printf '%s\t%s\n' $MAG $FILE
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
  cp BINNING/FINAL_MAGS/$MAG.fa NR_MAGS
done

# Manually replace by hybrid MAGs
## ILLUMINA MAG   NANOPORE MAG
## KUL_0013       KUL_0178
## KUL_0074       KUL_0182
## KUL_0098       KUL_0187
## KUL_0153       M11216_NANO_Bin_00060
## KWL_0052       M12208_NANO_Bin_00060
## KWL_0132       KWL_0361
## KWL_0136       KWL_0344
## KWL_0143       M12208_NANO_Bin_00062
## KWL_0150       KWL_0338
## KWL_0153       M12208_NANO_Bin_00116
## KWL_0206       M12208_NANO_Bin_00050
## KWL_0209       KWL_0328
## KWL_0212       M12208_NANO_Bin_00095
## KWL_0315       KWL_0352
```

### Do phylogenomic analysis with GTDB-Tk

```bash
GTDBTK_DATA_PATH=$HOME/GTDBTK_DATA_R95 # Change here to the location of the GTDB-Tk reference data in your system

gtdbtk classify_wf --genome_dir MAGS/NR_MAGS \
                   --out_dir MAGS/GTDB \
                   --extension fa \
                   --cpus $NTHREADS \
                   --pplacer_cpus $NTHREADS

gtdbtk ani_rep --genome_dir MAGS/NR_MAGS \
               --out_dir MAGS/GTDB/ANI \
               --extension fa \
               --cpus 40

gtdb_to_ncbi_majority_vote.py --gtdbtk_output_dir MAGS/GTDB \
                              --output_file MAGS/GTDB/NCBI_taxonomy.txt \
                              --ar122_metadata_file $GTDBTK_DATA_PATH/ar122_metadata_r95.tsv \
                              --bac120_metadata_file $GTDBTK_DATA_PATH/bac120_metadata_r95.tsv
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

SAMPLES=`cut -f 1 sample_metadata.tsv | sed '1d' | uniq`

for SAMPLE in $SAMPLES; do
  bowtie2 -1 POOLED_ILLUMINA/$SAMPLE.R1.fastq \
          -2 POOLED_ILLUMINA/$SAMPLE.R2.fastq \
          -S MAGS/MAPPING/$SAMPLE.sam \
          -x MAGS/MAPPING/NR_MAGs \
          --threads $NTHREADS \
          --no-unal

  samtools view -F 4 -bS MAGS/MAPPING/$SAMPLE.sam | samtools sort > MAGS/MAPPING/$SAMPLE.bam
  samtools index MAGS/MAPPING/$SAMPLE.bam
  samtools idxstats MAGS/MAPPING/$SAMPLE.bam > MAGS/MAPPING/$SAMPLE.reads.txt

  # Get number of mapped reads
  READS=`samtools view -F 4 MAGS/MAPPING/$SAMPLE.bam | wc -l`
  printf '%s\t%s\n' $SAMPLE $READS >> MAGS/MAPPING/mapped_reads.txt
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

### Get gene calls

```bash
anvi-export-gene-calls --contigs-db MAGS/CONTIGS.db \
                       --output-file MAGS/gene_calls.txt \
                       --gene-caller prodigal \
                       --skip-sequence-reporting

anvi-get-sequences-for-gene-calls --contigs-db MAGS/CONTIGS.db \
                                  --output-file MAGS/gene_calls.fa

anvi-get-sequences-for-gene-calls --contigs-db MAGS/CONTIGS.db \
                                  --output-file MAGS/gene_calls.faa \
                                  --get-aa-sequences
```

### Annotate MAGs against COG with DIAMOND

```bash
anvi-run-ncbi-cogs --contigs-db MAGS/CONTIGS.db \
                   --num-threads $NTHREADS
```

### Annotate MAGs against KEGG with DIAMOND

```bash
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
mkdir MAGS/KOFAM_TMP

KOFAM_DB_DIR=$HOME/KOfam # Change here to where you have downloaded and decompressed the KOfam database

# Run hmmer
while read knum threshold score_type profile_type F_measure nseq nseq_used alen mlen eff_nseq re_pos definition; do
  if [[ $score_type == "full" ]]; then
    hmmsearch -T $threshold --cpu $NTHREADS --tblout MAGS/KOFAM_TMP/$knum.hmm.txt $KOFAM_DB_DIR/profiles/$knum.hmm  MAGS/gene_calls.faa
  else
    hmmsearch --domT $threshold --cpu $NTHREADS --domtblout MAGS/KOFAM_TMP/$knum.hmm.txt $KOFAM_DB_DIR/profiles/$knum.hmm MAGS/gene_calls.faa
  fi
done < <(sed '1d' $KOFAM_DB_DIR/ko_list)

# Concatenate results
while read knum threshold score_type profile_type F_measure nseq nseq_used alen mlen eff_nseq re_pos definition; do
  if [[ $score_type == "full" ]]; then
    sed '/^#/d' MAGS/KOFAM_TMP/$knum.hmm.txt |
    tr -s ' ' '\t' |
    cut -f 1,3,5,6
  else
    sed '/^#/d' MAGS/KOFAM_TMP/$knum.hmm.txt |
    tr -s ' ' '\t' |
    cut -f 1,4,7,8
  fi
done < <(sed '1d' $KOFAM_DB_DIR/ko_list) > MAGS/gene_calls_KOFAM.txt

rm -r MAGS/KOFAM_TMP
```

### R analyses

The code below should be executed within R/Rstudio.

```r
# Load packages
library(tidyverse)


##### IMPORT AND PROCESS DATA #####

# Read metadata
metadata <- read_delim("sample_metadata.tsv", delim = "\t",
                       col_types = cols_only(Sample = col_character(), Site = col_factor(), Layer = col_factor(), Ecosystem = col_factor(levels = c("barren", "heathland", "meadow", "fen")))) %>%
  full_join(read_delim("~/DARKFUNCTIONS/POOLED_ILLUMINA/trimmed_reads.txt", delim = "\t", col_names = c("Sample", "TrimmedReads")), by = "Sample") %>% # Read number of trimmed reads
  full_join(read_delim("~/DARKFUNCTIONS/MAGS/MAPPING/mapped_reads.txt",  delim = "\t", col_names = c("Sample", "MappedReads")),  by = "Sample") %>%    # Read number of mapped reads
  mutate(PctReads = MappedReads / TrimmedReads)

# Create list of samples
SAMPLES <- metadata %>%
  pull(Sample)

# Read GTDB taxonomy
GTDB_tax <- bind_rows(read_delim("MAGS/GTDB/gtdbtk.bac120.summary.tsv", delim = "\t") %>%                        # Read data
                        select(user_genome, classification) %>%
                        rename(MAG = user_genome),
                      read_delim("MAGS/GTDB/gtdbtk.ar122.summary.tsv", delim = "\t") %>%
                        select(user_genome, classification) %>%
                        rename(MAG = user_genome)) %>%
  arrange(MAG) %>%                                                                                               # Reorder MAGs
  mutate(classification = gsub("[a-z]__", "", classification)) %>%                                               # Fix taxonomy
  separate(classification, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% # Split taxonomy
  filter(! MAG %in% c("KUL_0016", "KUL_0018"))                                                                   # Remove two eukaryotes MAGs

# Read GTDB trees
GTDB_tree <- list(bac_tree = ape::read.tree("MAGS/GTDB/gtdbtk.bac120.classify.tree"),
                  arc_tree = ape::read.tree("GTDB/gtdbtk.ar122.classify.tree"))

# Concatenate GTDB trees
GTDB_tree <- ape::bind.tree(GTDB_tree[["bac_tree"]], GTDB_tree[["arc_tree"]])

# Create list of MAGs
MAGs <- GTDB_tax %>%
  pull(MAG)

# Read anvio summary
anvio <- read_delim("MAGS/SUMMARY/bins_summary.txt", delim = "\t") %>%
  rename(MAG = bins) %>%
  filter(MAG %in% MAGs) %>%
  arrange(MAG)

# Read gene calls
gene_calls <- read_delim("MAGS/gene_calls.txt", delim = "\t") %>%
  select(-source, -version) %>%
  mutate(MAG = contig %>% str_extract("(KUL|KWL)+_[0-9]+")) %>%
  filter(MAG %in% MAGs)

# Read MAG detection
detection <- read_delim("MAGS/SUMMARY/bins_across_samples/detection.txt", delim = "\t") %>% # Read data
  rename_all(~str_to_lower(.)) %>%                                                          # Fix sample names
  rename(MAG = bins) %>%                                               
  filter(MAG %in% MAGs) %>%                                                                 # Remove eukaryotes
  arrange(MAG) %>%                                                                          # Reorder MAGs
  select(-all_of(SAMPLES), all_of(SAMPLES)) %>%                                             # Reorder columns
  mutate_at(SAMPLES, function(x) ifelse(x >= 0.5, 1, 0))                                    # Transform to presence/absence (detection > 0.5)

# Summarise detection as percentage of samples in each ecosystem
detection_eco <- detection %>%
  gather("Sample", "Detection", -MAG) %>%
  left_join(metadata, by = "Sample") %>%
  group_by(MAG, Ecosystem) %>%
  summarise(Detection = sum(Detection)) %>%
  ungroup %>%
  spread(Ecosystem, Detection) %>%
  mutate(barren = barren / 4) %>%
  mutate(heathland = heathland / 32) %>%
  mutate(meadow = meadow / 18) %>%
  mutate(fen = fen / 15) %>%
  gather("Ecosystem", "Detection", -MAG) %>%
  spread(Ecosystem, Detection)

# Consider a MAG detected in an ecosystem if it is detected in > 20% of the samples
detection_eco_pa <- detection_eco %>%
  gather("Ecosystem", "Detection", - MAG) %>%
  mutate(Detection = ifelse(Detection < 0.2, 0, 1)) %>%
  spread(Ecosystem, Detection)

# Read MAG coverage
coverage <- lapply(SAMPLES, function (SAMPLE) {
  read_delim(paste("MAGS/MAPPING/", SAMPLE, ".reads.txt", sep = ""), delim = "\t", col_names = c("contig", "length", "mapped", "unmapped")) %>%
    mutate(MAG = contig %>% str_extract("(KUL|KWL)+_[0-9]+")) %>%
    mutate(Sample = SAMPLE) %>%
    filter(MAG %in% MAGs) %>%
    select(MAG, Sample, mapped)
}) %>%
  bind_rows %>%
  group_by(MAG, Sample) %>%
  summarise(abundance = sum(mapped)) %>%
  ungroup %>%
  spread(Sample, abundance)

# Transform coverage to relative abundance (relative to number of trimmed reads)
coverage <- coverage %>%
  gather("Sample", "abundance", -MAG) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(abundance_rel = abundance / TrimmedReads) %>%
  select(MAG, Sample, abundance_rel) %>%
  spread(Sample, abundance_rel) %>%
  select(-all_of(SAMPLES), all_of(SAMPLES))

# Summarise means and standard deviations
summariseMeans <- function(x) {
  SAMPLES <- intersect(SAMPLES, x %>% names)

  DATA <- x %>%
    select(all_of(SAMPLES))

  bind_cols(Mean = apply(DATA, 1, mean),
            SD = apply(DATA, 1, sd))
}

coverage.sum <- list(All = coverage %>%
                       summariseMeans,
                     Barren = coverage %>%
                       select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>%
                       summariseMeans,
                     Heathland = coverage %>%
                       select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>%
                       summariseMeans,
                     Meadow = coverage %>%
                       select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "meadow") %>% pull(Sample)) %>%
                       summariseMeans,
                     Fen = coverage %>%
                       select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "fen") %>% pull(Sample)) %>%
                       summariseMeans)

# Load KEGG auxiliary files
keggR::loadKEGG("~/KEGG") # Change here to the location of the KEGG data in your system

# Read KEGG annotation
kegg <- keggR::readBlast("MAGS/gene_calls_KEGG.txt", evalue = 0.00001)

## Read KOFAM annotation
kofam <- keggR::readCustom("MAGS/gene_calls_KOFAM.txt", colnames = c("sequence", "KO", "evalue", "bitscore"))

## Assign KOs
kegg <- kegg %>%
  keggR::assignKEGG()

kofam <- kofam %>%
  keggR::assignKEGG()

## Split KO table by MAG
kegg <- lapply(MAGs, function(x) {
  GENE_CALLS <- gene_calls %>%
    filter(MAG == x) %>%
    pull(gene_callers_id)

  kegg %>%
    keggR::filterKOtable(GENE_CALLS)
}) %>% set_names(MAGs)

kofam <- lapply(MAGs, function(x) {
  GENE_CALLS <- gene_calls %>%
    filter(MAG == x) %>%
    pull(gene_callers_id)

  kofam %>%
    keggR::filterKOtable(GENE_CALLS)
}) %>% set_names(MAGs)


##### DETECTION ANALYSES #####

# Plot heatmap
plotHeatmap <- function(x, annotation_colors, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem, Layer) %>%
    pull(Sample)

  ANNOTATION_COL <- data.frame(Layer = metadata %>% select(Layer),
                               Ecosystem = metadata %>% select(Ecosystem),
                               row.names = metadata %>% pull(Sample))

  ANNOTATION_COLORS <- list(Ecosystem = c(barren = "#caa1f7", heathland = "#61b7d9", meadow = "#63d6aa", fen = "#f9b99f"),
                            Layer = c(mineral = "#b7d8ff", organic = "#98c699"))

  DATA <- data.frame(x %>%
                       select(all_of(ORDER)) %>%
                       sqrt, row.names = rownames(x))

  pheatmap::pheatmap(DATA, color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
                     border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
                     annotation_col = ANNOTATION_COL, annotation_colors = ANNOTATION_COLORS, ...)
}

## All MAGs
detection %>%
  column_to_rownames("MAG") %>%
  dist("euclidean") %>%                                                                            # Cluster MAGs by detection
  hclust("average") %>%
  pluck("order") %>%
  slice(detection, .) %>%
  select(-all_of(SAMPLES), metadata %>% arrange(Ecosystem, Layer) %>% pull(Sample)) %>%            # Reorder samples by ecosystem
  left_join(GTDB_tax %>% select(MAG, Domain, Phylum, Class), by = "MAG") %>%                       # Get taxonomy
  mutate(Phylum = ifelse(Phylum == "Proteobacteria", paste(Phylum, Class, sep = "_"), Phylum)) %>% # For Proteobacteria keep class level
  select(-Class) %>%
  arrange(Domain, Phylum) %>%                                                                      # Arrange by phylum
  plotHeatmap(annotation_row = data.frame(select(., Phylum), row.names = rownames(.)),
              labels_row = NA)

## Only MAGs detected across all ecosystems
detection_eco %>%
  filter(heathland & meadow & fen > 0) %>%
  select(MAG) %>%
  left_join(detection, by = "MAG") %>%
  left_join(GTDB_tax, by = "MAG") %>%
  unite(Taxonomy, Domain, Phylum, Class, Order, Family, Genus, Species, MAG) %>%
  arrange(Taxonomy) %>%
  plotHeatmap(labels_row = pull(., Taxonomy))

# Plot Venn diagram
VennDiagram::venn.diagram(x = list(detection_eco %>%
                                     filter(barren > 0) %>%
                                     pull(MAG),
                                   detection_eco %>%
                                     filter(fen > 0) %>%
                                     pull(MAG),
                                   detection_eco %>%
                                     filter(heathland > 0) %>%
                                     pull(MAG),
                                   detection_eco %>%
                                     filter(meadow > 0) %>%
                                     pull(MAG)),
                          category.names = c("barren" , "fen", "heathland", "meadow"),
                          filename = NULL,
                          imagetype = "png",
                          output = TRUE,
                          lwd = 2,
                          lty = 'blank',
                          fill = c("#caa1f7", "#f9b99f", "#61b7d9", "#63d6aa")) %>%
  grid::grid.draw()


##### PLOT ABUNDANCE BARPLOT #####

# Ten most abundant MAGs in each ecosystem
bind_rows(bind_cols(coverage, coverage.sum[["Barren"]]) %>%
            arrange(desc(Mean)) %>%
            slice(1:10),
          bind_cols(coverage, coverage.sum[["Heathland"]]) %>%
            arrange(desc(Mean)) %>%
            slice(1:10),
          bind_cols(coverage, coverage.sum[["Meadow"]]) %>%
            arrange(desc(Mean)) %>%
            slice(1:10),
          bind_cols(coverage, coverage.sum[["Fen"]]) %>%
            arrange(desc(Mean)) %>%
            slice(1:10)) %>%
  select(-Mean, -SD) %>%
  unique %>%
  left_join(GTDB_tax, by = "MAG") %>%
  unite(Taxonomy, Domain, Phylum, Class, Order, Family, Genus, Species, MAG) %>%
  gather("Sample", "Abundance", -Taxonomy) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Eco = Ecosystem %>% fct_recode("1" = "barren", "2" = "heathland", "3" = "meadow", "4" = "fen") %>% as.numeric) %>%
  mutate(Sample = reorder(Sample, Eco)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Ecosystem)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(Taxonomy), scales = "free", ncol = 4) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") +
  scale_fill_manual(name = "Ecosystem", values = c("barren" = "#caa1f7", "heathland" = "#61b7d9", "meadow" = "#63d6aa", "fen" = "#f9b99f")) +
  scale_y_continuous(labels = scales::percent)


##### SAVE R ENVIRONMENT #####
save.image("MAGS/MAGs.RData")
```

## Next step

Continue to [denitrifier MAGs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/07-denitrifier-MAGs.md).
