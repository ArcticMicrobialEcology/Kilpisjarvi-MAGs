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

### R analyses

The code below should be executed within R/Rstudio.

```r
# Load packages
library(tidyverse)
library(keggR)
library(vegan)


##### IMPORT AND PROCESS DATA #####

# Read metadata
metadata <- read_delim("sample_metadata_illumina.txt", delim = "\t") %>%
  select(Sample, Site, Layer, Ecosystem) %>%
  mutate(Ecosystem = as.factor(Ecosystem), Layer = as.factor(Layer)) %>%
  unique

# Create list of samples
SAMPLES <- metadata %>%
  pull(Sample)

# METAXA

## Read METAXA annotation (phylum level)
phylum <- read_delim("METAXA/summary_level_2.txt", delim = "\t") %>%                # Read data
  separate(Taxa, c("Domain", "Phylum"), ";") %>%                                    # Split taxonomy
  filter(Domain %in% c("Bacteria", "Archaea")) %>%                                  # Keep only bacteria and archaea
  select(-all_of(SAMPLES), all_of(SAMPLES))                                         # Reorder samples according to metadata

## Read METAXA annotation (genus level)
genus <- read_delim("METAXA/summary_level_6.txt", delim = "\t") %>%                    # Read data
  separate(Taxa, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), ";") %>%  # Split taxonomy
  filter(Domain %in% c("Bacteria", "Archaea")) %>%                                     # Keep only bacteria and archaea
  select(-all_of(SAMPLES), all_of(SAMPLES))                                            # Reorder samples according to metadata

## Transform to relative abundance
phylum.rel <- phylum %>%
  mutate_at(SAMPLES, function(x) x / sum(x))

genus.rel <- genus %>%
  mutate_at(SAMPLES, function(x) x / sum(x))

# KEGG

## Read KEGG annotation
kegg <- lapply(SAMPLES, function(SAMPLE) {
  readBlast(paste("KEGG/", SAMPLE, ".txt", sep = ""))
}) %>%
  set_names(SAMPLES)

## Assign KO
KOtable <- lapply(SAMPLES, function(SAMPLE) {
  kegg[[SAMPLE]] %>%
    assignKEGG
})

## Merge KO tables
KOtable <- lapply(SAMPLES, function(SAMPLE) {
  KOtable[[SAMPLE]] %>%
    getKOtable %>%
    select(KO, gene) %>%
    mutate(Sample = SAMPLE)
}) %>%
  bind_rows %>%
  group_by(KO, gene, Sample) %>%
  tally %>%
  ungroup %>%
  spread(Sample, n, fill = F)

## Normalise genes based on the abundance of the rpoB gene
KOtable.norm <- KOtable %>%
  select(-gene) %>%
  gather(key = "Sample", value = "Abundance", -KO) %>%
  spread(KO, Abundance) %>%
  mutate(rpoB = K03043) %>%
  gather(key = "KO", value = "Abundance", -Sample, -rpoB) %>%
  mutate(Abundance = Abundance) %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  select(-rpoB) %>%
  spread(Sample, Abundance) %>%
  left_join(.KO00000, by = "KO") %>%
  select(KO, gene, all_of(SAMPLES))

# Summarise means and standard deviations
summariseMeans <- function(x) {
  SAMPLES <- intersect(SAMPLES, x %>% names)

  DATA <- x %>%
    select(all_of(SAMPLES))

  bind_cols(Mean = apply(DATA, 1, mean),
            SD = apply(DATA, 1, sd))
}

phylum.sum <- list(All = phylum.rel %>%
                     summariseMeans,
                   Barren = phylum.rel %>%
                     select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>%
                     summariseMeans,
                   Heathland = phylum.rel %>%
                     select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>%
                     summariseMeans,
                   Meadow = phylum.rel %>%
                     select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "meadow") %>% pull(Sample)) %>%
                     summariseMeans,
                   Fen = phylum.rel %>%
                     select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "fen") %>% pull(Sample)) %>%
                     summariseMeans)

genus.sum <- list(All = genus.rel %>%
                    summariseMeans,
                  Barren = genus.rel %>%
                    select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>%
                    summariseMeans,
                  Heathland = genus.rel %>%
                    select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>%
                    summariseMeans,
                  Meadow = genus.rel %>%
                    select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "meadow") %>% pull(Sample)) %>%
                    summariseMeans,
                  Fen = genus.rel %>%
                    select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "fen") %>% pull(Sample)) %>%
                    summariseMeans)

KOtable.sum <- list(All = KOtable.norm %>%
                      summariseMeans,
                    Barren = KOtable.norm %>%
                      select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>%
                      summariseMeans,
                    Heathland = KOtable.norm %>%
                      select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>%
                      summariseMeans,
                    Meadow = KOtable.norm %>%
                      select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "meadow") %>% pull(Sample)) %>%
                      summariseMeans,
                    Fen = KOtable.norm %>%
                      select(-all_of(SAMPLES), metadata %>% filter(Ecosystem == "fen") %>% pull(Sample)) %>%
                      summariseMeans)


##### PLOT BARPLOTS #####

# METAXA
filterUnclassified <- function(x, y) {
  x %>%
    filter(!!rlang::sym(y) != "") %>%
    filter(!grepl("Unclassified", !!rlang::sym(y))) %>%
    filter(!grepl("Incertae Sedis", !!rlang::sym(y)))
}

## Plot the 5 most abundant genera in each ecosystem
bind_rows(bind_cols(genus.rel, genus.sum[["Barren"]]) %>%
            filterUnclassified("Genus") %>%
            arrange(desc(Mean)) %>% slice(1:5),
          bind_cols(genus.rel, genus.sum[["Heathland"]]) %>%
            filterUnclassified("Genus") %>%
            arrange(desc(Mean)) %>% slice(1:5),
          bind_cols(genus.rel, genus.sum[["Meadow"]]) %>%
            filterUnclassified("Genus") %>%
            arrange(desc(Mean)) %>% slice(1:5),
          bind_cols(genus.rel, genus.sum[["Fen"]]) %>%
            filterUnclassified("Genus") %>%
            arrange(desc(Mean)) %>% slice(1:5)) %>%
  select(-Mean, -SD) %>%
  unique %>%
  gather(key = "Sample", value = "Abundance", -Domain, -Phylum, -Class, -Order, -Family, -Genus) %>%
  mutate(Genus = ifelse(Phylum == "Proteobacteria", paste(Class, Genus, sep = ": "), paste(Phylum, Genus, sep = ": "))) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Eco = Ecosystem %>% fct_recode("1" = "barren", "2" = "heathland", "3" = "meadow", "4" = "fen") %>% as.numeric) %>%
  mutate(Sample = reorder(Sample, Eco)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Ecosystem)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(Genus), scales = "free", ncol = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(name = "Ecosystem", values = c("barren" = "#caa1f7", "heathland" = "#61b7d9", "meadow" = "#63d6aa", "fen" = "#f9b99f")) +
  scale_y_continuous(labels = scales::percent)

# KEGG

## Plot N genes
N.GENES <- c("nifH" = "K02588", "amoA" = "K10944", "nxrA_narG" = "K00370", "napA" = "K02567", "nirB" = "K00362", "nrfA" = "K03385", "nirK" = "K00368", "nirS" = "K15864", "norB" = "K04561", "nosZ" = "K00376")

KOtable.norm %>%
  filter(KO %in% N.GENES) %>%
  gather(key = "Sample", value = "Abundance", -KO, -gene) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Eco = Ecosystem %>% fct_recode("1" = "barren", "2" = "heathland", "3" = "meadow", "4" = "fen") %>% as.numeric) %>%
  mutate(Sample = reorder(Sample, Eco)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Ecosystem)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(gene), scales = "free", ncol = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(name = "Ecosystem", values = c("barren" = "#caa1f7", "heathland" = "#61b7d9", "meadow" = "#63d6aa", "fen" = "#f9b99f")) +
  scale_y_continuous(labels = scales::percent)


##### RICHNESS ANALYSES #####

computeRich <- function(x) {
  tibble(Sample = SAMPLES,
         Value = apply(x, 2, function (y) sum(y > 0))) %>%
    full_join(metadata, by = "Sample")
}

plotBoxplot <- function(DATA) {
  ggplot(DATA, aes(x = Ecosystem, y = Value)) +
    geom_boxplot(aes(fill = Ecosystem), outlier.shape = NA, position = position_dodge(preserve = "single")) +
    geom_point(aes(fill = Ecosystem), shape = 16, position = position_jitterdodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA),
          strip.background = element_rect(fill = "grey90")) +
    scale_fill_manual(name = "Ecosystem", values = c("barren" = "#caa1f7", "heathland" = "#61b7d9", "meadow" = "#63d6aa", "fen" = "#f9b99f"))
}

anovaRich <- function(x) {
  x %>%
    select(all_of(SAMPLES)) %>%
    computeRich %>%
    lm(Value ~ Ecosystem, .) %>%
    aov
}

# METAXA

## Plot boxplot
genus.rel %>%
  filterUnclassified("Genus") %>%
  select(all_of(SAMPLES)) %>%
  computeRich %>%
  plotBoxplot

## ANOVA (main test)
genus.rel %>%
  filterUnclassified("Genus") %>%
  anovaRich %>%
  summary

## ANOVA (pairwise test)
genus.rel %>%
  filterUnclassified("Genus") %>%
  anovaRich %>%
  multcomp::glht(linfct = multcomp::mcp(Ecosystem = "Tukey")) %>%
  summary

# KEGG

## Plot boxplot
KOtable %>%
  select(all_of(SAMPLES)) %>%
  computeRich %>%
  plotBoxplot

## ANOVA (main test)
KOtable %>%
  anovaRich %>%
  summary

## ANOVA (pairwise test)
KOtable %>%
  anovaRich %>%
  multcomp::glht(linfct = multcomp::mcp(Ecosystem = "Tukey")) %>%
  summary


##### MULTIVARIATE ANALYSES #####

plotOrdination <- function(DATA, XLIM, YLIM) {
  attach(metadata, warn.conflicts = F)

  COLORS <- c("#caa1f7", "#61b7d9", "#63d6aa", "#f9b99f")

  plot(DATA, display = "sites", type = "n", xlim = XLIM, ylim = YLIM)

  points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Ecosystem])
  ordiellipse(DATA, groups = Ecosystem, kind = "sd", draw = "polygon", col = COLORS)
  legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
  legend("bottomright", legend = levels(metadata$Ecosystem), bty = "n", col = COLORS, pch =15)

  detach(metadata)
}

plotDensity <- function(DATA, AXIS) {
  DATA %>%
    metaMDS(distance = "bray") %>%
    scores(display = "sites") %>%
    as.data.frame %>%
    rownames_to_column("Sample") %>%
    left_join(metadata, by = "Sample") %>%
    ggplot(aes(x = !!enquo(AXIS), fill = Ecosystem)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position = "none", panel.grid = element_blank()) +
    scale_fill_manual(values = c("#caa1f7", "#61b7d9", "#63d6aa", "#f9b99f"))
}

# METAXA

## Prepare data
genus.ord <- genus.rel %>%
  filterUnclassified("Genus") %>%
  select(Family, Genus, all_of(SAMPLES)) %>%
  mutate(Family_Genus = paste(Family, Genus, sep = ",")) %>%
  column_to_rownames("Family_Genus") %>%
  select(-Family, -Genus) %>%
  t %>%
  sqrt %>%
  as_tibble(rownames = NA)

## PERMANOVA
adonis(genus.ord ~ Ecosystem, metadata, permutations = 9999, method = "bray")
adonis(genus.ord ~ Layer, metadata, permutations = 9999, method = "bray")
adonis(genus.ord ~ Ecosystem*Layer, metadata, permutations = 9999, method = "bray")

## NMDS
genus.ord %>%
  metaMDS(distance = "bray") %>%
  plotOrdination(XLIM = c(-1, 2), YLIM = c(-1, 1.5))

## Density plots
genus.ord %>%
  plotDensity(NMDS1) +
  scale_x_continuous(limits = c(-1, 2), breaks = seq(-1, 2, by = 0.5))

genus.ord %>%
  plotDensity(NMDS2) +
  scale_x_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, by = 0.5))

# KEGG

## Prepare data
KOtable.ord <- KOtable.norm %>%
  select(KO, all_of(SAMPLES)) %>%
  column_to_rownames("KO") %>%
  t %>%
  sqrt %>%
  as_tibble(rownames = NA)

## PERMANOVA
adonis(KOtable.ord ~ Ecosystem, metadata, permutations = 9999, method = "bray")
adonis(KOtable.ord ~ Layer, metadata, permutations = 9999, method = "bray")
adonis(KOtable.ord ~ Ecosystem*Layer, metadata, permutations = 9999, method = "bray")

## NMDS
KOtable.ord %>%
  metaMDS(distance = "bray") %>%
  plotOrdination(XLIM = c(-0.2, 0.5), YLIM = c(-0.2, 0.2))

## Density plots
KOtable.ord %>%
  plotDensity(NMDS1) +
  scale_x_continuous(limits = c(-0.2, 0.5), breaks = seq(-0.2, 0.5, by = 0.1))

KOtable.ord %>%
  plotDensity(NMDS2) +
  scale_x_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, by = 0.1))


##### NITROGEN GENES #####

runAnova <- function(x, y) {
  df <- x %>%
    gather(key = "Sample", value = "Value") %>%
    full_join(metadata, by = "Sample")

  f <- as.formula(paste("Value", y, sep = " ~ "))

  lm <- lm(f, df)

  anova <- lm %>%
    aov %>%
    summary %>%
    unlist

  F.VALUE <- anova["F value1"]
  COEFFICIENT <- lm[["coefficients"]][2] %>% as.vector
  P.VALUE <- anova["Pr(>F)1"]

  res <- tibble(f.value = F.VALUE,
                coefficient = COEFFICIENT,
                p.value = P.VALUE)

  return(res)
}

# ANOVA
lapply(N.GENES, function(x) {
  GENE <- names(N.GENES)[N.GENES == x]

  KOtable %>%
    filter(KO == x) %>%
    select(all_of(SAMPLES)) %>%
    runAnova("Ecosystem") %>%
    mutate(KO = x, gene = GENE)
}) %>%
  bind_rows %>%
  mutate(p.adjust = p.adjust(p.value, method = "fdr"))

```


## Next step

Continue to [metagenome assembling](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/03-assembling.md).
