# Denitrifier MAGs 

### Find denitrifier MAGs in R

We will use R to go through the annotations and find MAGs encoding denitrifying reductases.  
We will export the gene calls for the key enzymes, which will then be checked against other databases.  
The code below should be executed within R/Rstudio.

```r
# Load packages
library(tidyverse)


##### LOAD PREVIOUS DATA #####

load("MAGS/MAGs.RData")


##### NITRATE REDUCERS #####

dir.create("DENITRIFIERS/NXRA_NARG")

# List genes involved in the pathway
nar <- list(KO = c("nxrA,narG" = "K00370", "nxrB,narH" = "K00371", "napA" = "K02567", "nrfA" = "K03385"))

# Find MAGs with key genes (narG or napA)
nar[["MAGs"]] <- lapply(MAGs, function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    select(KO, gene) %>%
    mutate(MAG = MAG) %>%
    unique
}) %>%
  bind_rows %>%
  filter(KO %in% c("K00370", "K02567")) %>%
  pull(MAG) %>%
  unique

# Get nxrA/narG gene calls
nar[["GENE_CALLS"]] <- lapply(nar[["MAGs"]], function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    filter(KO == "K00370") %>%
    pull(sequence) %>%
    unique %>%
    sort
}) %>%
  unlist

nar[["GENE_CALLS"]] %>%
  write_lines("DENITRIFIERS/NXRA_NARG/nxrA_narG_gene_calls.txt")


##### NITRITE REDUCERS #####

dir.create("DENITRIFIERS/NIRKS")

# List genes involved in the pathway
nir <- list(KO = c("nirK" = "K00368", "nirS" = "K15864"))

# Find MAGs with key genes (one of nirKS)
nir[["MAGs"]] <- lapply(MAGs, function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    select(KO, gene) %>%
    mutate(MAG = MAG) %>%
    unique
}) %>%
  bind_rows %>%
  filter(KO %in% c("K00368", "K15864")) %>%
  pull(MAG) %>%
  unique

# Get nirKS gene calls
nir[["GENE_CALLS"]] <- lapply(nir[["MAGs"]], function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    filter(KO %in% c("K00368", "K15864")) %>%
    pull(sequence) %>%
    unique %>%
    sort
}) %>%
  unlist

nir[["GENE_CALLS"]] %>%
  write_lines("DENITRIFIERS/NIRKS/nirKS_gene_calls.txt")


##### NITRIC OXIDE REDUCERS #####

dir.create("DENITRIFIERS/NORB")

# List genes involved in the pathway
nor <- list(KO = c("norB" = "K04561", "norC" = "K02305"))

# Find MAGs with key genes (norB)
nor[["MAGs"]] <- lapply(MAGs, function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    select(KO, gene) %>%
    mutate(MAG = MAG) %>%
    unique
}) %>%
  bind_rows %>%
  filter(KO == "K04561") %>%
  pull(MAG) %>%
  unique

# Get norB gene calls
nor[["GENE_CALLS"]] <- lapply(nor[["MAGs"]], function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    filter(KO == "K04561") %>%
    pull(sequence) %>%
    unique %>%
    sort
}) %>%
  unlist

nor[["GENE_CALLS"]] %>%
  write_lines("DENITRIFIERS/ORB/norB_gene_calls.txt")


##### NITROUS OXIDE REDUCERS #####

dir.create("DENITRIFIERS/NOSZ")

# List genes involved in the pathway
nos <- list(KO = c("nosZ" = "K00376"))

# Find MAGs with key genes (nosZ)
nos[["MAGs"]] <- lapply(MAGs, function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    select(KO, gene) %>%
    mutate(MAG = MAG) %>%
    unique
}) %>%
  bind_rows %>%
  filter(KO == "K00376") %>%
  pull(MAG) %>%
  unique

# Get nosZ gene calls
nos[["GENE_CALLS"]] <- lapply(nos[["MAGs"]], function(MAG) {
  bind_rows(kegg[[MAG]] %>% keggR::getKOtable(),
            kofam[[MAG]] %>% keggR::getKOtable()) %>%
    filter(KO == "K00376") %>%
    pull(sequence) %>%
    unique %>%
    sort
}) %>%
  unlist

nos[["GENE_CALLS"]] %>%
  write_lines("NOSZ/nosZ_gene_calls.txt")


##### SAVE R ENVIRONMENT #####
save.image("DENITRIFIERS/denitrifiers.RData")
```

### Annotate key genes against KEGG, RefSeq and SWISS

```bash
for GENE in NXRA_NARG/nxrA_narG NIRKS/nirKS NORB/norB NOSZ/nosZ; do
  # Get sequences
  seqtk subseq MAGS/gene_calls.faa DENITRIFIERS/"$GENE"_gene_calls.txt > DENITRIFIERS/$GENE.faa

  # DIAMOND against KEGG
  KEGG_DB_DIR=$HOME/kegg/genes # Change here to the location of the PROKARYOTES.pep.gz file in your system

  diamond blastp --query DENITRIFIERS/$GENE.faa  \
                 --out DENITRIFIERS/"$GENE"_blast_KEGG.txt \
                 --db $KEGG_DB_DIR/PROKARYOTES \
                 --outfmt 6 qseqid sseqid stitle pident qcovhsp evalue bitscore \
                 --threads $NTHREADS

  # BLAST against RefSeq
  blastp -query $GENE.faa \
         -out "$GENE"_blast_refseq.txt \
         -db refseq_protein \
         -outfmt "6 qseqid sseqid stitle pident qcovs evalue bitscore" \
         -num_threads $NTHREADS

  # BLAST against SWISS
  blastp -query $GENE.faa \
         -out "$GENE"_blast_swiss.txt \
         -db swiss \
         -outfmt "6 qseqid sseqid stitle pident qcovs evalue bitscore" \
         -num_threads $NTHREADS
done
```

### Parse BLAST results in R

Now we will import the KEGG, RefSeq and SWISS annotations to R and check manually the results.  
MAGs encoding denitrifying reductases with weak annotations will be removed from downstream analyses.  
The code below should be executed within R/Rstudio.

```r
# Load packages
library(tidyverse)


##### LOAD PREVIOUS DATA #####

load("DENITRIFIERS/denitrifiers.RData")


##### NITRATE REDUCERS #####

# Read BLAST results
nar[["KEGG"]] <- read_delim("DENITRIFIERS/NXRA_NARG/nxrA_narG_blast_KEGG.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  separate(stitle, into = c(NA, "stitle"), sep = "  ") %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  left_join(.PROKARYOTES.DAT, by = c(sseqid = "gene")) %>%
  left_join(.KO00000, by = "KO") %>%
  mutate(taxcode = separate(., sseqid, into = c("taxcode", NA), sep = ":") %>% pull(taxcode)) %>%
  left_join(.TAXONOMIC_RANK, by = "taxcode") %>%
  select(-taxcode, -domain, -phylum, -genus, -species)

nar[["REFSEQ"]] <- read_delim("DENITRIFIERS/NXRA_NARG/nxrA_narG_blast_refseq.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "\\[.+\\]")) %>%
  mutate(stitle = str_remove(stitle, " \\[.+\\]"))

nar[["SWISS"]] <- read_delim("DENITRIFIERS/NXRA_NARG/nxrA_narG_blast_swiss.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "OS=.*OX") %>% str_remove("OS=") %>% str_remove(" OX")) %>%
  mutate(gene = str_extract(stitle, "GN=.*PE") %>% str_remove("GN=") %>% str_remove(" PE")) %>%
  mutate(stitle = str_remove(stitle, " OS=.*"))

# Remove MAGs with weak annotations
nar[["MAGs"]] <- setdiff(nar[["MAGs"]], c("KWL_0111", "KWL_0223", "KWL_0351", "KWL_0150", "KWL_0020", "KWL_0185", "KWL_0134")) %>%
  sort


##### NITRITE REDUCERS #####

# Read BLAST results
nir[["KEGG"]] <- read_delim("DENITRIFIERS/NIRKS/nirKS_blast_KEGG.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  separate(stitle, into = c(NA, "stitle"), sep = "  ") %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  left_join(.PROKARYOTES.DAT, by = c(sseqid = "gene")) %>%
  left_join(.KO00000, by = "KO") %>%
  mutate(taxcode = separate(., sseqid, into = c("taxcode", NA), sep = ":") %>% pull(taxcode)) %>%
  left_join(.TAXONOMIC_RANK, by = "taxcode") %>%
  select(-taxcode, -domain, -phylum, -genus, -species)

nir[["REFSEQ"]] <- read_delim("DENITRIFIERS/NIRKS/nirKS_blast_refseq.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "\\[.+\\]")) %>%
  mutate(stitle = str_remove(stitle, " \\[.+\\]"))

nir[["SWISS"]] <- read_delim("DENITRIFIERS/NIRKS/nirKS_blast_swiss.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "OS=.*OX") %>% str_remove("OS=") %>% str_remove(" OX")) %>%
  mutate(gene = str_extract(stitle, "GN=.*PE") %>% str_remove("GN=") %>% str_remove(" PE")) %>%
  mutate(stitle = str_remove(stitle, " OS=.*"))

# Remove MAGs with weak annotations
nir[["MAGs"]] <- setdiff(nir[["MAGs"]], c("KWL_0109", "KWL_0148", "KWL_0351", "KWL_0002", "KWL_0138", "KWL_0357", "KWL_0028", "KWL_0004", "KWL_0027", "KWL_0283", "KWL_0061", "KWL_0062", "KWL_0213")) %>%
  sort


##### NITRIC OXIDE REDUCERS #####

# Read BLAST results
nor[["KEGG"]] <- read_delim("DENITRIFIERS/NORB/norB_blast_KEGG.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  separate(stitle, into = c(NA, "stitle"), sep = "  ") %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  left_join(.PROKARYOTES.DAT, by = c(sseqid = "gene")) %>%
  left_join(.KO00000, by = "KO") %>%
  mutate(taxcode = separate(., sseqid, into = c("taxcode", NA), sep = ":") %>% pull(taxcode)) %>%
  left_join(.TAXONOMIC_RANK, by = "taxcode") %>%
  select(-taxcode, -domain, -phylum, -genus, -species)

nor[["REFSEQ"]] <- read_delim("DENITRIFIERS/NORB/norB_blast_refseq.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "\\[.+\\]")) %>%
  mutate(stitle = str_remove(stitle, " \\[.+\\]"))

nor[["SWISS"]] <- read_delim("DENITRIFIERS/NORB/norB_blast_swiss.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "OS=.*OX") %>% str_remove("OS=") %>% str_remove(" OX")) %>%
  mutate(gene = str_extract(stitle, "GN=.*PE") %>% str_remove("GN=") %>% str_remove(" PE")) %>%
  mutate(stitle = str_remove(stitle, " OS=.*"))

# Remove MAGs with weak annotations
nor[["MAGs"]] <- setdiff(nor[["MAGs"]], c("KWL_0036", "KWL_0054", "KUL_0137", "KWL_0233", "KWL_0157", "KWL_0258", "KWL_0029")) %>%
  sort


##### NITROUS OXIDE REDUCERS #####

# Read BLAST results
nos[["KEGG"]] <- read_delim("DENITRIFIERS/NOSZ/nosZ_blast_KEGG.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  separate(stitle, into = c(NA, "stitle"), sep = "  ") %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  left_join(.PROKARYOTES.DAT, by = c(sseqid = "gene")) %>%
  left_join(.KO00000, by = "KO") %>%
  mutate(taxcode = separate(., sseqid, into = c("taxcode", NA), sep = ":") %>% pull(taxcode)) %>%
  left_join(.TAXONOMIC_RANK, by = "taxcode") %>%
  select(-taxcode, -domain, -phylum, -genus, -species)

nos[["REFSEQ"]] <- read_delim("DENITRIFIERS/NOSZ/nosZ_blast_refseq.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "\\[.+\\]")) %>%
  mutate(stitle = str_remove(stitle, " \\[.+\\]"))

nos[["SWISS"]] <- read_delim("DENITRIFIERS/NOSZ/nosZ_blast_swiss.txt", delim = "\t", col_names = c("gene_callers_id", "sseqid", "stitle", "pident", "qcov", "evalue", "bitscore")) %>%
  left_join(gene_calls %>% select(gene_callers_id, MAG), by = "gene_callers_id") %>%
  mutate(strain = str_extract(stitle, "OS=.*OX") %>% str_remove("OS=") %>% str_remove(" OX")) %>%
  mutate(gene = str_extract(stitle, "GN=.*PE") %>% str_remove("GN=") %>% str_remove(" PE")) %>%
  mutate(stitle = str_remove(stitle, " OS=.*"))

# Remove MAGs with weak annotations
nos[["MAGs"]] <- setdiff(nos[["MAGs"]], c("KWL_0276", "KWL_0370", "KUL_0133")) %>%
  sort


##### SUMMARISE DENITRIFIERS #####

denitrifiers <- tibble(MAG = c(nar[["MAGs"]], nir[["MAGs"]], nor[["MAGs"]], nos[["MAGs"]]) %>% unique) %>%
  left_join(tibble(MAG = nar[["MAGs"]], nar = 1), by = "MAG") %>%
  left_join(tibble(MAG = nir[["MAGs"]], nir = 1), by = "MAG") %>%
  left_join(tibble(MAG = nor[["MAGs"]], nor = 1), by = "MAG") %>%
  left_join(tibble(MAG = nos[["MAGs"]], nos = 1), by = "MAG") %>%
  replace(is.na(.), 0) %>%
  arrange(MAG)


##### CHECK ANNOTATIONS #####

checkAnnotation <- function(x) {
  # Check alignment coverage
  COV <- bind_rows(x[["KEGG"]] %>% mutate(source = "KEGG"),
                   x[["SWISS"]] %>% mutate(source = "SWISS"),
                   x[["REFSEQ"]] %>% mutate(source = "REFSEQ")) %>%
    filter(MAG %in% x[["MAGs"]]) %>%
    group_by(MAG, source) %>%
    arrange(desc(bitscore)) %>%
    slice(1) %>%
    group_by(source) %>%
    summarise(Min = min(qcov), Max = max(qcov), Mean = mean(qcov), SD = sd(qcov))

  # Check identity
  ID <- bind_rows(x[["KEGG"]] %>% mutate(source = "KEGG"),
                  x[["SWISS"]] %>% mutate(source = "SWISS"),
                  x[["REFSEQ"]] %>% mutate(source = "REFSEQ")) %>%
    filter(MAG %in% x[["MAGs"]]) %>%
    group_by(MAG, source) %>%
    arrange(desc(bitscore)) %>%
    slice(1) %>%
    group_by(source) %>%
    summarise(Min = min(pident), Max = max(pident), Mean = mean(pident), SD = sd(pident))

  # Check gene
  GENE <- bind_rows(x[["KEGG"]] %>% mutate(source = "KEGG"),
                    x[["SWISS"]] %>% mutate(source = "SWISS"),
                    x[["REFSEQ"]] %>% mutate(source = "REFSEQ")) %>%
    filter(MAG %in% x[["MAGs"]]) %>%
    group_by(MAG, source) %>%
    arrange(desc(bitscore)) %>%
    slice(1) %>%
    group_by(source, stitle) %>%
    tally %>%
    ungroup %>%
    arrange(source, desc(n))

  # Check strains
  STRAINS_RS <- x[["REFSEQ"]] %>%
    filter(MAG %in% x[["MAGs"]]) %>%
    group_by(MAG) %>%
    arrange(desc(bitscore)) %>%
    slice(1) %>%
    group_by(strain) %>%
    tally %>%
    ungroup %>%
    arrange(desc(n))

  # Return results
  res <- list(COV = COV, ID = ID, GENE = GENE, STRAINS_RS = STRAINS_RS)

  return(res)
}

# nar
checkAnnotation(nar)

# nir
checkAnnotation(nir)

# nor
checkAnnotation(nor)

# nos
checkAnnotation(nos)


##### PLOT TAXONOMIC COMPOSITION IN EACH DENITRIFICATION STEP #####

plotPie <- function(x, y) {
  pie(x$n, labels = x$Phylum, clockwise = T, init.angle = 90, main = y)
}

# Prepare data
pie.chart <- denitrifiers %>%
  left_join(detection_eco_pa, by = "MAG") %>%
  group_by(MAG) %>%
  select(MAG, nar, nir, nor, nos, meadow, fen) %>%
  gather("ecosystem", "detection", -MAG, -nar, -nir, -nor, -nos) %>%
  gather("gene", "presence", -MAG, -ecosystem, -detection) %>%
  filter(detection > 0) %>%
  filter(presence > 0) %>%
  mutate(detection = 1) %>%
  left_join(GTDB_tax, by = "MAG") %>%
  mutate(Phylum = ifelse(Phylum == "Proteobacteria", Class, Phylum))

pie.chart <- bind_rows(pie.chart %>%
                         filter(ecosystem == "meadow") %>%
                         group_by(gene, Phylum, ecosystem) %>%
                         tally,
                       pie.chart %>%
                         filter(ecosystem == "fen") %>%
                         group_by(gene, Phylum, ecosystem) %>%
                         tally)

# Plot
par(mfrow = c(4, 2))

pie.chart %>%
  filter(ecosystem == "meadow") %>%
  filter(gene == "nar") %>%
  arrange(desc(n)) %>%
  plotPie("NAR, meadow")

pie.chart %>%
  filter(ecosystem == "fen") %>%
  filter(gene == "nar") %>%
  arrange(desc(n)) %>%
  plotPie("NAR, fen")

pie.chart %>%
  filter(ecosystem == "meadow") %>%
  filter(gene == "nir") %>%
  arrange(desc(n)) %>%
  plotPie("NIR, meadow")

pie.chart %>%
  filter(ecosystem == "fen") %>%
  filter(gene == "nir") %>%
  arrange(desc(n)) %>%
  plotPie("NIR, fen")

pie.chart %>%
  filter(ecosystem == "meadow") %>%
  filter(gene == "nor") %>%
  arrange(desc(n)) %>%
  plotPie("NOR, meadow")

pie.chart %>%
  filter(ecosystem == "fen") %>%
  filter(gene == "nor") %>%
  arrange(desc(n)) %>%
  plotPie("NOR, fen")

pie.chart %>%
  filter(ecosystem == "meadow") %>%
  filter(gene == "nos") %>%
  arrange(desc(n)) %>%
  plotPie("NOS, meadow")

pie.chart %>%
  filter(ecosystem == "fen") %>%
  filter(gene == "nos") %>%
  arrange(desc(n)) %>%
  plotPie("NOS, fen")


##### PLOT MAG COMPLETION ACCORDING TO THE NUMBER OF DENITRIFICATION GENES #####

anvio %>%
  left_join(denitrifiers, by = "MAG") %>%
  filter(nar | nir | nor | nos > 0) %>%
  select(MAG, percent_completion, percent_redundancy, nar, nir, nor, nos) %>%
  group_by(MAG) %>%
  mutate(ngenes = sum(nar, nir, nor, nos)) %>%
  ungroup %>%
  mutate(type = ifelse(nar  == 1, "Nar", "NA")) %>%
  mutate(type = ifelse(nir  == 1, "Nir", type)) %>%
  mutate(type = ifelse(nor  == 1, "Nor", type)) %>%
  mutate(type = ifelse(nos  == 1, "Nos", type)) %>%
  mutate(type = ifelse(nar & nir == 1, "Nar_Nir", type)) %>%
  mutate(type = ifelse(nar & nor == 1, "Nar_Nor", type)) %>%
  mutate(type = ifelse(nar & nos == 1, "Nar_Nos", type)) %>%
  mutate(type = ifelse(nir & nor == 1, "Nir_Nor", type)) %>%
  mutate(type = ifelse(nir & nos == 1, "Nir_Nos", type)) %>%
  mutate(type = ifelse(nor & nos == 1, "Nor_Nos", type)) %>%
  mutate(type = ifelse(nar & nor & nos == 1, "Nar_Nor_Nos", type)) %>%
  mutate(ngenes = ngenes %>% as_factor) %>%
  ggplot(aes(x = ngenes, y = percent_completion, colour = type)) +
  geom_jitter() +
  facet_grid(cols = vars(ngenes), scales = "free_x")
```
