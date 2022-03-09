# Metagenome assembling

## Illumina data

Based on the taxonomic profiles, samples appear too different for a good co-assembly.  
We will then do two co-assemblies, one with the fen samples only and the other with the remaining (upland) samples.

### Define assembly and create list of file names

```bash
# Assembly of upland samples
ASSEMBLY=UPLAND_CO
SAMPLES=`awk -F '\t' '{if ($4 == "upland") {print $1}}' sample_metadata.tsv | uniq`
R1=`printf '%s\n' $SAMPLES | awk -v ORS="," '{print "POOLED_ILLUMINA/" $0 ".R1.fastq"}' | sed 's/,$/\n/'`
R2=`printf '%s\n' $SAMPLES | awk -v ORS="," '{print "POOLED_ILLUMINA/" $0 ".R2.fastq"}' | sed 's/,$/\n/'`

# Assembly of fen samples
ASSEMBLY=FEN_CO
SAMPLES=`awk -F '\t' '{if ($4 == "fen") {print $1}}' sample_metadata.tsv | uniq`
R1=`printf '%s\n' $SAMPLES | awk -v ORS="," '{print "POOLED_ILLUMINA/" $0 ".R1.fastq"}' | sed 's/,$/\n/'`
R2=`printf '%s\n' $SAMPLES | awk -v ORS="," '{print "POOLED_ILLUMINA/" $0 ".R2.fastq"}' | sed 's/,$/\n/'`
```

### Assemble reads with MEGAHIT

```bash
mkdir ASSEMBLIES

megahit -1 $R1 \
        -2 $R2 \
        --out-dir ASSEMBLIES/$ASSEMBLY \
        --min-contig-len 1000 \
        --k-min 57 \
        --k-max 157 \
        --k-step 10 \
        --memory 0.8 \
        --num-cpu-threads $NTHREADS
```

## Nanopore data

### Define assembly

```bash
# Sample m11216
ASSEMBLY=M11216_NANO
SAMPLE=m11216

# m12208
ASSEMBLY=M12208_NANO
SAMPLE=m12208
```

### Assemble reads with metaFlye

```bash
flye --nano-raw TRIMMED_NANOPORE/$SAMPLE.fastq \
     --out-dir ASSEMBLIES/$ASSEMBLY \
     --genome-size 5m \
     --threads $NTHREADS \
     --meta
```

### Polish contigs with pilon

```bash
# Map Illumina reads to the assembly with bowtie
bowtie2-build ASSEMBLIES/$ASSEMBLY/assembly.fasta \
              ASSEMBLIES/$ASSEMBLY/assembly

bowtie2 -1 POOLED_ILLUMINA/$SAMPLE.R1.fastq \
        -2 POOLED_ILLUMINA/$SAMPLE.R2.fastq \
        -S ASSEMBLIES/$ASSEMBLY/$SAMPLE.sam \
        -x ASSEMBLIES/$ASSEMBLY/assembly \
        --threads $NTHREADS \
        --no-unal

# Sort and index SAM files with SAMtools
samtools view -F 4 -bS ASSEMBLIES/$ASSEMBLY/$SAMPLE.sam | samtools sort > ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam
samtools index ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam

# Run pilon
java -Xmx128G -jar pilon-1.23.jar --genome ASSEMBLIES/$ASSEMBLY/assembly.fasta \
                                  --bam ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam \
                                  --outdir ASSEMBLIES/$ASSEMBLY \
                                  --output pilon \
                                  --threads $NTHREADS \
                                  --changes
```

## Check assemblies with metaQUAST

```bash
metaquast.py ASSEMBLIES/UPLAND_CO/final.contigs.fa ASSEMBLIES/FEN_CO/final.contigs.fa ASSEMBLIES/M11216_NANO/pilon.fasta ASSEMBLIES/M12208_NANO/pilon.fasta \
             -o ASSEMBLIES/METAQUAST \
             -t $NTHREADS \
             --fast
```

## Next step

Continue to [MAG binning](03-MAG-binning.md).
