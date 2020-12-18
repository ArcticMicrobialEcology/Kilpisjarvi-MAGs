# Metagenome assembling

## Illumina data

Based on the read-based analyses, samples appear too different for a good co-assembly.  
We will then do two co-assemblies, one with the fen samples only and the other with the remaining (upland) samples.  
This is done by assigning the variables $ASSEMBLY, $R1 and $R2.  

### Define assembly and create list of file names

```bash
# Assembly of upland samples
grep -v fen sample_metadata_illumina.txt | cut -f 1 > UPLAND_SAMPLES.txt

ASSEMBLY=UPLAND_CO

R1=`awk -v ORS="," '{print "TRIMMED_ILLUMINA/" $0 "R1.fastq"}' UPLAND_SAMPLES.txt | sed 's/,$/\n/'`
R2=`awk -v ORS="," '{print "TRIMMED_ILLUMINA/" $0 "R2.fastq"}' UPLAND_SAMPLES.txt | sed 's/,$/\n/'`

# Assembly of fen samples
grep fen sample_metadata_illumina.txt | cut -f 1 > FEN_SAMPLES.txt

ASSEMBLY=FEN_CO

R1=`awk -v ORS="," '{print "TRIMMED_ILLUMINA/" $0 ".R1.fastq"}' FEN_SAMPLES.txt | sed 's/,$/\n/'`
R2=`awk -v ORS="," '{print "TRIMMED_ILLUMINA/" $0 ".R2.fastq"}' FEN_SAMPLES.txt | sed 's/,$/\n/'`
```

### Assemble reads with MEGAHIT

```bash
mkdir ASSEMBLIES

megahit -1 `echo $R1` \
        -2 `echo $R2` \
        --out-dir ASSEMBLIES/$ASSEMBLY \
        --min-contig-len 1000 \
        --k-min 57 \
        --k-max 157 \
        --k-step 10 \
        --memory 0.8 \
        --num-cpu-threads $NTHREADS
```

## Nanopore data

Here we will do individual asssemblies for each sample.

### Define assembly

```bash
# Sample m11216
ASSEMBLY=M11216_NANO
SAMPLE=m11216

# m12208
ASSEMBLY=M12208_NANO
SAMPLE=m12208
```

### Assemble reads with Flye in metagenome mode

```bash
flye --nano-raw TRIMMED_NANOPORE/$SAMPLE.fastq \
     --out-dir ASSEMBLIES/$ASSEMBLY \
     --genome-size 5m \
     --threads $NTHREADS \
     --meta
```

### Polish contigs with pilon

```bash
# Map Illumina reads to assembly with bowtie
bowtie2-build ASSEMBLIES/$ASSEMBLY/assembly.fasta \
              ASSEMBLIES/$ASSEMBLY/assembly

bowtie2 -1 TRIMMED_ILLUMINA/$SAMPLE.R1.fastq \
        -2 TRIMMED_ILLUMINA/$SAMPLE.R2.fastq \
        -S ASSEMBLIES/$ASSEMBLY/$SAMPLE.sam \
        -x ASSEMBLIES/$ASSEMBLY/assembly \
        --threads $NTHREADS \
        --no-unal

# Sort and index SAM files
samtools view -F 4 -bS ASSEMBLIES/$ASSEMBLY/$SAMPLE.sam | samtools sort > ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam
samtools index ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam

# Run pilon
PILON_DIR=$HOME/pilon # Change here to the location of the pilon-1.23.jar file in your system

java -Xmx128G -jar $PILON_DIR/pilon-1.23.jar --genome ASSEMBLIES/$ASSEMBLY/assembly.fasta \
                                             --bam ASSEMBLIES/$ASSEMBLY/$SAMPLE.bam \
                                             --outdir ASSEMBLIES/$ASSEMBLY \
                                             --output pilon \
                                             --threads $NTHREADS \
                                             --changes
```

## Check assemblies with metaQUAST

```bash
metaquast.py ASSEMBLIES/UPLAND/final.contigs.fa ASSEMBLIES/FEN/final.contigs.fa ASSEMBLIES/m11216/pilon.fasta ASSEMBLIES/m12208/pilon.fasta \
             -o ASSEMBLIES/METAQUAST \
             -t $NTHREADS \
             --fast
```

## Next step

Continue to [MAG binning](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/master/04-MAG-binning.md).
