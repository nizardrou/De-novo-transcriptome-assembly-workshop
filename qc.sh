#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=04:00:00
#SBATCH --mem 100000

# Load modules for Trinity.
module purge
module load all
module load gencore/2
module load trinity/2.15.1
module load salmon/0.14.2
module load numpy/1.26.4
module load bowtie2/2.4.5
module load samtools/1.10
module load hisat2/2.2.1


## Perform de novo transcriptome assembly using Trinity, and supplying all the reads as input.
#Trinity \
#--seqType fq \
#--samples_file cond.txt \
#--CPU 28 \
#--output trinity_assembly \
#--max_memory 100G

# Index the assembly using hisat2.
hisat2-build -p 28 trinity_assembly.Trinity.fasta trinity_assembly

# Align the reads and capture the alignment stats
hisat2 \
    -x trinity_assembly \
    -p 28 \
    -1 fru_1_1.fastq,fru_2_1.fastq,glu_1_1.fastq,glu_2_1.fastq,pyr_1_1.fastq,pyr_2_1.fastq \
    -2 fru_1_2.fastq,fru_2_2.fastq,glu_1_2.fastq,glu_2_2.fastq,pyr_1_2.fastq,pyr_2_2.fastq \
    -S aln.sam

# Convert the SAM to a BAM, coordinate sort the BAM, index the sorted BAM, and gather the overall alignment stats.
samtools view -@ 28 -b -o aln.bam aln.sam
samtools sort -@ 28 -o aln.sorted.bam aln.bam
samtools index -@ 28 aln.sorted.bam
samtools flagstat aln.sorted.bam > stats.txt

# Calculate N50 stats
TrinityStats.pl trinity_assembly.Trinity.fasta > stats.txt

# Run BUSCO
module purge
module load all
module load gencore/2
export PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/bin:$PATH"
export PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/bin:$PATH"
export AUGUSTUS_CONFIG_PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/config/"
module load busco/5.2.0

busco \
    -i trinity_assembly.Trinity.fasta \
    --out busco \
    -m tran \
    --auto-lineage-prok \
    -c 28

# Concatenate all read1s together and all read2s together
cat *_1.fastq > read1.fastq
cat *_2.fastq > read2.fastq

# Run QUAST
module purge
module load gencore
module load Miniconda3/4.7.10
source activate /scratch/gencore/conda3/envs/quast_env
export PATH="/scratch/gencore/software/genemarks_t:$PATH"

quast.py \
    -t 28 \
    -f \
    -b \
    -o QUAST \
    --rna-finding \
    --bam aln.sorted.bam \
    -1 read1.fastq \
    -2 read2.fastq \
    trinity_assembly.Trinity.fasta


