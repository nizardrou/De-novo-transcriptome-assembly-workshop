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
Trinity \
--seqType fq \
--samples_file cond.txt \
--CPU 28 \
--output trinity_assembly \
--max_memory 100G


