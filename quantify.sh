#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=04:00:00
#SBATCH --mem 100000

# Load the software environment.
module purge
module load all
module load gencore/2
module load trinity/2.15.1
module load salmon/0.14.2
module load numpy/1.26.4
module load bowtie2/2.4.5
module load samtools/1.10
module load hisat2/2.2.1
module load rsem/1.3.3
module load r-base/4.2.3

# Estimate transcript abundances using RSEM/Bowtie2, the output will be in each folder, foreach sample.
align_and_estimate_abundance.pl \
--transcripts trinity_assembly.Trinity.fasta \
--seqType fq \
--samples_file cond.txt \
--est_method RSEM \
--thread_count 28 \
--gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
--prep_reference \
--aln_method bowtie2

# Create an expression matrix at the gene and isoform level from the RSEM estimations in the earlier step.
abundance_estimates_to_matrix.pl \
--est_method RSEM \
--gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--out_prefix RSEM.out \
fructose_rep1/RSEM.isoforms.results \
fructose_rep2/RSEM.isoforms.results \
glucose_rep1/RSEM.isoforms.results \
glucose_rep2/RSEM.isoforms.results \
pyruvate_rep1/RSEM.isoforms.results \
pyruvate_rep2/RSEM.isoforms.results

# Load the DESeq2 package (which also include an installation of EdgeR).
module load bioconductor-deseq2/1.32.0

# Run the Differential Gene Expression analysis with DESeq2 at the Gene level.
run_DE_analysis.pl \
--matrix RSEM.out.gene.counts.matrix \
--method DESeq2 \
 --samples_file dge.txt \
 --output deseq2_gene \
 --contrasts comp.txt

# Run the Differential Gene Expression analysis with DESeq2 at the Isoform level.
run_DE_analysis.pl \
--matrix RSEM.out.isoform.counts.matrix \
--method DESeq2 \
 --samples_file dge.txt \
 --output deseq2_isoform \
 --contrasts comp.txt

