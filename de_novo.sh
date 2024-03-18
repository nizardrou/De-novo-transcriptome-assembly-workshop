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

# Perform de novo transcriptome assembly using Trinity, and supplying all the reads as input.
Trinity \
--seqType fq \
--samples_file cond.txt \
--CPU 28 \
--output trinity_assembly \
--max_memory 100G

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


# Load the Trinotate software environment
module purge
module load gencore
module load Miniconda3/4.7.10
source  activate /scratch/gencore/conda3/envs/Trinotate/
export PATH="/scratch/gencore/software/Trinotate-Trinotate-v4.0.2:$PATH"
export TRINOTATE_DATA_DIR=/scratch/gencore/trinotate_addons/databases/4.0.2

# Make a copy of the assembly FASTA to prevent accidental overwrite.
cp trinity_assembly.Trinity.fasta transcripts.fasta

# Find the longest ORFs in the assembly.
TransDecoder.LongOrfs -t transcripts.fasta

# Blast the longest ORFs vs the uniprot/swissprot DBs.
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /scratch/gencore/trinotate_addons/databases/4.0.2/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 28 > blastp.outfmt6

# Run HMMscan on the longest ORFs vs the PFAM domain DB.
hmmscan --cpu 28 --domtblout pfam.domtblout /scratch/gencore/trinotate_addons/databases/4.0.2/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep

# Incorporate the results of the BLAST and HMMscan result into the best ORF predictions.
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

# Initialize the Trinotate SQLite DB to include the Genes/Transcript ID mappings, the assembly FASTA and the homology supported longest ORFs.
Trinotate \
--db myTrinotate.sqlite \
--init \
--gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
--transcript_fasta transcripts.fasta \
--transdecoder_pep transcripts.fasta.transdecoder.pep

# Run the annotations using Trinotate and perform all of the analyses except SignalP6, this will also load the results into the Trinotate SQLite DB.
Trinotate --db myTrinotate.sqlite --CPU 28 \
--transcript_fasta transcripts.fasta \
--transdecoder_pep transcripts.fasta.transdecoder.pep \
--trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2/ \
--run "swissprot_blastp swissprot_blastx pfam tmhmmv2 EggnogMapper" \
--use_diamond 

# Run Signalp6 separately (to override the default eukaryotic mode) and select gram positive.
/scratch/gencore/trinotate_addons/signalp-4.1/signalp \
-f short \
-t gram+ \
-n signalp.out \
transcripts.fasta.transdecoder.pep

# Load the SignalP6 results into the Trinotate SQLite DB.
Trinotate --db myTrinotate.sqlite --LOAD_signalp signalp.out

# Unload the Trinotate environment and load an older version of Trinotate (since the report generation is not working in the latest version)
conda deactivate
module purge
module load gencore/1
module load gencore_rnaseq/1.0
module load gencore_annotation
module load gencore_dev
module load gencore_trinity/1.0

# Generate the Trinotate tab-delimited report (with an excel extension) from the Trinotate SQLite DB.
Trinotate myTrinotate.sqlite report > trinotate_annotation_report.xls

# Use the annotation attributes for the transcripts as ‘names’ for the transcripts in the Trinotate SQLite DB.
import_transcript_names.pl myTrinotate.sqlite trinotate_annotation_report.xls

# Load the latest Trinity software environment.
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

