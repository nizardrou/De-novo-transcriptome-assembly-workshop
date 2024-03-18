#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=04:00:00
#SBATCH --mem 100000

# Load the software environment
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

