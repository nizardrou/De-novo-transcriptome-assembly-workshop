# De-novo-transcriptome-assembly-workshop
The repo is intended to accompany the NYUAD Core Bioinformatics hands-on workshop on de novo Transcriptome assembly and annotation.

Over the course of this 2-day workshop, participants will learn about:

- Assembling a de novo transcriptome using paired-end short reads (Illumina).
- Perform quality assessment on the resulting assembly.
- Annotating the assembly using ab initio gene finding as well as homology based searching.
- Perform RNAseq abundance estimation on a per sample basis using the assembly as a reference.

At the end of this 2-day workshop, it is expected that the participants should be able to replicate the steps involved using their own data as input. Although this workshop cannot possibly cover all the different approaches for assembling transcriptomes across a variety of species, it will provide a good reference to use for any transcriptome assembly and annotation analysis.

Although the analysis will focus on using the NYU Abu Dhabi High Performance Computing cluster (HPC), it can certainly be replicated on any other system provided that the neccessary software stack is preinstalled on that system.

We hope you find this workshop enjoyable and useful!



## Setting up the environment and copying the data
We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop, however, you can certainly run all of the analysis on any stand alone machine (server, personal laptop/Desktop etc.) provided that you have pre-installed the necessay software packages.

### Connecting to the HPC using a MAC/Linux machine and copying the data.
1. Open the "Terminal" app and type `ssh NetID@jubail.abudhabi.nyu.edu`. Enter your NYU password when prompted and hit enter.
2. Once logged in, navigate to your personal "SCRATCH" directory `cd $SCRATCH`.
3. Create a directory and change into it `mkdir de_novo_transcriptome && cd de_novo_transcriptome`.
4. Copy the data `cp -r /scratch/gencore/nd48/workshop_de_novo_trans/data .`.

### Connecting to the HPC using a Windows machine and copying the data.
1. Open the "Putty" app, and fill out the fields as follows **Host name**=jubail.abudhabi.nyu.edu, **Port**=22, and then click on "Open".
2. Enter your NetId, and your password when prompted.
3. Once logged in, navigate to your personal "SCRATCH" directory `cd $SCRATCH`.
4. Create a directory and change into it `mkdir de_novo_transcriptome && cd de_novo_transcriptome`.
5. Copy the data `cp -r /scratch/gencore/nd48/workshop_de_novo_trans/data .`.



## The software stack
As you might've imagined, this sort of analysis involves multiple steps, and multiple tools. The tools and their links are provided below in case you want to run this analysis on your own setup. We are using the HPC, so we don't have to install any of them since they are already there!

Just a quick note on installation. Whenever possible, we recommend using conda for installing and maintaining your software stack. Installing Bioinformatics software from source can be a painful experience sometimes, and conda takes care of most cases with relative ease.

- Trinity [https://github.com/trinityrnaseq/trinityrnaseq/wiki]
- Trinotate [https://github.com/Trinotate/Trinotate/wiki]
- Kallisto [https://github.com/pachterlab/kallisto]
- Augustus [https://bioinf.uni-greifswald.de/augustus/]
- Salmon [https://github.com/COMBINE-lab/salmon]
- rnaQuast [https://github.com/ablab/rnaquast]
- BUSCO [https://busco.ezlab.org/]
- SAMTools [https://www.htslib.org/]
- IGV [https://www.igv.org/]

## File Formats
We will be introducing a number of different file formats that are common in many types of omics analyses. Below is a list of these formats and links describing what they are.
It is assumed that you have some familiarity with these formats already.

- FASTA [https://en.wikipedia.org/wiki/FASTA_format]
- FASTQ [https://en.wikipedia.org/wiki/FASTQ_format]
- SAM [https://en.wikipedia.org/wiki/SAM_(file_format)] 
- BAM [https://en.wikipedia.org/wiki/Binary_Alignment_Map]

## The data
The datasets that we will be using for this workshop are publicly available sequencing data in the form of FASTQ sequencing files. The are publicly available to download from the Short Read Archive (SRA [https://www.ncbi.nlm.nih.gov/sra]) using the following accessions SRR28281136,SRR28281137,SRR28281138,SRR28281139,SRR28281140, and SRR28281141. There are 6 samples in total, organized as paired end sequences (so 12 files), and spanning 3 conditions with 2 replicates per condition. The organism is Staphylococcus aureus, and the 3 conditions represent bacterial growth under different carbon sources, Glucose, Fructose, and Pyruvate.


## Step 1: Assembling the reads using Trinity.
Almost all sequencing data analysis begin with quality checking and quality trimming of the data. This process ensures that we remove low quality sequenced nucleotides as well as other artefacts such as sequencing adapter contamination. And whilst this is also the first step in our analysis, we will not be addressing it here. However, if you would like to find out how that is performed, then please have a look at our Quality checking and Quality Trimming section of the Variant detection and Annotation Workshop here [https://github.com/nizardrou/Variant-Detection-and-Annotation-workshop/edit/main/README.md#step-1-data-quality-checking-and-quality-trimming].

We will be using Trinity for the purpose of assembling our transcriptome. Although there are other alternatives for transcriptome assembly, we have always produced good assemblies using this package. Another reason why we like Trinity is because it comes with a host of helper programs for post processing the assembly. It also performs well regardless of the organism size.
However, the golden rule here is that "no one size fits all"!
You should consider the fact that no tool is perfect and as such it is highly recommended for you to run multiple assemblies and assess the best performing tool.


When you copied the data to your $SCRATCH directory, included with the FASTQ files is the shell script "de_novo.sh". You will also find a FASTA file called "trinity_assembly.Trinity.fasta".
The "trinity_assembly.Trinity.fasta" is the pregenerated Trinity assembly. The reason that we have done that is due to time constraints. Assembly is a resources intensive (CPU and Memory) and time consuming process, and as such, it is not practicle for us to run the assembler and wait for 2+ hours for it to finish. That's not to say that we will not be showing you (below) how to run the assembly program.

The shell script that is provided "de_novo.sh" contains all the necessary analysis steps, which we will be addressing in this workshop. Let's start by examining it.

The first few lines are SLURM specific parameters. SLURM is a job scheduler and it will instruct the HPC where to submit our job. These instructions will be completely ignored if you are not submitting this script using "sbatch" so there is no need to delete them if you are running this on your own setup.

Since de novo assembly is a resource intensive task, we are asking for 100GB of memory and 28 CPUs. The larger the dataset, and the larger the transcriptome of an organism, the more resources you will need to allocate. It is a balancing act but it is always better to ask for more resources than what you need. Typically, you should estimate at least 1GB of memory for every 1 million reads, and in case you have too many reads/datasets, it is recommended that you run digital normalization on the reads before assembling [https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization].

### SLURM arguments
Here are the first few SLURM lines ,
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gencore
#SBATCH --cpus-per-task=28
#SBATCH --time=04:00:00
#SBATCH --mem 100000
```

### Loading the software stack
The next few lines load the neccessary software packages so that they will be available for us when we call them. These are specific to our setup at the NYUAD HPC, so you will need to change these according to your own environment. Also note that sometimes the parameters and arguments change from one software version to the next, so always check to see if the same parameters are applicable in case you are running a different version.

```
module purge
module load all
module load gencore/2
module load trinity/2.15.1
module load salmon/0.14.2
module load numpy/1.26.4
module load bowtie2/2.4.5
module load samtools/1.10
module load hisat2/2.2.1
```

### Running the Trinity assembly
The lines below will run the Trinity assembly. We are specifying that the input reads are in FASTQ format (they can also be FASTA), we are limiting the memory usage to 100GB and ask for 28 CPUs, that the output should be prefixed with "trinity_assembly", and that the samples their FASTQ files and their grouping is described in the tab-delimited file called "cond.txt".

```
Trinity \
--seqType fq \
--samples_file cond.txt \
--CPU 28 \
--output trinity_assembly \
--max_memory 100G
```

The "cond.txt" file is tab-delimited and should follow the order of,
```
CONDITION  CONDITION_REPLICATE_NUMBER  READ1_FASTQ_FILE_NAME  READ2_FASTQ_FILE_NAME
```

And that looks like this,
```
fructose    fructose_rep1   fru_1_1.fastq   fru_1_2.fastq
fructose    fructose_rep2   fru_2_1.fastq   fru_2_2.fastq
glucose glucose_rep1    glu_1_1.fastq   glu_1_2.fastq
glucose glucose_rep2    glu_2_1.fastq   glu_2_2.fastq
pyruvate    pyruvate_rep1   pyr_1_1.fastq   pyr_1_2.fastq
pyruvate    pyruvate_rep2   pyr_2_1.fastq   pyr_2_2.fastq

```

Although there are other ways to specify these inputs, we find this to be the best since it will allow us to reuse it later on when we quantify the expression levels for each sample after the assembly is complete.

### Submitting the assembly
All you have to do here is submit the script to SLURM using the command below,
```
sbatch de_novo.sh
```

And if you are running this on the command line rather than SLURM,
```
chmod 755 de_novo.sh
./de_novo.sh
```

Since this can take at least 2 hours, and our time is limited, we will be skipping this stage, hence why the commands for Trinity are commented out (lines starting with "#").


## Assessing the assembly Quality
After your assembly is complete, and before moving along with annotation and quantification, you want to make sure that your assembly performed well!
There are a few ways to assess assemblies, from stat based assessments, to the absence/presence of core genes, to alignment based metrics, and we will be covering these here.

### Assess using read mapping
In ideal situations, we want the reads that produced our assembly to map back to that assembly to the highest degree possible. A lower mapping percentages might indicate issues such as low sequencing yeild and not enough coverage, contamination, poor quality, or sequencing library quality issues. It might also indicate that our transcriptome is a lot more complex and that the assembler struggles to piece things together. It might also indicate that we have partially spliced genes so perhaps our extraction protocol needs to be investigated.

Mapping the reads back involves using a plice-aware read aligner such as Bowtie2, STAR, HISAT2 or similar. We will be using HISAT2 and SAMTools for this task.

The steps involved are below (included in the de_novo.sh script),

1. Index the transcriptome using HISAT2,
```
hisat2-build -p 28 trinity_assembly.Trinity.fasta trinity_assembly
```

2. Align the reads,
```
hisat2 \
-x trinity_assembly \
-p 28 \
-1 fru_1_1.fastq,fru_2_1.fastq,glu_1_1.fastq,glu_2_1.fastq,pyr_1_1.fastq,pyr_2_1.fastq \
-2 fru_1_2.fastq,fru_2_2.fastq,glu_1_2.fastq,glu_2_2.fastq,pyr_1_2.fastq,pyr_2_2.fastq \
-S aln.sam
```

3. Convert the SAM to a BAM, coordinate sort the BAM, index the sorted BAM, and gather the overall alignment stats,
```
samtools view -@ 28 -b -o aln.bam aln.sam
samtools sort -@ 28 -o aln.sorted.bam aln.bam
samtools index -@ 28 aln.sorted.bam
samtools flagstat aln.sorted.bam > stats.txt
```

Once these steps are complete, have a look at the "stats.txt" file and look for the overall alignment rates. Ideally we should see over 80% alignment rates, ours are ~93%  so we are showing excellent signs that most of our reads made it into the assembly!
```
61353459 reads; of these:
  61353459 (100.00%) were paired; of these:
    10998445 (17.93%) aligned concordantly 0 times
    36304307 (59.17%) aligned concordantly exactly 1 time
    14050707 (22.90%) aligned concordantly >1 times
    ----
    10998445 pairs aligned concordantly 0 times; of these:
      1097373 (9.98%) aligned discordantly 1 time
    ----
    9901072 pairs aligned 0 times concordantly or discordantly; of these:
      19802144 mates make up the pairs; of these:
        8727873 (44.08%) aligned 0 times
        7094790 (35.83%) aligned exactly 1 time
        3979481 (20.10%) aligned >1 times
92.89% overall alignment rate
```

There are dedicated alignment QC software that we can use in addition to what was described above, and get more detailed breakdown of the stats (such as Qualimap), but for this initial assessment, the approach above will suffice.

### Contig Nx and ExN50 Statistics
As we mentioned earlier, we like using Trinity because of the host of helper programs that come packaged with the software. These include a dedicated script to calculate Contig N50 values as well as the NxN50 that takes into account the read alignments and expression of the transcripts.
The N50 statistic is essentially the mean or median of lengths, and it is defined as the sequence length of the shortest contig at 50% of the total assembly length. Essentially, it tries to represent how contigeous our assembly is, with higher numbers" potentially" indicating a better assembly. We say potentially, because this statistic does nothing to express how well the assembly is put together, and so it should never be used as a sole metric, but rather an indicator.

To run the N50 stat,
```
TrinityStats.pl trinity_assembly.Trinity.fasta > assembly_stats.txt
```

We redirected the output of the command above ">" to a text file called "assembly_stats.txt". Let's examine it,
```
more assembly_stats.txt

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	3040
Total trinity transcripts:	3794
Percent GC: 34.52

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 14232
	Contig N20: 10849
	Contig N30: 8233
	Contig N40: 5784
	Contig N50: 4271

	Median contig length: 310
	Average contig: 1122.70
	Total assembled bases: 4259519


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 14232
	Contig N20: 10599
	Contig N30: 7935
	Contig N40: 5575
	Contig N50: 4175

	Median contig length: 302
	Average contig: 1097.71
	Total assembled bases: 3337029
```

We like to look at the number of genes and transcripts and the total transcriptome size in addition to the N50 stats. If we already have an expectation as to the total number of genes for our organism (from a closly related species perhaps), as well as the transcriptome size, this will give us an idea on how complete our transcriptome assembly is. At around 3K genes, we are in the right ballpark of what we expect. The total assembled bases is around 4.2Mb, which is more than what we would expect for this organism. Again, this might be because we have some background genomic (DNA) "contamination" in our sample, or because the sequencing has also captured plasmids.

Another more suitable statistic to calculate is the Ex90N50 and the Ex90 Gene Count. As explained in the Trinity pages, [https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#contig-ex90n50-statistic-and-ex90-gene-count],

_An alternative to the Contig Nx statistic that could be considered more appropriate for transcriptome assembly data is the ExN50 statistic. Here, the N50 statistic is computed as above but limited to the top most highly expressed genes that represent x% of the total normalized expression data. The gene expression is take as the sum of the transcript isoform expression and the gene length is computed as the expression-weighted mean of isoform lengths._

We will be generating this statistic later on as it requires us to align the samples back to our assembly and calculate expression values.


### Assessing the assembly using BUSCO

BUSCO stands for "Benchmarking Universal Single-Copy Orthologs", and it is used as a means to assess genome and transcriptome assemblies by looking for the presense of single-copy core genes within our assemblies. It can run on both eukaryotic and prokaryotic assemblies, and the higher the BUSCO scores, the better the assembly.
In contrast to what we have done previously, BUSCO is the first tool in our assessment list that looks at the content of our assembly as a means of assessing its quality.
It wraps around many packages in order to perform this task, and it can also run gene prediction software as part of the analysis. To know more about BUSCO follow the link at the "Sfotware Stack" section above.

To run BUSCO, we first need to load the appropriate package,
```
module purge
module load all
module load gencore/2
export PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/bin:$PATH"
export PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/bin:$PATH"
export AUGUSTUS_CONFIG_PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/config/"
module load busco/5.2.0
```

To run the BUSCO command,
```
busco \
-i trinity_assembly.Trinity.fasta \
--out BUSCO \
-m tran \
--auto-lineage-prok \
-c 28
```


In the command above the input is the Trinity assembly "-i trinity_assembly.Trinity.fasta", we specified the output folder to be "BUSCO" using the "--out BUSCO" flag, we also instructed BUSCO that we are running in Transcriptome mode and our organism is a prokaryote using the "-m tran" and "--auto-lineage-prok" flags respectively, and finally that we want BUSCO to ustilize 28 CPUs "-c 28".

The BUSCO output will be in the BUSCO folder.

### Assessing the assembly using QUAST
The final package that we will be using to assess the quality of our assembly is QUAST. QUAST stands for "Quality Assessment Tool for Genome Assemblies", and it is another common package that address the quality assessment of de novo genome and transcriptome assemblies.
In addition to statistics that we have generated earlier for N50, QUAST tries to capture potential misassemblies. Misassemblies are instances whereby contigs have been "stiched" together even though they shouldn't be as a result of the assembly program not being able to resolve them. These can occur around repetitive sequences and low complexity regions, or because we did not have sufficient coverage, read length etc. Misassemblies can be detected by mapping the reads back to the assembly and essentially looking for areas of 0 or very low coverage. If you think about it, the assembly should be a product of your reads, and if you have assembled a region that has no coverage by your reads, then it shouldn't exist!

QUAST can also wrap around ab initio gene finding and BUSCO, so it makes sense to execute QUAST at the same time.

We want to ask QUAST to map the reads back to the assembly, but we have 6 read1s and 6 read2s spanning all of our samples, and unlike HISAT2, we cannot supply them as a comma-separated list. The work around is to concatenate all of the read1s together and all of the read2s together.
```
cat *_1.fastq > read1.fastq
cat *_2.fastq > read2.fastq
```

Now, let's first load the neccessary packages,
```
module purge
module load gencore
module load Miniconda3/4.7.10
source activate /scratch/gencore/conda3/envs/quast_env
export PATH="/scratch/gencore/software/genemarks_t:$PATH"
```

The QUAST command should look something like this,
```
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
```

The command above can be explained as follows,
- -t 28 = number of CPUs to use.
- -f = enable gene finiding (gene predicition). This uses GeneMark for prokaryotes by default which is what we need.
- -b = run BUSCO.
- -o QUAST = The output folder name.
- --rna-finding = detect ribosomal RNA genes.
- --bam aln.sorted.bam = the BAM alignment file that we generated previously. This is used for structural variation detection and coverage histogram building.
- -1 read1.fastq -2 read2.fastq = the read1s that we concatenated earlier, which will be aligned to the assembly in order to detect misassemblies and coverage.
- trinity_assembly.Trinity.fasta = the assembly.

The output will be generated in a folder called "QUAST".

## Annotating the assembly using Trinotate
So now that we have checked the quality of our assembly and we are happy to proceed, the next logical step is to annotate the assembly. Broadly speaking, there are 2 ways that you can go about it.
1. Ab initio gene finding, which is what BUSCO and QUAST used in the previous steps to predict genes based on the sequence of our assembly (such as AUGUSTS, GeneMark, Glimmer etc).
2. Alignment and homology based searching, which is what Trinotate uses to annotate the genes.

It is important to note however that, Trinotate will take as input what you tell it are the genes! That is to say it will not try to predict where genes are, it will assume that you have done the work and you are providing it with a list of genes in FASTA format.

In our case, the assembly represents that gene list in FASTA format already, so we can go ahead and provide it to Trinotate.

So what does Trinotate do?

At its core, Trinotate provides the means for functional annotation of transcriptomes, particularly de novo transcriptomes. It does so by comparing against known homology databases such as PFAM and UniProt, GeneOntolgy, and EggNogg. Similar to Trinity, it also provides many helper scripts that allow for easy integration with Trinity and downstream expression quantification. 

Please note that as of March 2024, Trinotate is no longer under active development.



Let's start by loading the packages,
```
module purge
module load gencore
module load Miniconda3/4.7.10
source  activate /scratch/gencore/conda3/envs/Trinotate/
export PATH="/scratch/gencore/software/Trinotate-Trinotate-v4.0.2:$PATH"
export TRINOTATE_DATA_DIR=/scratch/gencore/trinotate_addons/databases/4.0.2
```

Let's make a copy of our FASTA assembly just to be on the safe side (in case we overwrite something, trust we've done it before!),
```
cp trinity_assembly.Trinity.fasta transcripts.fasta
```

Next, we need to run Transdecoder and identify the longest ORFs for our candidate genes (from the FASTA assembly),
```
TransDecoder.LongOrfs -t transcripts.fasta
```

Transdecoder will identify the longest Open Reading Frame from our transcripts, but this is inferred based on the sequence content. It will be better if we can add some supporting evidence for these ORFs. We can do this through Transdecoder but we must first carry out some homology searches on our longest ORFs.

To do this, we will search the uniprot and swissprot databases using BlastP, as well as against the PFAM domain database using HMMScan.
```
blastp \
-query transcripts.fasta.transdecoder_dir/longest_orfs.pep \
-db /scratch/gencore/trinotate_addons/databases/4.0.2/uniprot_sprot.pep \
-max_target_seqs 1 \
-outfmt 6 \
-evalue 1e-5 \
-num_threads 28 > blastp.outfmt6

hmmscan \
--cpu 28 \
--domtblout pfam.domtblout \
/scratch/gencore/trinotate_addons/databases/4.0.2/Pfam-A.hmm \
transcripts.fasta.transdecoder_dir/longest_orfs.pep
```

And then we will need to predict the ORFs again but only retaining the ones with hits,
```
TransDecoder.Predict \
-t transcripts.fasta \
--retain_pfam_hits pfam.domtblout \
--retain_blastp_hits blastp.outfmt6
```
The output from the steps above will be the following,
```
transcripts.fasta.transdecoder.bed
transcripts.fasta.transdecoder.cds
transcripts.fasta.transdecoder.gff3
transcripts.fasta.transdecoder.pep
```

Now that we are done with refining our ORF prediction, and we utilized an evidence based approach, we can carry on with the rest of the Trinotate steps.

In your data folder, you will find an SQLite database (myTrinotate.sqlite) that was created using Trinotate. This was initialized during the setup of Trinotate and it only needs to happen once, so you don't have to do this unless you are setting up Trinotate on your own setup.

The command below is what created this Database,
```
Trinotate --create \
--db myTrinotate.sqlite \
--trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2 \
--use_diamond
```

When the command above is initiated, it will download all of the databases that Trinotate relies on (such as PFAM, uniprot/swissprot, eggnogg etc.), it will build the databases, and initialize the database. Once this "myTrinotate.sqlite" database is created, you can copy it somewhere and reuse for future projects anytime you like. If at any point you need to update the databases, then all you have to do is run the above command and start using the newly create database.

Import the transcripts and their homolgy predicted protein translations into the database,
```
Trinotate \
--db myTrinotate.sqlite \
--init \
--gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
--transcript_fasta transcripts.fasta \
--transdecoder_pep transcripts.fasta.transdecoder.pep
```

Run the Trinotate command (more information below),
```
Trinotate --db myTrinotate.sqlite --CPU 28 \
--transcript_fasta transcripts.fasta \
--transdecoder_pep transcripts.fasta.transdecoder.pep \
--trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2/ \
--run "swissprot_blastp swissprot_blastx pfam tmhmmv2 EggnogMapper" \
--use_diamond
```

The command above will do the following,
- --CPU 28 = Use 28 CPUs.
- --transcript_fasta transcripts.fasta = Provide the transcripts (assembly) as input.
- --transdecoder_pep transcripts.fasta.transdecoder.pep = Provide the Transdecoder ORFs as input as well.
- --trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2/ = The path to the Trinotate database directory (the same one used when we initialized it above).
- --run "swissprot_blastp swissprot_blastx pfam tmhmmv2 EggnogMapper" = Instruct Trinotate, which analysis steps we want to perform.
- --use_diamond = Use "diamond blast" instead of "NCBI blast", which significantly speeds up the searches.

The command above will automatically load the results of the analysis into the SQLite database, so we don't have to do that ourselves (previous versions of Trinotate required that you do that after independently running the commands).

You might notice that we are missing "SignalP6" from our list above. SignalP6 predictis signal peptides and the location of their cleavage sites in proteins, and we can certainly ask Trinotate to run it at the same time as the other tools. However, by default, Trinotate will run SignalP6 in "eukaryotic" mode, which is not what we want, since we are dealing with a Gram positive bacteria.

But we can run SignalP6 independently ourselves and adjust the parameters to suite our needs like so,
```
/scratch/gencore/trinotate_addons/signalp-4.1/signalp \
-f short \
-t gram+ \
-n signalp.out \
transcripts.fasta.transdecoder.pep
```

Then all we have to do is load the results into the "myTrinotate.sqlite" database,
```
Trinotate --db myTrinotate.sqlite --LOAD_signalp signalp.out
```

At the end of this annotation process, it is useful to produce a tab-delimited file that includes all of the results of our analysis. This is achieved through the "Trinotate report" command.

**HOWEVER**, we have found out (and this might only be relevant to our own setup and installation), that the report generation in the latest version of Trinotate (the one that we have been using so far) does not work, and it just produces and empty report!

The work-around that we have found, is to generate the report using a previous version of Trinotate, but to do that, we have to unload our current environment and load the previous version.

```
conda deactivate
module purge
module load gencore/1
module load gencore_rnaseq/1.0
module load gencore_annotation
module load gencore_dev
module load gencore_trinity/1.0

Trinotate myTrinotate.sqlite report > trinotate_annotation_report.xls
import_transcript_names.pl myTrinotate.sqlite trinotate_annotation_report.xls
```

The work-around above is a perfect example of why setting and maintaining your software stack through conda is so useful! If you just had a source installation of Trinotate, it meant that you would now have to uninstall the current version and install the previous one.

And with that we are done with Trinotate, so let's take a few moments and examine the contents of the report (you can open it in EXCEL).

## Expression Quantification
So far, we have assembled our transcriptome, assessed the quality of our assembly, and annotated our assembly. The last step of this analysis involves us quantifying the expression values and comparing our conditions and replicates. This is similar to what we do with a "traditional" RNAseq analysis, and from this point on, you can certainly move ahead with your own prefered RNAseq methodolgy, however, since Trinity/Trinotate provide helper scripts to automate this and work intuitively with the outputs that we have generated so far, we will stick with that.

Trinity provides support for several differential expression analysis tools, currently including the following R packages:

- edgeR : [http://bioconductor.org/packages/release/bioc/html/edgeR.html]
- DESeq2: [http://bioconductor.org/packages/release/bioc/html/DESeq2.html]
- limma/voom: [http://bioconductor.org/packages/release/bioc/html/limma.html]
- ROTS: [http://www.btk.fi/research/research-groups/elo/software/rots/]

**PLEASE NOTE** that in our example data, we only have duplicates, not triplicates for each condition, and so we will not be able to perform statistical analysis on differentially expressed genes.

First, let's go ahead and switch back to loading Trinity,
```
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

```

Before proceeding with DE (Differential expression) analysis, we first need to compute a counts matrix. This will tell us how much of each assembled transcript, each one of our samples is "expressing". There are 2 main approaches to this, either alignment based (such as RSEM, HISAT2, BOWTIE2 etc.), or alignment-free methods that rely on k-mer estimations (such as Kallisto and Salmon). The pros/cons of each method are beyond the scope of this tutorial but here is a link to a paper that tries to benchmark them against each other [https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4869-5#:~:text=These%20alignment%2Dfree%20pipelines%20are,indexed%20transcript%20databases%20%5B4%5D.].

Trinity is packaged with a script "align_and_estimate_abundance.pl " that allows us to automoate this process.
```
align_and_estimate_abundance.pl \
--transcripts trinity_assembly.Trinity.fasta \
--seqType fq \
--samples_file cond.txt \
--est_method RSEM \
--thread_count 28 \
--gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
--prep_reference \
--aln_method bowtie2
```

In the command above, we are instructing the script to estimate transcript abundance using RSEM [https://github.com/deweylab/RSEM], and utilize BOWTIE2 for this purpose. We prefer running the alignment before-hand and providing the BAM file(s) to RSEM, but we wanted to show you the process within Trinity.
Adding the "--prep_reference" will also index the assembly using BOWTIE2 as a first step, in case you already indexed the assembly, you can skip this flag.
The output will be generated on a per-sample basis, that is, a folder for each replicate will be created with RSEM outputs in each folder.

Let's look at the output of one of these samples for example,
```
ls -1 fructose_rep1/
bowtie2.bam
bowtie2.bam.for_rsem.bam
bowtie2.bam.ok
RSEM.genes.results
RSEM.isoforms.results
RSEM.isoforms.results.ok
RSEM.stat
```

We have the BAM alignment files, but in addition we have the 2 files that we are most interested in, the "RSEM.genes.results", and the "RSEM.isoforms.results", representing transcript quantification at the gene level and at the transcript level respectively.

Now that we have the individual results for each replicate, we will need to create an expression matrix for the whole experiment (all replicates). We can do that using the "abundance_estimates_to_matrix.pl " script from Trinity.

```
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
```

This step will compute the expression matrix for both genes and isoforms, so you don't need to run separately for RSEM.isoforms.results and RSEM.genes.results files (in fact if you supply the genes then it won't work and give you an error).

The files that are created are,
```
ls -1 RSEM.out.*
RSEM.out.gene.counts.matrix
RSEM.out.gene.TMM.EXPR.matrix
RSEM.out.gene.TPM.not_cross_norm
RSEM.out.gene.TPM.not_cross_norm.runTMM.R
RSEM.out.gene.TPM.not_cross_norm.TMM_info.txt
RSEM.out.isoform.counts.matrix
RSEM.out.isoform.TMM.EXPR.matrix
RSEM.out.isoform.TPM.not_cross_norm
RSEM.out.isoform.TPM.not_cross_norm.runTMM.R
RSEM.out.isoform.TPM.not_cross_norm.TMM_info.txt
```

Now that we have our expression matrices generated, the next logical step is to run DGE (Differential gene expression analysis). At this point, you can branch off and do your own preferred method, by for the purpose of this workshop, we are going to stay within the Trinity package.

Trinity provides different options for running DGE as highlighted earlier. We are going to proceed with DESeq2.

To do this, first we need to run the script "run_DE_analysis.pl",
```
module load bioconductor-deseq2/1.32.0

run_DE_analysis.pl \
--matrix RSEM.out.gene.counts.matrix \
--method DESeq2 \
 --samples_file dge.txt \
 --output deseq2 \
 --contrasts comp.txt
```

An alternative would be to use our own instance of NASQAR2 [http://nasqar2.abudhabi.nyu.edu/], and select the DESeq2 app [http://nasqar2.abudhabi.nyu.edu/deseq2shiny/].

The output folder is called "deseq2", so let's examine the content!

