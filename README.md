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
TrinityStats.pl trinity_assembly.Trinity.fasta
```

This will produce,
```
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

We like to look at the number of genes and transcripts and the total transcriptome size in addition to the N50 stats. If we already have an expectation as to the total number of genes for our organism (from a closly related species perhaps), as well as the transcriptome size, this will give us an idea on how complete our transcriptome assembly is. At around 3K genes, we are in the right ballpark of what we expect.




