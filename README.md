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
