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

