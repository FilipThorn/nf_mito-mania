# nf_mito-mania
Nextflow pipeline that assembles mitochondria scaffolds using MITObim, checks the scaffold by mapping reads to it using bwa mem, and calling variants in freebayes to obtain a consensus sequence. 

## Workflow

1) Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04) \
   Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10) \
   Install [`mitobim`](https://github.com/chrishah/MITObim) \
   Install [`mira`](https://sourceforge.net/projects/mira-assembler/files/MIRA/stable/) 

2) Download git clone of this repository:
   ```bash
   git clone https://github.com/FilipThorn/nf_mito-mania
   ```
3) Edit nextflow.config file:
   ```bash
   mitobim = "/PATH/TO/MITObim.pl"                                    #path to MITObim script
   mitobimRef = "/PATH/TO/lycPyr_mtDNA.fa"                            #refernce for mitobim and mira
   ref_strain = "lycPyr_mtDNA"                                        #name of refernce for mitobim and mira
   mira =  "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/mira"     #path to mira
   mira_dir = "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/"      #path to mira dir
   ```
4) Run Mitomania workflow:
   ```bash
   nextflow run mitomania.nf --reads_MB /PATH/TO/Indivxxx_L001_U.fastq.gz --reads_PE '/PATH/TO/*_R{1,2}.fastq.gz' --reads_SE '/PATH/TO/*_U.fastq.gz' --outdir /PATH/TO/RESULTS
   ```
&nbsp;
&nbsp;
&nbsp;

*Currently requires both paired and unpaired reads to run as well as four seperate libraries per individual. Updates to generalize the script are ongoing.*  \
Input file structure
```bash
usr:~data/$ ls
Indivxxx_L001_U.fastq.gz Indivxxx_L001_R2.fastq.gz Indivxxx_L001_R1.fastq.gz 
Indivxxx_L002_U.fastq.gz Indivxxx_L002_R2.fastq.gz Indivxxx_L002_R1.fastq.gz
Indivxxx_L003_U.fastq.gz Indivxxx_L003_R2.fastq.gz Indivxxx_L003_R1.fastq.gz
Indivxxx_L004_U.fastq.gz Indivxxx_L004_R2.fastq.gz Indivxxx_L004_R1.fastq.gz 
```
&nbsp;
&nbsp;

## HPC enviroment
Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
 
