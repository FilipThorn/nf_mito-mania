# nf_mito-mania
Nextflow pipeline that assembles mitochondria scaffolds using mitobim and checked scaffold by mapping reads to it using bra mem and calls variants in freebayes for a consensus sequence. 

## Workflow

1) Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04)\n
   Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10)\n 
   Install [`mitobim`](https://github.com/chrishah/MITObim)\n
   Install [`mira`](https://sourceforge.net/projects/mira-assembler/files/MIRA/stable/)\n

2) Download git clone of this repository:
   ```bash
   git clone https://github.com/FilipThorn/nf_mito-mania
   ```
3) Edit nextflow.config file:
   ```bash
   mitobim = "/PATH/TO/MITObim.pl"
   mitobimRef = "/PATH/TO/lycPyr_mtDNA.fa"
   ref_strain = "lycPyr_mtDNA"
   mira =  "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/mira"
   mira_dir = "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/"
   ```
4) Make mitobim scaffold:
   ```bash
   nextflow run mitobim_reference.nf --reads_MB /PATH/TO/Indivxxx_L001_U.fastq.gz --outdir /PATH/TO/RESULTS
   ```
   *Check mitobim scaffolds* 
   
5) Map reads to mitobim scaffold, call varient sites and build consensus:
   ```bash
   nextflow run map_reads.nf --ref '/RESULTS/3.IndexRefs/Indivxxx_L001_U/*.fasta' --outdir /RESULTS --reads_PE '/*_R{1,2}.fastq.gz' --reads_SE '/*_U.fastq.gz'
   ```
 
 ## HPC enviroment
Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
 
