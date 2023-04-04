# nf_mito-mania
Nextflow pipeline that assembles mitochondria scaffolds using MITObim, checks the scaffold by mapping reads to it using bwa mem, and calling variants in freebayes to obtain a consensus sequence. Workflow is based on a [`Mozes Blom`](https://github.com/MozesBlom/mitogenome) pipeline.

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
   mitobim = "/PATH/TO/MITObim.pl"                                     #path to MITObim script
   mitobimRef = "/PATH/TO/MTDNASEED_ref.fa"                            #refernce for mitobim and mira
   ref_strain = "MTDNASEED_ref"                                        #name of refernce for mitobim and mira
   mira =  "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/mira"      #path to mira
   mira_dir = "/PATH/TO/mira_4.0.2_linux-gnu_x86_64_static/bin/"       #path to mira dir
   ```
4) Input tab separated file:
  ```bash 
  Individual   Unpaired_fn R1_fn R2_fn n<br>
ind1 /Absolute/PATH/ind1_unpaired_lib1_reads.fa.gz /Absolute/PATH/ind1_lib1_R1.fa.gz /Absolute/PATH/ind1_lib1_R2.fa.gz  n<br>
ind1 /Absolute/PATH/ind1_unpaired_lib2_reads.fa.gz /Absolute/PATH/ind1_lib2_R1.fa.gz /Absolute/PATH/ind1_lib2_R2.fa.gz  n<br>
ind2 /Absolute/PATH/ind2_unpaired_lib1_reads.fa.gz /Absolute/PATH/ind2_lib1_R1.fa.gz /Absolute/PATH/ind2_lib1_R2.fa.gz  n<br>
ind2 /Absolute/PATH/ind2_unpaired_lib2_reads.fa.gz /Absolute/PATH/ind2_lib2_R1.fa.gz /Absolute/PATH/ind2_lib2_R2.fa.gz  n<br>


5) Run Mitomania workflow:
   ```bash
   nextflow run mitomania.nf --input_tsv_fn input_tsv --outdir /PATH/TO/RESULTS
   ```
&nbsp;
&nbsp;

## HPC enviroment
Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
 
