#!/usr/bin/env nextflow

if (params.help) {
    log.info """\
            __________________
            |                |
            | |```````| |`````
            | |____   | |
            |     |   | |
            | |````   | |
            | |ilip   | |hörn     
            –––––––––––––––––––––––––––––––––––––––
            Mitomania
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            Modified 17th of March 2022

            'USAGE'
            nextflow run mitomania.nf --reads_MB /PATH/TO/Indivxxx_L001_U.fastq.gz --reads_PE '/PATH/TO/*_R{1,2}.fastq.gz' --reads_SE '/PATH/TO/*_U.fastq.gz' --outdir /PATH/TO/RESULTS         
            
            'Mandatory arguments:'
            --reads_MB               Path to unpaired reads for MITObim backbone. i.e Indivxxx_L001_U.fastq.gz
            --reads_PE               Path to paired reads 
            --reads_SE               Path to unpaired reads
            --outdir                 Path to output directory
            
            'OPTIONS'
            --help                   Outputs this help log      
            -resume                  Nextflow cmd to resume modified workflow
            --subset_size   INTEGER


            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 


            For mitobim see https://github.com/chrishah/MITObim
            For bwa see http://bio-bwa.sourceforge.net
            For freebayes see https://github.com/freebayes/freebayes
            For samtools see http://www.htslib.org
            For bcftools see https://samtools.github.io/bcftools/bcftools.html
            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         Mitomania 
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         reads_MB     : ${params.reads_MB}
         reads_PE     : ${params.reads_PE}
         reads_SE     : ${params.reads_SE}
         outdir       : ${params.outdir}
         subset_size  : ${params.subset_size}
         """
         .stripIndent()


channel
        .fromPath(params.input_tsv_fn)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.Individual, file(row.Unpaired_fn)) }
        .groupTuple(by: 0)
        .into { mitobim_ch; SE_ch }

channel
        .fromPath(params.input_tsv_fn)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.Individual, file(row.R1_fn), file(row.R2_fn)) }
        .groupTuple(by: 0)
        .set { PE_ch }

process Subsample {

    tag "$id"

    publishDir "${params.outdir}/1.SubsetReads/$id", mode:'copy'

    input:
    tuple val(id), file(bb) from mitobim_ch

    output:
    tuple val(id), file("${id}_subset.fastq" ) into subset_ch.into{ subset_ch; subset_comb; subset_comb_mira}

    script:
    """
    seqtk sample -s1234 $read $params.subset_size > ${id}_subset.fastq
    """
}

process manifest {
    
    tag "$id"
    
    label 'FAST'

    publishDir "${params.outdir}/2.Mitobim/$id", mode:'copy'
    
    input:
    tuple val(id), file(sub) from subset_ch

    output:
    tuple val(id), file("manifest.conf") into mani_ch.combine(subset_comb, by:0)

    script:
    """
    manifest.sh $id $params.mitobimRef $params.ref_strain "${params.outdir}/1.SubsetReads/$id/${sub}" manifest.conf
    """
}


process mira {
    
    tag "$id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/2.Mitobim/$id", mode:'copy'
    
    input:
    tuple val(id), file(mani) , file(sub) from mani_ch

    output:
    file ("${id}_assembly/${id}_d_chkpt/*")
    file ("${id}_assembly/${id}_d_info/*")
    file ("${id}_assembly/${id}_d_results/*")
    file ("${id}_assembly/${id}_d_tmp/*")
    tuple val(id), file("${id}_assembly/${id}_d_results/${id}_out.maf") into mira_ch.combine(subset_comb_mira, by:0)

    script:
    """
    $params.mira $mani
    """
}


process mitobim {
    
    tag "$id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/3.Index/$id", mode:'copy'
    publishDir "${params.outdir}/2.Mitobim/$id", pattern: '*.log', mode:'copy'

    input:
    tuple val(id), file(maf), file(sub) from reads_comb_mira
    
    output:
    file("${id}.log")
    tuple val(id), file("${id}_noIUPAC_bb.fa") into bb_ch.into{ bb_ch; bb2_ch; bb3_ch; bb4_ch }

    script:
    """
    $params.mitobim -start 1 -end 100 -sample $id -ref $params.ref_strain -readpool "${id}_subset.fastq" -maf "${id}.maf" --clean --mirapath $params.mira_dir &> ${id}.log

    cp ./iteration*/*_noIUPAC.fasta ./${id}_noIUPAC_bb.fa
    """
}

                 

process BWA_Index {

    tag "$id"

    publishDir "${params.outdir}/3.Index/$id", mode:'copy'

    input:
    tuple val(id), file(ref) from bb_ch

    output:
    tuple val("${id}"), file("${id}.fa"), path("*.fa.{amb,ann,bwt,pac,sa}") into index.into{index_SE; index_PE}

    script:
    """
    bwa index $ref -p ${id}.fa
    """
}


process BWA_mem_SE {    

    tag "$id"    
    input:
    tuple val(id), file(ref), file(index), file(reads) from index_SE.combined(SE_ch, by:0)

    output:
    tuple val(id), file("${id}_mito_sorted_SE.bam") into aln_SE

    script:
    """
    bwa mem $ref $reads -t ${task.cpus} -R "@RG\\tID:${id}\\tSM:${sample}\\tPL:ILLUMINA" > ${id}_mito_SE.sam 
    samtools view -bS ${id}_mito_SE.sam | samtools sort -@ ${task.cpus} -O bam -o ${id}_mito_sorted_SE.bam
    """
}


process BWA_mem_PE {

    tag "$id"
    input:
    tuple val(id), val(sample), file(ref), file(index), file(R1), file(R2) from index_PE.combined(PE_ch, by:0)

    output:
    tuple val(id), file("${id}_mito_sorted_PE.bam") into aln_PE

    script:
    """
    bwa mem $ref $R1 $R2 -t ${task.cpus} -R "@RG\\tID:${id}\\tSM:${sample}\\tPL:ILLUMINA" > ${id}_mito_PE.sam
    samtools view -bS ${id}_mito_PE.sam | samtools sort -@ ${task.cpus} -O bam -o ${id}_mito_sorted_PE.bam
    """
}


process SamtoolsSortMerge{
        
    tag "$id"
    input:
    tuple val(id), file(aln_SE), file(aln_PE) from aln_SE.combine(aln_PE, by:0)
    
    output:
    tuple val(id), file("${id}_mito_sorted.bam") into sam_ch.into{ free_ch; mask_ch}
    file("${id}_mito_sorted.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${id}_mito_sorted.bam $aln_PE $aln_SE 
    samtools index -b ${id}_mito_sorted.bam
    """
}


process VariantCall{

    publishDir "${params.outdir}/5.Variantcall/$id", mode:'copy'

    tag "$id"

    input:
    tuple val(id), file(bam), file(ref) from free_ch.combine(bb2_ch, by:0)

    output:
    file("${id}_ploidy2.vcf.gz")
    tuple val(id), file("${id}.vcf.gz") into vcf_ch

    script:
    """   
    freebayes -f $ref $bam > ${id}_ploidy2.vcf
    vcffilter -f "TYPE = snp & ( AB = 0 | AB < 0.1 )" ${id}_ploidy2.vcf > ${id}.vcf
    bgzip -c ${id}_ploidy2.vcf > ${id}_ploidy2.vcf.gz
    bgzip -c  ${id}.vcf > ${id}.vcf.gz
    """
}


process mask {

    publishDir "${params.outdir}/6.masks/$id", mode:'copy'
    
    tag "$id"

    input:
    tuple val(id), file(bam), file(ref) from mask_ch.combine(bb3_ch, by:0)
    
    output:
    tuple val(id), file("${id}_mask_maxDP_*.txt") into bed_ch
    file("*")
    
    script:
    """
    samtools depth -aa $bam -o ${id}_depth.txt

    mask.sh ${id}_depth.txt $id   
    """
}


process CallConsensus {

    publishDir "${params.outdir}/7.consensus/$id", mode:'move'
    
    tag "$id"

    input:
    tuple val(id), file(MB), file(vcf), file(mask) from bb3_ch.combine( vcf_ch, by: 0).combine( bed_ch, by: 0)

    output:
    file("*")

    script:
    """
    bcftools index $vcf
    bcftools consensus -f $MB -m $mask -o ${id}_consensus_mito.fa -s ${id} $vcf
    sed -i "s/>.*/>${id}/" ${id}_consensus_mito.fa
    """
}
