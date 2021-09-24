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


// Channel for MITObim reads
   Channel.fromPath( params.reads_MB)
         .map { it -> [it.name -  ~/\.fastq\.gz/, it] }
         .set { reads_MB_ch }

// Channels for bwa reads
reads_PE = Channel.fromPath(params.reads_PE)
                   .map {it -> [it.name.tokenize("_")[0,1].join('_'), it] }
                   .groupTuple(sort: true)



reads_SE = Channel.fromPath( params.reads_SE)
                .map {it -> [it.name.tokenize("_")[0,1].join('_'), it]}
                .groupTuple()   

process Subsample {

    tag "$sample_id"

    publishDir "${params.outdir}/1.SubsetReads/$sample_id", mode:'copy'

    input:
    tuple val(sample_id), file(read) from reads_MB_ch

    output:
    tuple val(sample_id), file( "${sample_id}_subset.fastq" ) into subset_ch

    script:
    """
    seqtk sample -s1234 $read $params.subset_size > ${sample_id}_subset.fastq
    """
}

subset_ch.into{ subset_ch; subset_comb; subset_comb_mira}

process manifest {
    
    tag "$sample_id"
    
    label 'FAST'

    publishDir "${params.outdir}/2.Mitobim/$sample_id", mode:'copy'
    
    input:
    tuple val(sample_id), file(sub) from subset_ch

    output:
    tuple val(sample_id), file("manifest.conf") into mani_ch

    script:
    """
    manifest.sh $sample_id $params.mitobimRef $params.ref_strain "${params.outdir}/1.SubsetReads/$sample_id/${sub}" manifest.conf
    """
}

subset_comb.combine( mani_ch, by: 0 )
    .set { reads_comb } 

process mira {
    
    tag "$sample_id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/2.Mitobim/$sample_id", mode:'copy'
    
    input:
    tuple val(sample_id), file("${sample_id}_subset.fastq") , file("manifest.conf") from reads_comb 

    output:
    file ("${sample_id}_assembly/${sample_id}_d_chkpt/*")
    file ("${sample_id}_assembly/${sample_id}_d_info/*")
    file ("${sample_id}_assembly/${sample_id}_d_results/*")
    file ("${sample_id}_assembly/${sample_id}_d_tmp/*")
    tuple val(sample_id), file("${sample_id}_assembly/${sample_id}_d_results/${sample_id}_out.maf") into mira_ch

    script:
    """
    $params.mira "manifest.conf"
    """
}

subset_comb_mira.combine( mira_ch, by: 0 )
    .set { reads_comb_mira } 

// define the regexp pattern
ID_post = ~/_L00\d_U/

process mitobim {
    
    tag "$sample_id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/3.Index/$Id", mode:'copy'
    publishDir "${params.outdir}/2.Mitobim/$sample_id", pattern: '*.log', mode:'copy'
    input:
    tuple val(sample_id), file("${sample_id}_subset.fastq") , file("${sample_id}.maf") from reads_comb_mira
    
    output:
    file("${sample_id}.log")
    tuple val(Id), file("${sample_id}-lycPyr_mtDNA-it*_noIUPAC.fasta") into mitobim_ch

    script:
    Id = ("${sample_id}" - ID_post)
    """
    $params.mitobim -start 1 -end 100 -sample $sample_id -ref $params.ref_strain -readpool "${sample_id}_subset.fastq" -maf "${sample_id}.maf" --clean --mirapath $params.mira_dir &> ${sample_id}.log

    cp ./iteration*/*_noIUPAC.fasta ./
    """
}

// pipe ref into seperate channels 
mitobim_ch.into {ref; ref2; ref3}
                 

// define the regexp pattern
lib = ~/_L00\d/

process BWA_Index {

    tag "$sample_id"

    publishDir "${params.outdir}/3.Index/$sample_id", mode:'copy'

    input:
    tuple val(sample_id), file(ref) from ref

    output:
    tuple val("${sample_id}_L001"), file("${sample_id}"), path("*.{amb,ann,bwt,pac,sa}") into index1
    tuple val("${sample_id}_L002"), file("${sample_id}"), path("*.{amb,ann,bwt,pac,sa}") into index2
    tuple val("${sample_id}_L003"), file("${sample_id}"), path("*.{amb,ann,bwt,pac,sa}") into index3
    tuple val("${sample_id}_L004"), file("${sample_id}"), path("*.{amb,ann,bwt,pac,sa}") into index4


    script:
    """
    bwa index $ref -p ${sample_id}
    cp $ref $sample_id
    """
}

//mix
index1.mix(index2,index3,index4).into{index_SE; index_PE}

//combine

bwa_SE = reads_SE.combine(index_SE, by:0)

bwa_PE = reads_PE.combine(index_PE, by:0)


process BWA_mem_SE {    

    tag "$sample_id"    
    
    input:
    tuple val(sample_id), file(sample), file(ref), file(index) from bwa_SE

    output:
    tuple val(sample_id), file("*") into aln_SE

    script:
    RG_id = ("${sample_id}" - lib)
    """
    bwa mem $ref $sample -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${RG_id}\\tPL:ILLUMINA" > ${sample_id}_mito_SE.sam 
    """
}


process BWA_mem_PE {

    tag "$sample_id"

    input:
    tuple val(sample_id), file(pairs) ,file(ref), file(index) from bwa_PE

    output:
    tuple val(sample_id), file("*") into aln_PE

    script:
    RG_id = ("${sample_id}" - lib)
    """
    bwa mem $ref $pairs -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${RG_id}\\tPL:ILLUMINA" > ${sample_id}_mito_PE.sam
    """
}

aln_ch = aln_SE.combine(aln_PE, by:0)

process SamtoolsSortMerge{
        
    tag "$sample_id"
    
    input:
    tuple val(sample_id), file(aln_SE), file(aln_PE) from aln_ch
    
    output:
    tuple val(merge_id), file("${sample_id}_mito_sorted.bam") into sam_ch
    file("${sample_id}_mito_sorted.bam.bai")

    script:
    merge_id = ("${sample_id}" - lib)
    """
    samtools view -bS $aln_SE | samtools sort -@ ${task.cpus} -O bam -o ${sample_id}_mito_sorted_SE.bam 
    samtools view -bS $aln_PE | samtools sort -@ ${task.cpus} -O bam -o ${sample_id}_mito_sorted_PE.bam 
    samtools merge -@ ${task.cpus} ${sample_id}_mito_sorted.bam ${sample_id}_mito_sorted_PE.bam ${sample_id}_mito_sorted_SE.bam 
    samtools index ${sample_id}_mito_sorted.bam
    """
}

//groupTuple for library merging
sam_all = sam_ch.groupTuple(by: 0, sort: true, size: 4)

process LibMerge{

    publishDir "${params.outdir}/4.MappedMergedReads/$sample_id", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(bams) from sam_all

    output:
    file("*")
    tuple val(sample_id), file("${sample_id}_mito_sorted.bam") into mapped_ch

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}_mito_sorted.bam $bams 
    samtools index ${sample_id}_mito_sorted.bam
    """
}

free_ch = ref2.combine(mapped_ch, by:0)

process VariantCall{

    publishDir "${params.outdir}/5.Variantcall/$sample_id", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(ref), file(bams) from free_ch

    output:
    tuple val(sample_id), file("${sample_id}.vcf.gz") into vcf_ch
    
    script:
    """
    freebayes -f $ref $bams --ploidy 1 | vcffilter -f "QUAL > 20" > ${sample_id}.vcf
    bgzip -c ${sample_id}.vcf > ${sample_id}.vcf.gz
    """
}


cons_ch = ref3.combine( vcf_ch, by: 0)

process CallConsensus {

    publishDir "${params.outdir}/6.consensus/$sample_id", mode:'copy'

    input:
    tuple val(sample_id), file(MB), file(vcf) from cons_ch

    output:
    file("*")

    script:
    """
    bcftools index $vcf
    bcftools consensus -f $MB -o ${sample_id}_consensus_mito.fa -s ${sample_id} $vcf
    sed -i "s/>.*/>${sample_id}/" ${sample_id}_consensus_mito.fa
    """
}
