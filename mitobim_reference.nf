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
            Extract mitogenomes
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            nextflow run mitobim_reference.nf --reads_MB /PATH/TO/Indivxxx_L001_U.fastq.gz --outdir /PATH/TO/RESULTS
         
            'Mandatory arguments:'
            --reads_MB               Merged cleanded reads. i.e Indivxxx_L001_U.fastq.gz
            --outdir                 Path to output directory
            
            'OPTIONS'
            --help                   Outputs this help log      
            -resume                  Nextflow cmd to resume modified workflow
            --subset_size   INTEGER


            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 


            For mitobim see https://github.com/chrishah/MITObim
            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         Extract mitogenomes 
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         reads_MB     : ${params.reads_MB}
         outdir       : ${params.outdir}
         subset_size  : ${params.subset_size}
         """
         .stripIndent()


// Channel for mitobim reads
   Channel.fromPath( params.reads_MB)
         .ifEmpty { error "Cannot find any path matching: ${params.reads_MB}" }
         .map { it -> [it.name -  ~/\.fastq\.gz/, it] }
         .set { reads_MB_ch }

process Subsample {

	tag "$sample_id"

    publishDir "${params.outdir}/1.SubsetReadsMB/$sample_id", mode:'copy'

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

    publishDir "${params.outdir}/2.MitobimRef/$sample_id", mode:'copy'
    
    input:
    tuple val(sample_id), file(sub) from subset_ch

    output:
    tuple val(sample_id), file("manifest.conf") into mani_ch

    script:
    """
    manifest.sh $sample_id $params.mitobimRef $params.ref_strain "${params.outdir}/1.SubsetReadsMB/$sample_id/${sub}" manifest.conf
    """
}

subset_comb.combine( mani_ch, by: 0 )
    .set { reads_comb } 

process mira {
	
	tag "$sample_id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/2.MitobimRef/$sample_id", mode:'copy'
    
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

process mitobim {
	
	tag "$sample_id"
    
    label 'HIGH_RAM'

    publishDir "${params.outdir}/3.IndexRefs/$sample_id", mode:'copy'
    publishDir "${params.outdir}/2.MitobimRef/$sample_id", pattern: '*.log', mode:'copy'
    input:
    tuple val(sample_id), file("${sample_id}_subset.fastq") , file("${sample_id}.maf") from reads_comb_mira
    
    output:
    file("${sample_id}.log")
    file("${sample_id}-lycPyr_mtDNA-it*_noIUPAC.fasta")

    script:
    """
    $params.mitobim -start 1 -end 100 -sample $sample_id -ref $params.ref_strain -readpool "${sample_id}_subset.fastq" -maf "${sample_id}.maf" --clean --mirapath $params.mira_dir &> ${sample_id}.log

    cp ./iteration*/*_noIUPAC.fasta ./
    """
}



