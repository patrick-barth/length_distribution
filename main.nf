#!/usr/bin/env nextflow

import groovy.json.JsonOutput // used for parameter output

nextflow.enable.dsl=2

include{
    collect_metadata
    get_md5sum
    multiqc
    collect_versions
} from './modules/default_processes.nf'

include{
    quality_control
    quality_control_2
    adapter_removal
    filter_bacterial_contamination
} from './modules/read_processing.nf'

include{
    gb_to_fasta
} from './modules/alignment.nf'

if (params.aligner == "bowtie2"){
    include{
        build_index_bowtie
        mapping_bowtie
    } from './modules/alignment.nf'
} else if (params.aligner == "star"){
    include{
        build_index_STAR
        mapping_STAR
    } from './modules/alignment.nf'
}


/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: Pipeline that reports the length distribution of all reads alignablke to a provided reference.
                |       
                |Required arguments:
                |   --reads         Location of the input file file (FASTQ).
                |
                |Optional arguments:
                |   --min_length    Minimum length for reads after adapter trimming.
                |                   [default: ${params.min_length}]
                |   --min_qual      Minimum base quality.
                |                   [default: ${params.min_qual}]
                |   --min_percent_qual_filter   Minimum percentage of bases within a read that need to
                |                               be above the quality threshold
                |                               [default: ${params.min_percent_qual_filter}]
                |   --aligner       States which alignment tool is used. Currently available are: 
                |                   'bowtie2' and 'star'
                |  -w            The NextFlow work directory. Delete the directory once the process
                |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//preparation for workflow

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        Reads        : ${params.input_file}
        Reference    : ${params.reference}
        Annotation   : ${params.annotation}
        output to    : ${params.output_dir}
        --
        Minmum read length              : ${params.min_length}
        Filter bacterial contamination  : ${params.filter_bacterial_contamination}
        Kraken DB directory             : ${params.kraken_db}
        --
        run as       : ${workflow.commandLine}
        started at   : ${workflow.start}
        config files : ${workflow.configFiles}
        """
        .stripIndent()

/*
 * Input
 */

//essential input files
input_reads     = Channel.fromPath( params.reads )
reference       = Channel.fromPath( params.reference )

//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

// Collect all input files
input_files = input_reads.concat(Channel.of(annotation))
                    .concat(reference)
                    .flatten().toList()

reference_extension = reference.getExtension()

/*
 * Starting subworkflow descriptions
 */
workflow preprocessing {
    take: 
        input_reads
    main:
        quality_control(input_reads)
        adapter_removal(input_reads)
        if(filter_bacterial_contamination) { 
            filter_bacterial_contamination(adapter_removal.out.fastq_trimmed)
        }
        processed_reads = filter_bacterial_contamination ? filter_bacterial_contamination.out.fastq_filtered : adapter_removal.out.fastq_trimmed
        quality_control_2(processed_reads)

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out
        multiqc_quality_control_post_preprocessing  = quality_control_2.out
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming

        // data for downstream processes
        fastq_reads_quality_filtered                = processed_reads
}

workflow alignment {
    take:
        reference
        annotation
        reads

    main:
        if(params.aligner == "bowtie2"){
            build_index_bowtie(reference)
            mapping_bowtie(build_index_bowtie.out.index.first(),
                            reads)

            alignments_tmp          =   mapping_bowtie.out.bam_alignments
            version_index_tmp       =   build_index_bowtie.out.version
            version_align_tmp       =   mapping_bowtie.out.version
            report_tmp              =   mapping_bowtie.out.report
        } else if (params.aligner == "star"){
            build_index_STAR(reference,
                            annotation)
            mapping_STAR(reads
                        .combine(build_index_STAR.out.index))

            alignments_tmp          =   mapping_STAR.out.bam_alignments
            version_index_tmp       =   build_index_STAR.out.version
            version_align_tmp       =   mapping_STAR.out.version
            report_tmp              =   mapping_STAR.out.report
        } 

    emit:
        version_index   =   version_index_tmp
        version_align   =   version_align_tmp
        reports         =   report_tmp

        alignments      =   alignments_tmp
}

workflow count_features{
    take:
        //TODO: Inputs
    main:
        //TODO: Logic
    emit:
        //TODO: Outputs
}

workflow read_extraction{
    take:
        //TODO: Inputs
    main:
        //TODO: Logic
    emit:
        //TODO: Outputs
}

workflow length_distribution{
    take:
        //TODO: Inputs
    main:
        //TODO: Logic
    emit:
        //TODO: Outputs
}

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    preprocessing(input_reads)
    

    //Alignment
    if(reference_extension == 'gb'){
        gb_to_fasta()//TODO:input from reference as gb
    }
    reference_collect = reference_extension == 'gb' ? gb_to_fasta.out.reference : reference
    alignment() //TODO: input from preprocessing
    //feature_splitting
    if(reference_extension == 'gb'){
        gb_to_gtf()//TODO: input from reference as gb
    }
    if(){ // split after features
        count_features() //TODO: input from 
    }
    read_extraction() //TODO: Input from alignments and processed reads


    // Collect metadata
    collect_metadata()
    get_md5sum(input_files)
    collect_versions(collect_metadata.out.version
                        .concat(get_md5sum.out.version)
                        .unique()
                        .flatten().toList()
    )
}



/*
 * Prints complection status to command line
 */
workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}