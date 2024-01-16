#!/usr/bin/env nextflow

import groovy.json.JsonOutput // used for parameter output

nextflow.enable.dsl=2

include{
    collect_metadata
    get_md5sum
    multiqc
    collect_versions
} from './modules/default_processes.nf'

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: A description of your script and maybe some examples of how
                |                to run the script
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
         input from   : ${params.input_file}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         """
         .stripIndent()

//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

// Collect all input files
input_files = input_reads.concat(Channel.of(annotation))
                    .concat(reference)
                    .flatten().toList()

/*
 * Starting subworkflow descriptions
 */
workflow preprocessing {
    take: 
        input_reads
    main:
        quality_control(input_reads)
        adapter_removal(input_reads)
        quality_filter(adapter_removal.out.fastq_trimmed)
        quality_control_2(quality_filter.out.fastq_quality_filtered)

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out
        multiqc_quality_control_post_preprocessing  = quality_control_2.out
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_quality_filter                      = quality_filter.out.report_quality_filter

        // data for downstream processes
        fastq_reads_quality_filtered                = quality_filter.out.fastq_quality_filtered
}



/*
 * Actual workflow connecting subworkflows
 */
workflow {
    preprocessing(input_reads)

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