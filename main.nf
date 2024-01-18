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

include{
    gb_to_gtf
    count_features
    feature_splitting
} from './modules/feature_splitting.nf'

include{
    extract_read_names
    collect_reads
    extract_reads
} from './modules/read_extraction.nf'

include{
    count_length_distribution
    calculate_length_percentage
} from './modules/length_distribution.nf'


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
        Reads        : ${params.reads}
        Reference    : ${params.reference}
        Annotation   : ${params.annotation}
        output to    : ${params.output_dir}
        --
        Minmum read length              : ${params.min_length}
        Filter bacterial contamination  : ${params.filter_bacterial_contamination}
        Kraken DB directory             : ${params.kraken_db_dir}
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

kraken_db = file(params.kraken_db_dir).toAbsolutePath()

// Collect all input files
input_files = input_reads.concat(Channel.of(annotation))
                    .concat(reference)
                    .concat(Channel.of(kraken_db))
                    .flatten().toList()

reference_extension = file(params.reference).getExtension()

/*
 * Starting subworkflow descriptions
 */
workflow preprocessing {
    take: 
        input_reads
        kraken_db
    main:
        quality_control(input_reads)
        adapter_removal(input_reads)
        if(params.filter_bacterial_contamination) { 
            filter_bacterial_contamination(adapter_removal.out.fastq_trimmed,
                kraken_db)
        }
        processed_reads = params.filter_bacterial_contamination ? filter_bacterial_contamination.out.fastq : adapter_removal.out.fastq_trimmed
        quality_control_2(processed_reads)

        // Collect versions
        versions = quality_control.out.version.first()
                        .concat(quality_control_2.out.version.first())
                        .concat(adapter_removal.out.version.first())
        versions = params.filter_bacterial_contamination ? versions.concat(filter_bacterial_contamination.out.version.first()) : versions

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out.output
        multiqc_quality_control_post_preprocessing  = quality_control_2.out.output
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_bac_contamination                   = params.filter_bacterial_contamination ? filter_bacterial_contamination.out.report : Channel.empty()

        versions = versions

        // data for downstream processes
        fastq_reads                                 = processed_reads
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

workflow count_split_features{
    take:
        reference
        alignments
        annotation

    main:
        if(reference_extension == 'gb'){
            gb_to_gtf(reference)
        }
        feature_info = reference_extension == 'gb' ? gb_to_gtf.out.gtf : annotation
        count_features(alignments.flatten()
                        .combine(feature_info))
        feature_splitting(count_features.out.feature_alignments)

    emit:
        read_names_split    =   feature_splitting.out.read_names
}

workflow read_extraction{
    take:
        alignments
        reads
        split_features

    main:
        extract_read_names(alignments)
        collect_reads(reads.collect())
        
        extract_read_names.out.names
            .concat(split_features).
            set{collected_read_names}

        extract_reads(collected_read_names.flatten()
                        .combine(collect_reads.out.reads))

    emit:
        extracted_reads = extract_reads.out.reads
}

workflow length_distribution{
    take:
        reads
    main:
        count_length_distribution(reads)
        calculate_length_percentage(count_length_distribution.out.distribution)
    emit:
        length_percentage = calculate_length_percentage.out.percentage
}

/*
 * Actual workflow
 */
workflow {

    preprocessing(input_reads,
        kraken_db)
    
    //Alignment
    if(reference_extension == 'gb'){
        gb_to_fasta(reference)
    }
    reference_collect = reference_extension == 'gb' ? gb_to_fasta.out.reference : reference
    alignment(reference_collect,
        annotation,
        preprocessing.out.fastq_reads) 
    //feature_splitting
    if(params.split_features){
        count_split_features(reference,
            alignment.out.alignments,
            annotation) 
    }
    collect_split = params.split_features ? count_split_features.out.read_names_split : Channel.empty()
    read_extraction(alignment.out.alignments,
        preprocessing.out.fastq_reads,
        collect_split)

    length_distribution(read_extraction.out.extracted_reads)



    /*
     * Collect metadata
     */
    
    multiqc(preprocessing.out.multiqc_adapter_removal,
        preprocessing.out.multiqc_quality_control,
        preprocessing.out.multiqc_quality_control_post_preprocessing,
        alignment.out.reports,
        preprocessing.out.multiqc_bac_contamination
    )
    collect_metadata()
    get_md5sum(input_files)

    /*
    * Collect versions of all invoked processes
    */
    collected_versions = preprocessing.out.versions
                        .concat(alignment.out.version_index.first())
                        .concat(alignment.out.version_align.first())
                        .concat(extract_read_names.out.version.first())
                        .concat(extract_reads.out.version.first())
                        .concat(count_length_distribution.out.version.first())
                        .concat(calculate_length_percentage.out.version.first())
                        .concat(collect_metadata.out.version)
                        .concat(get_md5sum.out.version)

    collected_versions = reference_extension == 'gb' ? collected_versions.concat(gb_to_fasta.out.version.first()) : collected_versions
    collected_versions = reference_extension == 'gb' && params.split_features ? collected_versions.concat(gb_to_gtf.out.version.first()) : collected_versions
    collected_versions = params.split_features ? collected_versions.concat(count_features.out.version.first()).concat(feature_splitting.out.version.first()) : collected_versions

    collect_versions(collected_versions
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