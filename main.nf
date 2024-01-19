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

if(params.nucleotide_distribution){
    include{
        filter_for_length
        //extract_read_names
        //collect_reads
        //extract_reads
        extract_sequences_only
        calculate_nucleotide_distribution
    } from './modules/nucleotide_distribution.nf'
}

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: Pipeline that reports the length distribution of all reads alignablke to a provided reference.
                |       
                |Required arguments:
                |   --reads         Location of the input file file (FASTQ).
                |   --reference     Reference file. Can be provided as .fasta or .gb. If provided as .gb the 
                |                   reference will also be used as annotation file for feature splitting
                |
                |Optional arguments:
                |   --annotation    Annotation file. Will be used for alignment and to split read extraction
                |                   according to features. Can be provided as .gtf or .gff
                |   --min_length    Minimum length for reads after adapter trimming.
                |                   [default: ${params.min_length}]
                |   --aligner       States which alignment tool is used. Currently available are: 
                |                   'bowtie2' and 'star'
                |                   [default: ${params.aligner}]
                |   --report_all_alignments       Reports all potantial alignments per read. Overwrites --max_alignments
                |                                 [default: ${params.report_all_alignments}]
                |   --max_alignments        Maximum number of alignments returned per read.
                |                           [default: Default value of chosen aligner (--aligner)]
                |   --split_features        Split reads according to the feature they aligned to. Will create one output
                |                           per feature (given that at least one read aligned to the feature) additional
                |                           to the output for the complete sample
                |                           [default: ${params.split_features}]
                |   --filter_bacterial_contamination    Filters reads for bacterial contaminations.
                |                                       [default: ${params.filter_bacterial_contamination}]
                |   --kraken_db_dir     Database used to filter out bacterial contamination
                |                       [default: ${params.kraken_db_dir}]
                |                       Default only works with access to the computational system of the computational resources
                |                       of the Bioinformatics & Systems Biology group of the Justus Liebig University in Giessen:
                |                       https://www.uni-giessen.de/de/fbz/fb08/Inst/bioinformatik
                |   --nucleotide_distribution   Calculates nucleotide distribution for every position of reads and alignments of
                |                               a specified length
                |                               [default: ${params.nucleotide_distribution}]
                |   --nucleotide_length         Read/Alignment lengths to be investigated for nucleotide distribution
                |                               [default: ${params.nucleotide_length}]
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
        Aligner                         : ${params.aligner}
        Report all alignments           : ${params.report_all_alignments}
        Maximum alignments per read     : ${params.max_alignments}
        Split features                  : ${params.split_features}
        Nucleotide distribution         : ${params.nucleotide_distribution}
        Nucleotide length               : ${params.nucleotide_length}
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
            version_index           =   build_index_bowtie.out.version
            version_align           =   mapping_bowtie.out.version
            report_tmp              =   mapping_bowtie.out.report
        } else if (params.aligner == "star"){
            build_index_STAR(reference,
                            annotation)
            mapping_STAR(reads
                        .combine(build_index_STAR.out.index))

            alignments_tmp          =   mapping_STAR.out.bam_alignments
            version_index           =   build_index_STAR.out.version
            version_align           =   mapping_STAR.out.version
            report_tmp              =   mapping_STAR.out.report
        } 
        //collect versions
        versions = version_index.first()
                    .concat(version_align.first())

    emit:
        versions        =   versions
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

        //Collect versions
        versions = count_features.out.version.first()
                    .concat(feature_splitting.out.version.first())
        versions = reference_extension == 'gb' ? versions.concat(gb_to_gtf.out.version.first()) : versions

    emit:
        read_names_split    =   feature_splitting.out.read_names
        versions            =   versions
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

        //Collect versions
        versions = extract_read_names.out.version.first()
                    .concat(extract_reads.out.version.first())

    emit:
        extracted_reads = extract_reads.out.reads
        versions        = versions
}

workflow length_distribution{
    take:
        reads
    main:
        count_length_distribution(reads)
        calculate_length_percentage(count_length_distribution.out.distribution)

        //Collect versions
        versions = count_length_distribution.out.version.first()
                    .concat(calculate_length_percentage.out.version.first())

    emit:
        length_percentage   = calculate_length_percentage.out.percentage
        versions            = versions
}

workflow nucleotide_distribution{
    take:
        processed_reads
        alignments
        //lengths //TODO: needs to be implemented in away that several lengths parameters are allowed
    main:
        filter_for_length(processed_reads)
        extract_read_names(alignments)
        collect_reads(filter_for_length.out.reads)
        extract_reads(extract_read_names.out.names.flatten()
                        .combine(collect_reads.out.reads))
        // Create new channels with read origins to help distinguish them
        Channel.of('all')
            .combine(filter_for_length.out.reads)
            .set{all_reads}
        Channel.of('alignments')
            .combine(extract_reads.out.reads)
            .set{alignable_reads}
        extract_sequences_only(all_reads.concat(alignable_reads))
        calculate_nucleotide_distribution(extract_sequences_only.out.sequences
                                            .map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/, file)}
                                            .groupTuple())

        //Collect versions
        versions = filter_for_length.out.version.first()
                    .concat(extract_read_names.out.version.first())
                    .concat(extract_reads.out.version.first())
                    .concat(extract_sequences_only.out.version.first())
                    .concat(calculate_nucleotide_distribution.out.version.first())

    emit:
        distribution    = calculate_nucleotide_distribution.out.distribution
        versions        = versions
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
    if(params.nucleotide_distribution){
        nucleotide_distribution(preprocessing.out.fastq_reads,
            alignment.out.alignments)
    }

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
                        .concat(alignment.out.versions)
                        .concat(read_extraction.out.versions)
                        .concat(length_distribution.out.versions)
                        .concat(collect_metadata.out.version)
                        .concat(get_md5sum.out.version)

    collected_versions = params.split_features ? collected_versions.concat(count_split_features.out.versions) : collected_versions
    collected_versions = params.nucleotide_distribution ? collect_versions.concat(nucleotide_distribution.out.versions) : collected_versions

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