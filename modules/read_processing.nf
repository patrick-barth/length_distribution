/*
 * Checks reads for general metrics
 * Input: [FASTQ] Unpreprocessed reads 
 * Output: [HTML] General report for unpreprocessed reads 
 */
process quality_control {
	tag {query.simpleName}
	
	input:
	path query

	output:
	path "${query.baseName}*"

	"""
	fastqc ${query} -o .
	"""
}

/*
 * Checks preprocessed reads for general metrics
 * Input: [FASTQ] Preprocessed reads 
 * Output: [HTML] General report for preprocessed reads
 */
process quality_control_2 {
	tag {query.simpleName}
	
	input:
	path query

	output:
	path "quality-control-2*"

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .
	"""
}

/*
 * Removes adapters from reads
 * Input: [FASTQ] Read file
 * Params: params.min_length -> Minimum length reads need after trimming to not be omitted
 * Output:  fastq_trimmed 	-> [FASTQ] Read file with afdapters trimmed
 *			report_trimming -> [TXT] Report adapter trimming
 */
process adapter_removal {
	tag {query.simpleName}

	input:
	path query

	output:
	path "${query}_trimmed.fq", emit: fastq_trimmed 
	path "${query}_trimming_report.txt", emit: report_trimming 
	tuple path("${task.process}.version.txt"), path("${task.process}.version2.txt"), 	emit: version

	"""
	trim_galore --cores ${task.cpus} --basename ${query} -o . --length ${params.min_length} ${query} --quality 0

	echo -e "${task.process}\ttrim_galore\t\$(trim_galore -v | head -4 | tail -1 | sed -e 's/^[ \t]*//' | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" > ${task.process}.version2.txt
	"""
}

process filter_bacterial_contamination {
	tag {query.simpleName}
	conda 'kraken2 bracken' //TODO: Create Docker container
	publishDir "${params.output_dir}/statistics/bacterial_contamination_filter", mode: 'copy', pattern: "${query.simpleName}.bac-filter.report"

	input:
	file(query)

	output:
	path("${query.simpleName}.bac-filter.report"), emit: report
	path("${query.simpleName}.bac-filtered.fastq"), emit: fastq

	"""
	kraken2	--use-names \
		--threads ${task.cpus} \
		--db ${params.kraken_db} \
		--fastq-input \
		--report ${query.simpleName}.bac-filter.report \
		--unclassified-out ${query.simpleName}.bac-filtered.fastq \
		${query} \
		> ${query.simpleName}.kraken //TODO: Check which reports are actually important
	"""
}