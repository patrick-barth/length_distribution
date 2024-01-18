process quality_control {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics/qc-preprocessing", mode: 'copy', pattern: "${query.baseName}_fastqc.{html,zip}"

	
	input:
	path query

	output:
	path "${query.baseName}_fastqc.{html,zip}", 				emit: output
	path("${task.process}.version.txt"), 	emit: version

	"""
	fastqc ${query} -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process quality_control_2 {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics/qc-postprocessing", mode: 'copy', pattern: "${query.simpleName}_2_fastqc.{html,zip}"

	
	input:
	path query

	output:
	path "${query.simpleName}_2_fastqc.{html,zip}", emit: output
	path("${task.process}.version.txt"), 			emit: version

	"""
	cat ${query} > ${query.simpleName}_2.fastq
	fastqc ${query.simpleName}_2.fastq -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
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

//TODO: Create Docker container
//TODO: Check which reports are actually important
//TODO: Get version of bracken
process filter_bacterial_contamination {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics/bacterial_contamination_filter", mode: 'copy', pattern: "${query.simpleName}.bac-filter.report"

	input:
	file(query)
	path(db)

	output:
	path("${query.simpleName}.bac-filter.report"), emit: report
	path("${query.simpleName}.bac-filtered.fastq"), emit: fastq
	path("${task.process}.version.txt"), 	emit: version

	"""
	kraken2	--use-names \
		--threads ${task.cpus} \
		--db ${db} \
		--report ${query.simpleName}.bac-filter.report \
		--unclassified-out ${query.simpleName}.bac-filtered.fastq \
		${query} \
		> ${query.simpleName}.kraken

	echo -e "${task.process}\tkraken2\t\$(kraken2 --version | head -1 | cut -d' ' -f3)" > ${task.process}.version.txt
	"""
}