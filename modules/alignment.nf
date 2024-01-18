//TODO: add container
//TODO: add version seqret
process gb_to_fasta {
	
	input:
	path(ref)

	output:
	path("${ref.simpleName}.fasta"), emit: reference

	"""
	seqret -sequence ${ref} -outseq ${ref.simpleName}.fasta -osformat fasta
	"""
}

process build_index_bowtie {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*"), emit: index
	path("${task.process}.version.txt"), 	emit: version

	"""
	bowtie2-build ${ref} ${ref}
	
	echo -e "${task.process}\tbowtie2\t\$(bowtie2-build --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process mapping_bowtie{
	tag {query.simpleName}
	publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "${query.simpleName}.bam"
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query.simpleName}.statistics.txt"

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path "${query.simpleName}.bam", 			emit: bam_alignments
	path "${query.simpleName}.statistics.txt", 	emit: report
	path("${task.process}.version.txt"), 		emit: version

	script:
	def all_alignments = params.report_all_alignments ? '-a' : ''
	def some_alignments = params.max_alignments && !params.report_all_alignments ? "-k " + params.max_alignments : ''

	"""
	bowtie2 --no-unal \
		--very-sensitive \
		-L 10 \
		-q \
		${all_alignments} \
		${some_alignments} \
		-p ${task.cpus} \
		--seed 0 \
		-U ${query} \
		-x ${ref} \
		2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.simpleName}.bam

	echo -e "${task.process}\tbowtie2\t\$(bowtie2 --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

process build_index_STAR {

	input:
	path(referenceGenome)
	path(gtf)

	output:
	path(index), 							emit: index
	path("${task.process}.version.txt"), 	emit: version


	script:
	annotation_file = !params.annotation == 'NO_FILE' ? '--sjdbGTFfile ' + ${gtf} : ''

	"""
	mkdir index
	STAR --runThreadN ${task.cpus} \
		--runMode genomeGenerate \
		--genomeDir ./index \
		--genomeFastaFiles ${referenceGenome} \
		${annotation_file}
	
	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}

process mapping_STAR{
	tag {query.simpleName}
	publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "${query.baseName}.Aligned.sortedByCoord.out.bam"
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query.baseName}.Log.*"

	input:
	tuple path(query), path(indexDir)

	output:
	path("${query.baseName}.Aligned.sortedByCoord.out.bam"), 	emit: bam_alignments
	path("${query.baseName}.Log.*"), 							emit: report
	path("${task.process}.version.txt"), 						emit: version

	script:
	def all_alignments = params.report_all_alignments ? '--outSAMmultNmax -1' : ''
	def some_alignments = params.max_alignments && !params.report_all_alignments ? "--outSAMmultNmax " + params.max_alignments : ''

	"""
	STAR --runThreadN ${task.cpus} \
		--genomeDir ${indexDir} \
		--readFilesIn ${query} \
		--outFileNamePrefix ${query.baseName}. \
		${all_alignments} \
		${some_alignments} \
		--outSAMtype BAM SortedByCoordinate

	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}