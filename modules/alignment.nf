/*
 * Prepares index for bowtie2
 * Input: [FASTA] Reference file 
 * Output: Tuple of [FASTA] reference file and [BT2] index files generated for the reference   
 */
process build_index_bowtie {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*")

	"""
	bowtie2-build ${ref} ${ref}
	"""
}

/*
 * Alignes reads to a reference via bowtie2, filters out unaligned sequeces
 * and converts output to BAM. Reference needs to be indexed.
 * Input: Tuple of [FASTA] reference file and [BT2] index files generated for the reference
 * 		[FASTA] Read files to be aligned to the reference
 * Params:  params.report_all_alignments    -> Reports all possible alignments and dismisses params.max_alignments
 *          params.max_alignments           -> Amount of alignemnts to report. Only used when params.report_all_alignments is false
 * Output: bam_alignments 	-> [BAM] Aligned sequences
 * 		  report_alignments -> [TXT] Alignment reports
 */
process mapping_bowtie{
	tag {query.simpleName}

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path "${query.baseName}.bam", emit: bam_alignments
	path "${query.simpleName}.statistics.txt", emit: report_alignments

	script:
	if(params.report_all_alignments)
		"""
		bowtie2 --no-unal -q -a -p ${task.cpus} --seed 0 -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
		"""
    else if(params.max_alignments)
        """
        bowtie2 --no-unal -q -p ${task.cpus} --seed 0 -k ${params.max_alignments} -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
        """
	else
		"""
		bowtie2 --no-unal -q -p ${task.cpus} --seed 0 -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
		"""
}

/*
 * Prepares index for STAR. Includes alignment if one is given
 * Input: [FASTA] Reference sequence
 *		[GTF]|[GFF3] Annotation file - if none is given it should say [NO_FILE]
 * Output: [DIR] Directory with index generated for given reference file
 */
process build_index_STAR {

	input:
	path(referenceGenome)
	path(gtf)

	output:
	path(index)

	script:
	if(params.annotation == 'NO_FILE')
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} 
		"""
	else
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${gtf}
		"""
}

/*
 * Alignes reads to reference via STAR. Reference needs to be index before.
 * Suppresses clipping at 5' end (--alignEndsType Extend5pOfRead1)
 * Output is provided as BAM file that is sorted by coordinates
 * Input: Tuple of [FASTQ] Read files to be aligned and [DIR] Directory containing the STAR index
  * Params:  params.report_all_alignments    -> Reports all possible alignments and dismisses params.max_alignments
 *          params.max_alignments           -> Amount of alignemnts to report. Only used when params.report_all_alignments is false
 * Output: bam_alignments -> [BAM] Aligned reads sorted by coordinates 
 *		report_alignments -> [TXT] Alignment report
 */
process mapping_STAR{
	tag {query.simpleName}

	input:
	tuple path(query), path(indexDir)

	output:
	path("${query.baseName}.Aligned.sortedByCoord.out.bam"), emit: bam_alignments
	path("${query.baseName}.Log.*"), emit: report_alignments

	script:
	if(params.report_all_alignments)
		"""
		STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --outSAMmultNmax -1 --outSAMtype BAM SortedByCoordinate
		"""
    else if(params.max_alignments)
        """
        STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --outSAMmultNmax ${params.max_alignments} --outSAMtype BAM SortedByCoordinate
        """
	else
		"""
		STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --outSAMmultNmax ${params.max_alignments} --outSAMtype BAM SortedByCoordinate
		"""
}