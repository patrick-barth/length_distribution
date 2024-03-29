//TODO: Add container
//TODO: Add version
process filter_for_length {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.simpleName}.filtered.fastq"), emit: reads
	path("${task.process}.version.txt"), 		emit: version


	"""
	cutadapt -m ${params.nucleotide_length} -M ${params.nucleotide_length} -o ${query.simpleName}.filtered.fastq ${query}

	echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" > ${task.process}.version.txt
	"""
}

//TODO: Add container
process extract_read_names {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.txt"), 		emit: names
    path("${task.process}.version.txt"), 	emit: version

    """
	samtools view ${query} | cut -f1 > ${query.simpleName}.txt

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

process collect_reads {
    input:
    path(query)

    output:
    path('collected_reads.fastq'), emit: reads

    """
	cat ${query} > collected_reads.fastq
	"""
}

//TODO: Add container
process extract_reads {
    tag{names.simpleName}

    input:
    tuple path(names), path(reads)

    output:
    path("${names.simpleName}.extracted.fastq"), emit: reads
    path("${task.process}.version.txt"), 	emit: version

    """
	seqtk subseq ${reads} ${names} > ${names.simpleName}.extracted.fastq 
    
    echo -e "${task.process}\seqtk\t\$(seqtk 2>&1 | head -3 | tail -1 | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}

//TODO: Add container
//TODO: Add version
process extract_sequences_only {
	tag {query.simpleName}

	input:
	tuple val(origin), path(query)

	output:
	path("${query.simpleName}.${origin}.txt"), emit: sequences
	path("${task.process}.version.txt"), 	emit: version


	"""
	awk '{if(NR%4==2) print \$1}' ${query} > ${query.simpleName}.${origin}.txt

	echo -e "${task.process}\tawk\t\$(awk -W version | head -1)" > ${task.process}.version.txt
	"""
}

//TODO: Add container
//TODO: Add version
process calculate_nucleotide_distribution{
	tag {name}

	publishDir "${params.output_dir}/nucleotide_distribution", mode: 'copy'

	input:
	tuple val(name), path(sequences)

	output:
	path("${name}.${params.nucleotide_length}.nuc-distribution.tsv"), emit: distribution
	path("${task.process}.version.txt"), 	emit: version

	"""
	calc-nucleotide-distribution.py --sequences_all ${name}.all.txt --sequences_alignment ${name}.alignments.txt --length ${params.nucleotide_length} --output ${name}.${params.nucleotide_length}.nuc-distribution.tsv

	echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}