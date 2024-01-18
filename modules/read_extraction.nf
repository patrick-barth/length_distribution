//TODO: Add container
process extract_read_names {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.txt"), emit: names
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
    path("${task.process}.version.txt"), 	emit: version

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