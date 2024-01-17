//TODO: Add container
//TODO: add version samtools
process extract_read_names {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.txt"), emit: names

    """
	samtools view ${query} | cut -f1 > ${query.simpleName}.txt
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
//TODO: add version seqtk
process extract_reads {
    tag{names.simpleName}

    input:
    tuple path(names), path(reads)

    output:
    path("${names.simpleName}.extracted.fastq"), emit: reads

    """
	seqtk subseq ${reads} ${names} > ${names.simpleName}.extracted.fastq 
	"""
}