//TODO: add version awk
process count_length_distribution {
    tag{query.simepleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.length-distribution.txt"), emit: distribution

    """
	awk '{if(NR%4==2) print length(\$1)}' ${query} | sort -n | uniq -c > ${query.simpleName}.length-distribution.txt
	"""
}

//TODO: Add container
//TODO: add version python
process calculate_length_percentage {
    tag{query.simpleName}
    publishDir "${params.output_dir}/length-distribution", mode: 'copy', pattern: "${query.simpleName}.perc.txt"

    input:
    path(query)

    output:
    path("${query.simpleName}.perc.txt"), emit: percentage

    script:
	"""
	percentage-for-length-distribution.py --input ${query} --output ${query.simpleName}.perc.txt
	"""
}