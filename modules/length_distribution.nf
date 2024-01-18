process count_length_distribution {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.length-distribution.txt"), emit: distribution
    path("${task.process}.version.txt"), 	emit: version

    """
	awk '{if(NR%4==2) print length(\$1)}' ${query} | sort -n | uniq -c > ${query.simpleName}.length-distribution.txt

    echo -e "${task.process}\tawk\t\$(awk -W version | head -1)" > ${task.process}.version.txt
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
    path("${task.process}.version.txt"), 	emit: version

    script:
	"""
	percentage-for-length-distribution.py --input ${query} --output ${query.simpleName}.perc.txt

    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}