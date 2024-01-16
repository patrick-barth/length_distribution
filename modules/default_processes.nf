process collect_metadata{
    output:
    path("metadata.txt"),                   emit: output
    path("${task.process}.version.txt"),    emit: version

    script:
    """
    cat <<EOF > metadata.txt
    Author: ${params.manifest.author}
    Pipeline version: ${params.manifest.version}
    Working directory: ${workflow.workDir}
    User name: ${workflow.userName}
    Command line: ${workflow.commandLine}
    Container engine: ${workflow.containerEngine}    
    Containers used: ${workflow.container}
    Git repository: ${workflow.repository}
    Repository revision: ${workflow.revision}
    EOF

    echo -e "${task.process}\tcat\t\$(cat --version | head -1 | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
    """
}

process get_md5sum {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "md5sums.txt"

    input:
    path(query)

    output:
    path("md5sums.txt"),                    emit: metadata_output
    path("${task.process}.version.txt"),    emit: version

    script:
    """
    for i in ${query}
	do
		if test -f \$i; then
			md5sum \$i >> md5sums.txt
		fi
	done

    echo -e "${task.process}\tmd5sum\t\$(md5sum --version | head -1 | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
    """
}

process multiqc{
	publishDir "${params.output_dir}/statistics", mode: 'copy'

	input:
	path(adapter)
	path(qual)
	path(qc1)
	path(qc2)
	path(mapping)

	output:
	path("multiqc_*"),                      emit: output
    path("${task.process}.version.txt"),    emit: version
    
    script:
	"""
	multiqc .

    
    echo -e "${task.process}\tmultiqc\t\$(multiqc --version |  cut -f3 -d' ')" > ${task.process}.version.txt
	"""
}

process collect_versions {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "tool_versions.txt"

    input:
    path(query)

    output:
    path('tool_versions.txt')

    script:
    """
    echo -e "Process\ttool\tversion" > tool_versions.txt
    for i in ${query}
	do
		cat \$i >> tool_versions.txt
	done
    """
}
