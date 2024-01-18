//TODO: Create Docker container
process gb_to_gtf {
    
    input:
    path(ref)

    output:
    path("${ref.simpleName}.gtf"), emit: gtf
    path("${task.process}.version.txt"), 	emit: version

    """
    seqret -sequence ${ref} -outseq ${ref.simpleName}.gff -osformat gff3 -feature
	sed -i 's/ID=.*;label=\\(.*\\);/ID=\\1;label=\\1;/' ${ref.simpleName}.gff
	# remove fasta from gff
	sed -i -n '/##FASTA/q;p' ${ref.simpleName}.gff

    # convert gff to gtf
	gffread ${ref.simpleName}.gff -o ${ref.simpleName}.gtf -T
	sed -i 's/transcript_id/gene_id/g' ${ref.simpleName}.gtf

    echo -e "${task.process}\tseqret\t\$(seqret --version 2>&1)" > ${task.process}.version.txt
    echo -e "${task.process}\tgffread\t\$(gffread --version 2>&1)" >> ${task.process}.version.txt
    """
}

//TODO: Add container
process count_features {
    tag{query.simpleName}

    input:
    tuple path(query), path(annotation)

    output:
    path("${query}.featureCounts.bam"), emit: feature_alignments
    path("${task.process}.version.txt"), 	emit: version

    """
    featureCounts -a ${annotation} -o featureCounts.tsv ${query} -R BAM

    echo -e "${task.process}\tfeatureCounts\t\$(featureCounts -v 2>&1 | head -2 | tail -1 | cut -d' ' -f2)" >> ${task.process}.version.txt
    """
}

//TODO: Add container
process feature_splitting {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("*.txt"), emit: read_names
    path("${task.process}.version.txt"), 	emit: version

    """
	#samtools view ${query} | awk 'BEGIN {FS=OFS="\\t"} {if(!/^@/ && !/^\$/){\$3 = \$22}; print}' | sed 's/XT:Z://g' | samtools view -b -o ${query.baseName}_feature.bam -
	# This line splits read names according to the feature the reads were aligned to
	samtools view ${query} | cut -f1,22 | sed 's/XT:Z://g' | awk  'BEGIN {FS=OFS="\\t"} {sub("\$", ".txt", \$2); \$2 = "${query.simpleName}_"\$2; print}' | awk '{print \$1 > \$2}'

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
    echo -e "${task.process}\tawk\t\$(awk -W version | head -1)" > ${task.process}.version.txt
	"""
}