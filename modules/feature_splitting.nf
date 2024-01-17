//TODO: Create Docker container
//TODO: Add version seqret
//TODO: Add version gffread
//TODO: Add version sed
process gb_to_gtf {
    
    input:
    path(ref)

    output:
    path("${ref.simpleName}.gtf"), emit: gtf

    """
    seqret -sequence ${ref} -outseq ${ref.simpleName}.gff -osformat gff3 -feature
	sed -i 's/ID=.*;label=\\(.*\\);/ID=\\1;label=\\1;/' ${ref.simpleName}.gff
	# remove fasta from gff
	sed -i -n '/##FASTA/q;p' ${ref.simpleName}.gff

    # convert gff to gtf
	gffread ${ref.simpleName}.gff -o ${ref.simpleName}.gtf -T
	sed -i 's/transcript_id/gene_id/g' ${ref.simpleName}.gtf
    """
}

//TODO: Add container
//TODO: add version featureCounts
process count_features {
    tag{query.simpleName}

    input:
    tuple path(query), path(annotation)

    output:
    path("${query}.featureCounts.bam"), emit: feature_alignments

    """
    featureCounts -a ${annotation} -o featureCounts.tsv ${query} -R BAM
    """
}

//TODO: Add container
//TODO: add version samtools
process feature_splitting {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("*.txt"), emit: read_names

    """
	#samtools view ${query} | awk 'BEGIN {FS=OFS="\\t"} {if(!/^@/ && !/^\$/){\$3 = \$22}; print}' | sed 's/XT:Z://g' | samtools view -b -o ${query.baseName}_feature.bam -
	# This line splits read names according to the feature the reads were aligned to
	samtools view ${query} | cut -f1,22 | sed 's/XT:Z://g' | awk  'BEGIN {FS=OFS="\\t"} {sub("\$", ".txt", \$2); \$2 = "${query.simpleName}_"\$2; print}' | awk '{print \$1 > \$2}'
	"""
}