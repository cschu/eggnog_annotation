#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process eggnog_mapper {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(genome_id), path(proteins)
	path(eggnog_db)

	output:
	tuple val(genome_id), path("${genome_id}/${genome_id}.emapper.annotations"), emit: eggnog

	script:
	"""
	mkdir -p ${genome_id}
	gzip -dc ${genome_id}.faa.gz > ${genome_id}.faa
	emapper.py -i ${genome_id}.faa --data_dir ${eggnog_db} --output ${genome_id}/${genome_id} -m diamond --cpu $task.cpus --dbmem
	rm -f ${genome_id}.faa
	"""
}


process prodigal {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(genome_id), path(genome_fna)

	output:
	tuple val(genome_id), path("${genome_id}/${genome_id}.faa.gz"), emit: proteins
	tuple val(genome_id), path("${genome_id}/${genome_id}.ffn.gz"), emit: genes
	tuple val(genome_id), path("${genome_id}/${genome_id}.gff.gz"), emit: genome_annotation

	script:
	def gunzip_cmd = (genome_fna.name.endsWith(".gz")) ? "gzip -dc ${genome_fna} > \$(basename ${genome_fna} .gz)" : ""
	"""
	mkdir -p ${genome_id}
	${gunzip_cmd}
	prodigal -i \$(basename ${genome_fna} .gz) -f gff -o ${genome_id}/${genome_id}.gff -a ${genome_id}/${genome_id}.faa -d ${genome_id}/${genome_id}.ffn
	gzip ${genome_id}/*
	"""
}


workflow {
	genomes_ch = Channel
		.fromPath(params.input_dir + "/" + params.file_pattern)
		.map { file ->
			def genome_id = file.name.replaceAll(suffix_pattern, "")
			return tuple(genome_id, file)
		}
		.groupTuple(sort: true)

	/* Classify and annotate the input genomes and filter by annotation */

	prodigal(genomes_ch)

	eggnog_mapper(prodigal.out.proteins, params.eggnog_db)
}
