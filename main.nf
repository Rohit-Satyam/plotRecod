#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if( params.help ) {

log.info """
* -------------------------------------------------
 *  RNAgrinder@KAUST: Analyzing RNASeq Dataset
 * -------------------------------------------------
Usage:
	nextflow run main.nf --input "${params.input}" --outdir ${params.outdir} --ref ${params.ref} \
	--feature ${params.gtf} --mode ${params.mode} 
Input:
	#### Mandatory Arguments ####
	* --cpus: No of threads to run the pipeline. Default [${params.cpus}]
	* --feature: Path to snapGene Feature file. Default [${params.feature}]
	* --input: Path to FastqQ files. Default [${params.input}]
	* --mode: If data is Paired-end pass "PE" else "SE". Default [${params.mode}]
	* --outdir: Path/Name of the output directory. Default [${params.outdir}]
	* --ref: Path to reference fasta file. Default [${params.ref}]
	* --geneID: Provide the gene ID of the recodonised genes. This will be used as header of reference fasta and in naming feature file. Default [${params.geneID}]
	
	Default [${params.index_dir}]

	#### Parameters to pass additional Arguments ####
	* --bwa_ext: Additional arguments to pass to BWA. Default [${params.bwa_ext}]
	* --fastp_ext: Additional arguments to pass to FASTP. Default [${params.fastp_ext}]
	* --samtools_ext: Additional arguments to pass to Samtools. Default [${params.samtools_ext}]

	#### Parameters to Skip certain Steps ####
	* --skipTrim: Set this "true" to skip Trimming Step. Default [${params.skipTrim}]
"""

exit 0
}

include {FASTP} from './modules/01_fastp'
include {BWA2;modifyFeature} from './modules/02_align'
include {BAMSTATS} from './modules/03_bamstats'

params.help= false
params.input = false
params.outdir= false

params.mode = false


workflow{

if (params.input != false){
	if (params.mode == "PE"){
		Channel.fromFilePairs(params.input, checkIfExists: true )
		.set { input_fastqs }
		} else if (params.mode == "SE") {
			Channel.fromPath(params.input, checkIfExists: true ).map { file -> tuple(file.simpleName, file)}
			.set { input_fastqs }
	}
}


if (params.skipTrim == false){
        FASTP(input_fastqs)
	BWA2(FASTP.out[0],Channel.fromPath(params.ref))
	BAMSTATS(BWA2.out[0])
} else {
	BWA2(input_fastqs,Channel.fromPath(params.ref))
	BAMSTATS(BWA2.out[0])
}

if (params.feature !=""){
modifyFeature(Channel.fromPath(params.feature))
}


}
