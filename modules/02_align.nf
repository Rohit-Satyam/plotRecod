params.memory = "3g"
params.cpus = 1
params.outdir = "."
params.jobs = 1


process BWA2{
        publishDir "${params.outdir}/02_alignment/", mode: 'copy'
        cpus params.cpus
        maxForks params.jobs

input:
        tuple val(sid), path(reads)
        each path(ref)

output:
tuple val(sid), path ("${sid}.sorted.bam")
path("${sid}.sorted.bam.bai")
path("*.fa")

shell:
if ("${params.mode}" == "PE")

'''
## Read group info
id=$(zcat !{reads[0]} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//')

## Renaming the fasta header with Gene ID
seqkit seq !{ref} -u | seqtk seq -C - | seqtk rename - | sed 's/>1/>!{params.geneID}/g'  > !{params.geneID}.fa

# Index reference
bwa-mem2 index !{params.geneID}.fa

# Align the reads
bwa-mem2 mem -M -R "$(echo "@RG\\tID:${id}\\tSM:!{sid}\\tPL:ILLUMINA")" \
-t !{task.cpus} !{params.geneID}.fa !{params.bwa_ext}\
!{reads[0]} !{reads[1]} | samtools sort -@ !{task.cpus} !{params.samtools_ext} -o !{sid}.sorted.bam  -

samtools index !{sid}.sorted.bam
'''

}


process modifyFeature{
        publishDir "${params.outdir}/02_alignment/", mode: 'copy'
        cpus params.cpus
        maxForks params.jobs
input:
path(feature)

output:
path("*.tsv")


shell:


'''
# Preparing minimal feature file after removing duplicate regions
!{projectDir}/bin/reformat.sh !{feature.toRealPath()} > !{params.geneID}.feature.tsv
'''

}
