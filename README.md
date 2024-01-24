![image](https://github.com/Rohit-Satyam/plotRecod/assets/55479480/75f3357b-1f47-4450-8292-b79e264b8821)# plotRecod
Pipeline to plot coverage of recodonized genes from RNAseq data.

## Input required
**For Nextflow Pipeline**
1. Fasta file from SNAPGene software of the recodonized gene. Since the pipeline uses BWA which does not perform spliced alignment, users are expected to provide the recodonized fasta with exons only. But feel free to add other regions such as loxPint and UTRs.
2. (optional) A feature file from SNAPgene containing information about the regions and their corresponding information. This file is used for visualization purposes only and is modified and renamed when given along with the nextflow command. Users are encouraged to remove unnecessary regions like primers, exons, mRNA, genes etc. which are redundant and span the same region. If not the pipeline tries to get rid of features having exact coordinates (start and end) and remove primer entries automatically. The `geneID.tsv` file so produced can be further modified before proceeding to the visualization step.

**For plotting coverage**
1. A `.csv` file containing the following fields:
Path1: Path of control sample mosdepth file i.e. `./results/03_bamStats`
File1: Name of the control sample mosdepth file i.e. `control.per-base.bed.gz`
Path2: Path of treatment sample mosdepth file i.e. `./results/03_bamStats`
File2: Path of treatment sample mosdepth file i.e.`treatment.per-base.bed.gz`
Cond: Defining the status of the sample. eg. `DMSO treated` and `RAPA treated`

Note: The row order and column names are not important but the column order must be preserved.
 

## Output
1. Adapter Trimmed fastq files.
2. Aligned and indexed bam files. The alignment folder also contains the deduplicated feature file that will be required for making a coverage plot. Users can add/delete rows from this file.
3. Alignment statistics. The mosdepth `per-base-coverage` file that is used by the R script to create coverage plots is produced as well.
  
