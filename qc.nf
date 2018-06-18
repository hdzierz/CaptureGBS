#!/usr/bin/env nextflow

// Paired end reads. Configure this in nextflow.conf or via command line option --reads. The files specified below are test files
params.reads = "$baseDir/KiwiTestData/*.R{1,2}.fq.gz"
params.genome = "$baseDir/KiwiTestData/kiwitest.fasta.gz"
params.targets_file = "$baseDir/targets.txt"

params.trimmomaticv = "0.36"
params.trimmomatico = "SLIDINGWINDOW:5:20 MINLEN:70"
params.pearv= "0.9.10"
params.pearo = "-v 5 -t 70"
params.fastqcv = "0.11.2"
params.multiqcv = "1.2"

/*
 * The reference genome file
 */
genome_file = file(params.genome)
targets_file = file(params.targets_file)

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }
  


process trimm{
    module "Trimmomatic/$params.trimmomaticv"
    cpus 1

    input:
    file targets_file from targets_file
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, "*.fastq.gz" into trimmed

    publishDir "001.trimm", mode: 'copy', overwrite: true


    """
    java -jar -Xms8096m -Xmx16192m \$TRIMMOMATIC PE -threads 1 $reads \
    ${pair_id}.trim.R1.fastq.gz ${pair_id}.trimUn.R1.fastq.gz ${pair_id}.trim.R2.fastq.gz  ${pair_id}.trimUn.R2.fastq.gz  \
    ${params.trimmomatico}

    """

}

process merge_reads{
    module "pear/$params.pearv"

    publishDir "002.merge_reads", mode: 'copy', overwrite: true

    input:
    set pair_id, file(reads) from trimmed

    output:
    set pair_id, "${pair_id}.merge.*.gz" into merged

    """
    pear -f ${pair_id}.trim.R1.fastq.gz -r ${pair_id}.trim.R2.fastq.gz -o ${pair_id}.merge.fastq ${params.pearo}
    gzip *.fastq
    """
}


process qc{
    module "FastQC/$params.fastqcv"

    input:
    set pair_id, file(reads) from merged

    output:
    set pair_id, "${pair_id}*" into qc 

    publishDir "002.merge_reads", mode: 'copy', overwrite: true

    """
    fastqc  *assembled.fastq.gz
    """

}


process multiqc{
    module "MultiQC/${params.multiqcv}"

    input:
    set pair_id, file(qc) from qc

    output:
    set pair_id, "*.html" into mqc

    publishDir "002.merge_reads", mode: 'copy', overwrite: true

    """
    multiqc .
    """
}

