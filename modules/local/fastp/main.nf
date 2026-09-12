// Adapter/quality trimming with fastp (paired-end).
process FASTP {
    tag "${meta.id}"
    label 'qc'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_1.fastp.fastq.gz"),
                     path("${meta.id}_2.fastp.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_1.fastp.fastq.gz -O ${meta.id}_2.fastp.fastq.gz \\
        ${task.ext.args} \\
        -h ${meta.id}.fastp.html \\
        -j ${meta.id}.fastp.json \\
        -R "${meta.id} fastp report" \\
        -w ${task.cpus}
    """
}