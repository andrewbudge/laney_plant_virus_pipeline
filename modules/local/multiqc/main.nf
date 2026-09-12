// Aggregate run QC report from fastp JSON + SortMeRNA logs.
process MULTIQC {
    label 'collect'

    input:
    path report_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc --force -n multiqc_report.html .
    """
}