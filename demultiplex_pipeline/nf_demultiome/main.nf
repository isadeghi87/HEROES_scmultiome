params.outdir = './results'
params.cellrangerInput = '' // Path to Cell Ranger input CSV
params.demultiplexInput = '' // Path to demultiplexing input CSV
params.steps = '' // Empty by default, implies running all steps

boolean runMapping = params.steps.isEmpty() || params.steps.contains('mapping')
boolean runDemultiplexing = params.steps.isEmpty() || params.steps.contains('demultiplexing')


// Check if the cellrangerInput parameter is provided and valid
if (!params.cellrangerInput) {
    error "No Cell Ranger input CSV file provided. Use --cellrangerInput to specify."
}
// Prepare the cellrangerInputs channel
Channel
    .fromPath(params.cellrangerInput)
    .splitCsv(header: true, sep: ',', strip: true)
    .groupBy { it.sample }
    .map { sample, entries ->
        def libCsv = "${workDir}/${sample}_libraries.csv"
        libCsv.withWriter { w ->
            w.writeLine('fastqs,sample,library_type')
            entries.each { entry ->
                w.writeLine("${entry.fastqs},${sample},${entry.library_type}")
            }
        }
        return [sample, file(libCsv)]
    }
    .set { cellrangerInputs }

// Channel for CellRanger output (mapping step)
Channel cellRangerOutCh = Channel.empty()

// Channel for Demultiplexing input
Channel demultiplexInputCh = Channel.empty()

if (runMapping) {
    // Process to run CellRangerARC
   process CellRangerARC {
    tag "${sample}"
    publishDir "${params.outputDir}/CellRanger_ARC", mode: 'copy'

    input:
    tuple val(sample), path(libCsv) from cellrangerInputs

    script:
    """
    module add cellranger-arc/2.0.2
    cellranger-arc count --id="$sample" \\
      --reference="${params.reference}" \\
      --libraries="$libCsv" \\
      --localcores=\$task.cpus
    """
}
}

if (runDemultiplexing) {
    if (params.steps == 'demultiplexing' && params.demultiplexInput) {
        // When only running demultiplexing, expect 5-column CSV
        demultiplexInputCh = Channel
            .fromPath(params.demultiplexInput)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row -> tuple(row.sample, row.pool, file(row.vcfPath), row.bamPath, row.barcodePath) }
    } else if (runMapping) {
        // When running from mapping, construct input channel from CellRanger output
        cellRangerOutCh
            .map { sample, outputDir ->
                // Infer bamPath and barcodePath from CellRanger output directory
                tuple(sample, sample, file("${outputDir}/vcfPath.vcf"), "${outputDir}/bamPath.bam", "${outputDir}/barcodePath.txt")
            }
            .set { demultiplexInputCh }
    }
}

// Filter VCF
process FilterVCF {
    tag "${pool}"

    input:
    tuple val(pool), val(sample), path(vcf) from demultiplexInputCh

    output:
    path "${params.outdir}/filtered/${sample}.filtered.vcf.gz" into filteredVcfCh

    script:
    """
    # Add your command to filter VCF here
    """
}

process MergeVCF {
    output:
    path "${params.outdir}/merged/merged.vcf.gz" into mergedVcfCh

    script:
    """
    # Add your command to merge VCF here
    """
}

process CellSNP {
    tag "${pool}"

    input:
    tuple val(pool), path(vcf) from mergedVcfCh
    tuple val(sample), path(bam), path(barcode) from bamBarcodeCh

    output:
    path "${params.outdir}/cellSNP/${pool}" into cellsnpOutCh

    script:
    """
    # Add your command for CellSNP here
    """
}

process Vireo {
    tag "${pool}"

    input:
    tuple val(pool), path(cellsnpDir) from cellsnpOutCh

    output:
    path "${params.outdir}/Vireo/${pool}" into vireoOutCh

    script:
    """
    # Add your command for Vireo here
    """
}
