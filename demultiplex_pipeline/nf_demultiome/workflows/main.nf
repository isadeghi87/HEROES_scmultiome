/*
 * Workflow logic
 */
workflow {
    cellrangerInputCh = Channel.empty()
    demultiplexInputCh = Channel.empty()
    bamBarcodeCh = Channel.empty()

    if (runMapping) {
        cellrangerInputCh = Channel
            .fromPath(params.cellrangerInput)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(row.sample, row.fastqs, row.library_type) }
            .flatMap { sample, fastqs, library_type ->
                // Generate libraries CSV for each sample
                def libCsv = "${workDir}/${sample}_libraries.csv"
                new File(libCsv).withWriter { w ->
                    w.writeLine("fastqs,sample,library_type")
                    w.writeLine("$fastqs,$sample,$library_type")
                }
                return tuple(sample, file(libCsv))
            }

        CellRangerARC(cellrangerInputCh)
        bamBarcodeCh = cellrangerOutCh.map { sample, outDir ->
            def bamPath = "${outDir}/${sample}/outs/possorted_bam.bam"
            def barcodePath = "${outDir}/${sample}/outs/barcoded_fastqs.tsv"
            return tuple(sample, file(bamPath), file(barcodePath))
        }
    }

    if (runDemultiplexing) {
        demultiplexInputCh = Channel
            .fromPath(params.demultiplexInput)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row -> tuple(row.pool, row.sample, row.vcfPath, row.bamPath, row.barcodePath) }

        FilterVCF(demultiplexInputCh)
        MergeVCF(filteredVcfCh)
        CellSNP(mergedVcfCh, bamBarcodeCh)
        Vireo(cellsnpOutCh)
    }
}
