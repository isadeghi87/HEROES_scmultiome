#!/bin/bash
module load python/3.10.1

bam=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool12/outs/individual_bams/I070_032_atac_sorted.bam 

cnvpytor -root file.pytor -rd $bam
cnvpytor -root file.pytor -his 1000 10000 100000
cnvpytor -root file.pytor -partition 1000 10000 100000
cnvpytor -root file.pytor -call 1000 10000 100000