setwd("/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover")

files = list.files(path = "/b06x-isi/b062/g-i/HEROES-AYA/INFORM_HEROES_SNP_Genotyping")

## exclude INF files
files = files[grep('INF',files,invert = T)]

## list VCF files with liftover done

done = list.files(path = "/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files",
                  pattern = "*_hg38.vcf.gz$")

done = gsub('_hg38.vcf.gz','',done)

## obtain those without liftover
new = setdiff(files,done)

## write a data frame for liftover 
ids =  gsub("(^I[0-9]+_[0-9]+).*", "\\1", new)

df = data.frame(new = new, id = ids)

write.table(df,
            file = "/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files/vcf_input.txt",quote = F,
            col.names = F,row.names = F,sep = "\t")
