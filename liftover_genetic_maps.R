library(bigsnpr)
library(bigstatsr)
library(dplyr)
library(bigreadr)
library(data.table)
for(chr in 1:22){
basename <- paste0("chr", chr, ".OMNI.interpolated_genetic_map")
info_snp<-fread(glue::glue("/project/mathilab/data/1kg/20130502_phase3_final/{basename}")) #path where hg19 map files can be foun
setnames(info_snp, c('SNP', 'pos', 'Dist'))
info_snp[, chr:=chr]
snp_modifyBuild(info_snp, from="hg19", to="hg38", liftOver="/project/mathilab/bin/liftOver")[,.(SNP, pos, Dist)] %>% fwrite2(file=glue::glue("maps_hg38/chr{chr}.OMNI.interpolated_genetic_map"), col.names=F, sep="\t")
cat(chr, '\n')
}
