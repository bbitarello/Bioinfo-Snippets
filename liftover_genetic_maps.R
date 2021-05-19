liftover_genetic_map<-function(dir="~data/1kg/20130502_phase3_final/", outdir="maps_hg38/", from="hg19", to="hg38", base=".OMNI.interpolated_genetic_map"){
require(bigsnpr)
require(bigstatsr)
require(dplyr)
require(bigreadr)
require(data.table)
assertdir("maps/")
for(chr in 1:22){
basename <- glue::glue("chr{chr}{base}")
info_snp<-fread(glue::glue("{dir}/{basename}")) #path where hg19 map files can be foun
setnames(info_snp, c('SNP', 'pos', 'Dist'))
info_snp[, chr:=chr]
snp_modifyBuild(info_snp, from="hg19", to="hg38", liftOver="~/bin/liftOver")[,.(SNP, pos, Dist)] %>% fwrite2(file=glue::glue("{outdir}{basename}"), col.names=F, sep="\t")
cat(chr, '\n')
}
}
