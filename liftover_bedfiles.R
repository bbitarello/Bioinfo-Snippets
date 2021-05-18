liftover_BED<-function(x="/home/bbita/R/x86_64-pc-linux-gnu-library/4.0.2/XPASS/extdata/EAS_fourier_ls-all.bed", hg_old="hg19", hg_new="hg38",pop='EAS', outfile=glue::glue("XPASS/inst/extdata/{pop}_fourier_ls-all")){
require(data.table)
require(tidyverse)
options(scipen=999)
temp<-tempfile(tmpdir="tmp-dir")
hg_old<-gsub("H", "h",hg_old)
hg_new<-gsub("h", "H",hg_new)
url <- glue::glue("ftp://hgdownload.cse.ucsc.edu/goldenPath/{hg_old}/liftOver/{hg_old}To{hg_new}.over.chain.gz")
chain <- tempfile(fileext = ".over.chain.gz")
download.file(url, destfile=chain)
x1<-fread(x, na.strings="None")
x1<-na.omit(x1)
x1$start<-as.numeric(x1$start)
x1$stop<-as.numeric(x1$stop)
x1_old<-x1
n1<-nrow(x1_old)
print(paste0("Number of lines is: ", n1))
x1[,start2:=start+1][,id:=paste0(chr, "_", stop)]
x1[,stop2:=stop+1]
x1[,chr2:=as.numeric(gsub("chr","", chr))]
fwrite(x1[order(chr2, start, start2)][, .(chr, start, start2, id)],col.names=F, file=glue::glue("{temp}_start.BED"), sep="\t")
fwrite(x1[order(chr2, stop, stop2)][, .(chr, stop, stop2, id)],col.names=F, file=glue::glue("{temp}_stop.BED"), sep="\t")
system(glue::glue("liftOver {temp}_start.BED {chain} {temp}_start_n.bed {temp}_start_u.bed"))
system(glue::glue("liftOver {temp}_stop.BED {chain} {temp}_stop_n.bed {temp}_stop_u.bed"))
#}
new_start<-fread(glue::glue("{temp}_start_n.bed"))
bad_start<-fread(glue::glue("grep -v '^#' {temp}_start_u.bed"))
new_stop<-fread(glue::glue("{temp}_stop_n.bed"))
bad_stop<-fread(glue::glue("grep -v '^#' {temp}_stop_u.bed"))
all_bad<-unique(sort(c(bad_start$V4, bad_stop$V4, new_start[!(V1 %in% paste0("chr",1:22))]$V4, new_stop[!(V1 %in% paste0("chr",1:22))]$V4)))
print(all_bad)
n2<-nrow(merge(new_start[!(V4 %in% all_bad)], new_stop[!(V4 %in% all_bad)], by=c("V1","V4")))
merge(new_start[!(V4 %in% all_bad)], new_stop[!(V4 %in% all_bad)], by=c("V1", "V4")) %>% mutate(chr=as.numeric(gsub("chr","", V1)), start=as.numeric(V2.x), stop=as.numeric(V2.y), id=V4) %>% arrange(chr, start, stop) %>% select(chr, start, stop) %>% as.data.table %>% fwrite2(file=glue::glue("{outfile}.bed"), sep="\t", col.names=T)
return(print(glue::glue("Wrote lifted {type} file to {outfile}_{hg_new}.bed. {n2} out of {n1} were lifted.")))
}
