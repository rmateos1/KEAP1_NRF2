#ALL THE GENES BY USING SJ COORDINATES LUAD
library(data.table)
SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip.txt")

SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip = data.frame(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip, stringsAsFactors = F)
#appply scale makes things faster but removes the SJ names. We keep them here.
SJnames = colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip)[substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip),1,2) == "SJ"]
my_palette = colorRampPalette(c("red","white", "green"))(256)
numericSJLUAD = subset(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip,
                       select =  -c(TCGA_ID, isNFE2L2mut, ismutexon2, isKEAP1mut, isCUL3mut , mutlabels , NRE2L2_CNA ,exon2_skipping))
#whichismut= SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip$isNFE2L2mut | SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip$isKEAP1mut 
whichismut= SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip$mutlabels != "None"
#numericSJLUAD =apply(numericSJLUAD, 2, as.numeric)
pvalue = c()
numericSJLUAD = apply(log(numericSJLUAD+ 0.001), 1, scale)
numericSJLUAD = t(numericSJLUAD)
colnames(numericSJLUAD) <- SJnames
#numericSJLUAD =t(scale(t(log(numericSJLUAD+ 0.001))))
for (cont in 1:dim(numericSJLUAD)[2]) {
  pvalue = c(pvalue,
             wilcox.test(numericSJLUAD[, cont][whichismut], numericSJLUAD[, cont][!whichismut], alternative = "greater")$p.value)
}
library(gplots)
library(reshape2)
library(ggplot2)
library(biomaRt)
cont1 = 1

p = 1e-02

pe = (2:8)


#pdf("heatmapsbypvalubothmutated_SJ_LUAD_Normalizedbycol.pdf")
#for (cont1 in 1: 1 )#length(pe))
#{
numericSJLUAD = numericSJLUAD[,order(pvalue, decreasing = F) ]
pvalue = pvalue[order(pvalue, decreasing = F) ]
psignificantisthereanymutinLUAD_forplot =   data.frame(subset(
  SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip,
  select =  c(TCGA_ID, isNFE2L2mut, ismutexon2, isKEAP1mut, isCUL3mut , mutlabels , NRE2L2_CNA ,exon2_skipping)) ,numericSJLUAD[,pvalue < 10 ^ -5])









library(stringr)
#SJnames = colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip)[substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip),1,2) == "SJ"]
#colnames(numericSJLUAD) <- SJnames
#esjeis = SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip[,substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip),1,2) == "SJ"]
#esjeismutlabels = SJandNFE2L2andKEAP1andCUL3mutationinLUAD_CNA_ex2skip$mutlabels != "None"
#SJmutornot = data.frame(esjeismutlabels,esjeis)


significantSJfullfoundinLUAD = colnames(numericSJLUAD[,pvalue < 10 ^ -5])

significantSJfullfoundinLUADcoord = str_split_fixed(significantSJfullfoundinLUAD,"_", 2)[,2]
significantSJfullfoundinLUADcoord = str_split_fixed(significantSJfullfoundinLUADcoord,"\\.", 3)

significantSJfullfoundinLUADcoord = data.frame(significantSJfullfoundinLUADcoord, stringsAsFactors = F)
colnames(significantSJfullfoundinLUADcoord) <- c("chromosome", "start", "end")
significantSJfullfoundinLUADcoord$start = as.numeric(significantSJfullfoundinLUADcoord$start)
significantSJfullfoundinLUADcoord$end = as.numeric(significantSJfullfoundinLUADcoord$end)

genome = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/fullgenomecoordinatesGENCODE", stringsAsFactors = F)
genome = genome[,1:4]
colnames(genome) = c("chrom", "txStart", "txEnd", "ENST_ID")
genome$chrom <- substring(genome$chrom,4)
genome = data.frame(genome)

geneswithSJfullfoundinLUAD= rep("", length(significantSJfullfoundinLUAD))
for(cont in 1:dim(significantSJfullfoundinLUADcoord)[1]){
  caso <- genome %>% 
    filter(chrom == significantSJfullfoundinLUADcoord$chromosome[cont]) %>%
    filter((txStart <= significantSJfullfoundinLUADcoord$start[cont] & txEnd  >= significantSJfullfoundinLUADcoord$start[cont] ) | (txStart <= significantSJfullfoundinLUADcoord$end[cont] & txEnd  >= significantSJfullfoundinLUADcoord$end[cont] )) %>%
    dplyr::select(ENST_ID)
  if(dim(caso)[1] == 0){
    print(cont) 
    print(significantSJfullfoundinLUADcoord[cont,]) 
  }
  geneswithSJfullfoundinLUAD[cont]= paste0(unique(caso$ENST_ID), collapse = ",")
}

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
SJanditsENSTLUAD = data.frame(SJlist = significantSJfullfoundinLUAD, ensembl_transcript_id_version =  geneswithSJfullfoundinLUAD, stringsAsFactors = F)

geneIDlistLUAD = unique(unlist(str_split(geneswithSJfullfoundinLUAD, ",")))

ENSTandGeneIDLUAD = getBM(
  values = geneIDlistLUAD,
  filters = c("ensembl_transcript_id_version"),
  attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version"),
  mart = mart
)


GeneIDforeachENSTofaSJ<- function(listofSJ, reference){
  require(stringr)
  theENST=data.frame(ensembl_transcript_id_version =  unlist(str_split(listofSJ, ",")), stringsAsFactors = F)
  return(paste0(inner_join( theENST , reference, by = "ensembl_transcript_id_version")$ensembl_gene_id_version, collapse = ","))
}

##external_gene_name
eachsampleexternalgenenameLUAD = unlist(lapply(SJanditsENSTLUAD$ensembl_transcript_id_version,GeneIDforeachENSTofaSJ,  reference = ENSTandGeneIDLUAD  ))
SJandENSTandexternamegenenameLUAD = data.frame(SJanditsENSTLUAD, ensembl_gene_id_version =eachsampleexternalgenenameLUAD )

fwrite(SJandENSTandexternamegenenameLUAD, "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandENSTandexternamegenenameLUAD.txt", sep = "\t")
