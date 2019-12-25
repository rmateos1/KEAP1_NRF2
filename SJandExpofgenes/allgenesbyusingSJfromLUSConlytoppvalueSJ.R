#ALL THE GENES BY USING SJ COORDINATES LUSC
library(data.table)
SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip.txt")

SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip = data.frame(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip, stringsAsFactors = F)
#appply scale makes things faster but removes the SJ names. We keep them here.
SJnames = colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip)[substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip),1,2) == "SJ"]
my_palette = colorRampPalette(c("red","white", "green"))(256)
numericSJLUSC = subset(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip,
                       select =  -c(TCGA_ID, isNFE2L2mut, ismutexon2, isKEAP1mut, isCUL3mut , mutlabels , NRE2L2_CNA ,exon2_skipping))
#whichismut= SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$isNFE2L2mut | SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$isKEAP1mut 
whichismut= SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$mutlabels != "None"
#numericSJLUSC =apply(numericSJLUSC, 2, as.numeric)
pvalue = c()
numericSJLUSC = apply(log(numericSJLUSC+ 0.001), 1, scale)
numericSJLUSC = t(numericSJLUSC)
colnames(numericSJLUSC) <- SJnames
#numericSJLUSC =t(scale(t(log(numericSJLUSC+ 0.001))))
for (cont in 1:dim(numericSJLUSC)[2]) {
  pvalue = c(pvalue,
             wilcox.test(numericSJLUSC[, cont][whichismut], numericSJLUSC[, cont][!whichismut], alternative = "greater")$p.value)
}
library(gplots)
library(reshape2)
library(ggplot2)
library(biomaRt)
cont1 = 1

p = 1e-02

pe = (2:8)


#pdf("heatmapsbypvalubothmutated_SJ_LUSC_Normalizedbycol.pdf")
#for (cont1 in 1: 1 )#length(pe))
#{
numericSJLUSC = numericSJLUSC[,order(pvalue, decreasing = F) ]
pvalue = pvalue[order(pvalue, decreasing = F) ]
psignificantisthereanymutinLUSC_forplot =   data.frame(subset(
  SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip,
  select =  c(TCGA_ID, isNFE2L2mut, ismutexon2, isKEAP1mut, isCUL3mut , mutlabels , NRE2L2_CNA ,exon2_skipping)) ,numericSJLUSC[,pvalue < 10 ^ -5])









library(stringr)
#SJnames = colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip)[substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip),1,2) == "SJ"]
#colnames(numericSJLUSC) <- SJnames
#esjeis = SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip[,substring(colnames(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip),1,2) == "SJ"]
#esjeismutlabels = SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$mutlabels != "None"
#SJmutornot = data.frame(esjeismutlabels,esjeis)


significantSJfullfoundinLUSC = colnames(numericSJLUSC[,pvalue < 10 ^ -5])

significantSJfullfoundinLUSCcoord = str_split_fixed(significantSJfullfoundinLUSC,"_", 2)[,2]
significantSJfullfoundinLUSCcoord = str_split_fixed(significantSJfullfoundinLUSCcoord,"\\.", 3)

significantSJfullfoundinLUSCcoord = data.frame(significantSJfullfoundinLUSCcoord, stringsAsFactors = F)
colnames(significantSJfullfoundinLUSCcoord) <- c("chromosome", "start", "end")
significantSJfullfoundinLUSCcoord$start = as.numeric(significantSJfullfoundinLUSCcoord$start)
significantSJfullfoundinLUSCcoord$end = as.numeric(significantSJfullfoundinLUSCcoord$end)

genome = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/fullgenomecoordinatesGENCODE", stringsAsFactors = F)
genome = genome[,1:4]
colnames(genome) = c("chrom", "txStart", "txEnd", "ENST_ID")
genome$chrom <- substring(genome$chrom,4)
genome = data.frame(genome)

geneswithSJfullfoundinLUSC= rep("", length(significantSJfullfoundinLUSC))
for(cont in 1:dim(significantSJfullfoundinLUSCcoord)[1]){
  caso <- genome %>% 
    filter(chrom == significantSJfullfoundinLUSCcoord$chromosome[cont]) %>%
    filter((txStart <= significantSJfullfoundinLUSCcoord$start[cont] & txEnd  >= significantSJfullfoundinLUSCcoord$start[cont] ) | (txStart <= significantSJfullfoundinLUSCcoord$end[cont] & txEnd  >= significantSJfullfoundinLUSCcoord$end[cont] )) %>%
    dplyr::select(ENST_ID)
  if(dim(caso)[1] == 0){
    print(cont) 
    print(significantSJfullfoundinLUSCcoord[cont,]) 
  }
  geneswithSJfullfoundinLUSC[cont]= paste0(unique(caso$ENST_ID), collapse = ",")
}

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
SJanditsENSTLUSC = data.frame(SJlist = significantSJfullfoundinLUSC, ensembl_transcript_id_version =  geneswithSJfullfoundinLUSC, stringsAsFactors = F)

geneIDlistLUSC = unique(unlist(str_split(geneswithSJfullfoundinLUSC, ",")))

ENSTandGeneIDLUSC = getBM(
  values = geneIDlistLUSC,
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

eachsampleexternalgenenameLUSC = unlist(lapply(SJanditsENSTLUSC$ensembl_transcript_id_version,GeneIDforeachENSTofaSJ,  reference = ENSTandGeneIDLUSC  ))
SJandENSTandexternamegenenameLUSC = data.frame(SJanditsENSTLUSC, ensembl_gene_id_version =eachsampleexternalgenenameLUSC )

fwrite(SJandENSTandexternamegenenameLUSC, "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandENSTandexternamegenenameLUSC.txt", sep = "\t")
