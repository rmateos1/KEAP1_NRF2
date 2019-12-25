##write.table(SJandNFE2L2andKEAP1mutationinLUSC , "SJandNFE2L2andKEAP1mutationinLUSC.txt", sep = "\t", col.names = T, row.names = F, quote = F )
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




hLUSC = psignificantisthereanymutinLUSC_forplot[,substring(colnames(psignificantisthereanymutinLUSC_forplot),1,2) == "SJ"]
mutationsandIDLUSC = psignificantisthereanymutinLUSC_forplot[,substring(colnames(psignificantisthereanymutinLUSC_forplot),1,2) != "SJ"]
hLUSC = as.matrix(hLUSC)
rownames(hLUSC) = mutationsandIDLUSC$mutlabels
q = rownames(hLUSC)
colores = c()
for (cont in 1:length(q)) {
  if (q[cont] == "None") {
    colores = c(colores, "White")
  } else if (q[cont] == "KEAP1") {
    colores = c(colores, "Orange")
  } else if (q[cont] == "NFE2L2_Other_Exon") {
    colores = c(colores, "Gray")
  } else if (q[cont] == "NFE2L2_Exon_2") {
    colores = c(colores, "Blue")
  } else if (q[cont] == "CUL3") {
    colores = c(colores, "Brown")
  } else if (q[cont] == "CNA") {
    colores = c(colores, "Purple")
  } else if (q[cont] == "Exon2_Skipping") {
    colores = c(colores, "Green")
  } else{
    colores = c(colores, "Red")
  }
}

#loghh =t(scale(t(log(h+ 0.001))))
loghh = hLUSC


loghhh = loghh
loghhh[loghhh > 3] <- 3
loghhh[loghhh < -3] <- -3
my_palette = colorRampPalette(c("white","white", "green"))(256)
pdf("C:/Users/Raul/Documents/NRF2/Features added version 201912/heatmapsbypvalubothmutated_SJ_LUSC_logandscaledbyrowcolorthresholds3minus3normalizedpretestgreater_NEWFEATURESwithNAD.pdf")
heatmap.2(
  loghhh,
  na.rm = T,
  cexCol = 0.6,
  cexRow = 0.01,
  margins = c(8, 8),
  trace = "none",
  col = my_palette, 
  RowSideColors = colores,
  #dendogram = NULL,
  #Colv = NULL,
  main = paste0("LUSC: pvalue = 0." , paste0(rep(0, 4), collapse = ""), 1),
)

dev.off()

