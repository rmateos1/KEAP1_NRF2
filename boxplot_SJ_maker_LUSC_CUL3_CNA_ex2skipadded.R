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



substring(colnames(psignificantisthereanymutinLUSC_forplot),1,2) == "SJ"
###boxplots###
#50 mas los labels, que son los que no empiezan por SJ
pp = psignificantisthereanymutinLUSC_forplot[,1:(25+sum(substring(colnames(psignificantisthereanymutinLUSC_forplot),1,2) != "SJ"))]

q = melt(
  pp,
  id.vars =colnames(pp)[substring(colnames(pp),1,2) != "SJ"]
)
#q$value <- log(q$value + 0.001)

q$mutlabels = factor(q$mutlabels, levels = c("None" , "NFE2L2_Other_Exon" ,"NFE2L2_Exon_2", "KEAP1" , "CUL3"  ,"CNA" ,"Exon2_Skipping", "Multihit"))



myColors <- brewer.pal(8,"Set1")
names(myColors) <- levels(q$mutlabels)
ggbx = ggplot(data = q,  aes(  x = variable,    y = value,    fill = mutlabels ) )+
  geom_boxplot(position = position_dodge(1),outlier.size=0.5)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1) ) + ggtitle("LUSC normalized greater test") +
 scale_fill_manual(values = myColors)
#geom_dotplot(binaxis='y', stackdir='center',
#               position=position_dodge(1), binwidth = 0.1)


ggsave(
  plot = ggbx,
  file = "C:/Users/Raul/Documents/NRF2/Features added version 201912/ggboxlogpvaluesiskeap1mutinLUSC_forplottop25ormalizedandgreaterNEWFEATURESwithNAD.pdf",
  height = 10,
  width = 50 ,
  units = "cm",
  device =
    "pdf"
)


