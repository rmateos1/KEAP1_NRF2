#pdf("C:/Users/Raul/Documents/NRF2/tsne_LUSC_LUSC_bothmutations_exon2NORMALIZED.pdf",height = 12, width =12)

SJselection = fread("C:/Users/Raul/Documents/NRF2//Features added version 201912/SJselectionLUSCallgenes.txt",check.names = F, stringsAsFactors = F)
SJselection = tableSJselection
colnames(SJselection)[1]<- "TCGA_ID"
SJselection$TCGA_ID = apply(str_split_fixed(SJselection$TCGA_ID, "-", 5)[,1:4], 1, paste0 , collapse = "-")
#  SJselection = read.table("C:/Users/Raul/Documents/NRF2/allgenesSJ2/SJselectionLUSCallgenes.txt",check.names = F, stringsAsFactors = F)
#SJselection = data.frame(TCGA_ID = colnames(SJselection) , t(SJselection), stringsAsFactors = F)
labelling = data.frame(numero = c(0,1,2,3,4,5), color = c("#333333", "Red", "Orange" , "#339933" , "Magenta" , "#3300FF"), Gene = c("None","NFE2L2_Exon_2", "NFE2L2_Other_Exon", "KEAP1", "CUL3", "Multihit"), stringsAsFactors = F)



LUSC = read.table("C:/Users/Raul/Documents/NRF2/Pancancer_LUSC.txt", header = T,  stringsAsFactors = F, sep = "\t")
#LUSC = LUSC[LUSC$Variant_Classification != "Silent",]
LUSC = LUSC[LUSC$Variant_Classification != "3UTR",]
LUSC = LUSC[LUSC$Variant_Classification != "5UTR",]
LUSC = LUSC[LUSC$Variant_Classification != "Intron",]
LUSC$TCGA_ID = apply(str_split_fixed(LUSC$Tumor_Sample_Barcode, "-", 5)[,1:4], 1, paste0, collapse = "-")
#LUSC NFE2L2

NFE2L2inLUSC = unique(data.frame(TCGA_ID = LUSC$TCGA_ID , isNFE2L2mut = LUSC$Hugo_Symbol == "NFE2L2", stringsAsFactors = F))
NFE2L2inLUSC = data.frame(TCGA_ID  = names((table(NFE2L2inLUSC)>0)[,2]),isNFE2L2mut =(table(NFE2L2inLUSC)>0)[,2],stringsAsFactors = F)
NFE2L2inLUSC$ismutexon2 <- "None"
#this deals with multiple NFE2L2 mutations in the same sample
mutationsinNFE2L2 = data.frame(LUSC$TCGA_ID[LUSC$TCGA_ID %in% NFE2L2inLUSC$TCGA_ID[NFE2L2inLUSC$isNFE2L2mut] & LUSC$Hugo_Symbol == "NFE2L2"],LUSC$Exon_Number[LUSC$TCGA_ID %in% NFE2L2inLUSC$TCGA_ID[NFE2L2inLUSC$isNFE2L2mut] & LUSC$Hugo_Symbol == "NFE2L2"])
mutationsinNFE2L2 = data.frame(TCGA_ID = rownames(table(mutationsinNFE2L2)),ismutexon2 = data.frame(table(mutationsinNFE2L2)>0)$X2.5)
NFE2L2inLUSC$ismutexon2[match(mutationsinNFE2L2$TCGA_ID , NFE2L2inLUSC$TCGA_ID)] <-  mutationsinNFE2L2$ismutexon2
####
NFE2L2inLUSC$ismutexon2[NFE2L2inLUSC$ismutexon2=="FALSE"]<- "NFE2L2_Other_Exon"
NFE2L2inLUSC$ismutexon2[NFE2L2inLUSC$ismutexon2=="TRUE"]<- "NFE2L2_Exon_2"
NFE2L2inLUSC$ismutexon2 <- factor(NFE2L2inLUSC$ismutexon2, levels = c("None", "NFE2L2_Exon_2", "NFE2L2_Other_Exon" ))
genesandNFE2L2mutationinLUSCNORMALIZED = inner_join(SJselection, NFE2L2inLUSC)

 
#NFE2L2 and KEAP1 IN LUSC

KEAP1inLUSC = unique(data.frame(TCGA_ID = LUSC$TCGA_ID , isKEAP1mut = LUSC$Hugo_Symbol == "KEAP1",stringsAsFactors = F))
KEAP1inLUSC = data.frame(TCGA_ID  = names((table(KEAP1inLUSC)>0)[,2]),isKEAP1mut =(table(KEAP1inLUSC)>0)[,2],stringsAsFactors = F)
SJandNFE2L2andKEAP1mutationinLUSC = inner_join(genesandNFE2L2mutationinLUSCNORMALIZED, KEAP1inLUSC )



#NFE2L2 and KEAP1 and CUL3 IN LUSC

CUL3inLUSC = unique(data.frame(TCGA_ID = LUSC$TCGA_ID , isCUL3mut = LUSC$Hugo_Symbol == "CUL3",stringsAsFactors = F))
CUL3inLUSC = data.frame(TCGA_ID  = names((table(CUL3inLUSC)>0)[,2]),isCUL3mut =(table(CUL3inLUSC)>0)[,2],stringsAsFactors = F)
SJandNFE2L2andKEAP1andCUL3mutationinLUSC = inner_join(SJandNFE2L2andKEAP1mutationinLUSC, CUL3inLUSC )



SJandNFE2L2andKEAP1andCUL3mutationinLUSC$mutlabels <- SJandNFE2L2andKEAP1andCUL3mutationinLUSC$ismutexon2
levels(SJandNFE2L2andKEAP1andCUL3mutationinLUSC$mutlabels) <- c("None"  , "NFE2L2_Exon_2" , "NFE2L2_Other_Exon", "KEAP1", "CUL3", "Multihit")
SJandNFE2L2andKEAP1andCUL3mutationinLUSC$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isKEAP1mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isNFE2L2mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isCUL3mut] <- "KEAP1"
SJandNFE2L2andKEAP1andCUL3mutationinLUSC$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isCUL3mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isNFE2L2mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isKEAP1mut] <- "CUL3"
SJandNFE2L2andKEAP1andCUL3mutationinLUSC$mutlabels[(SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isKEAP1mut + SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isNFE2L2mut + SJandNFE2L2andKEAP1andCUL3mutationinLUSC$isCUL3mut)>1] <- "Multihit"
#tsne = Rtsne(SJandNFE2L2andKEAP1mutationinLUSC[,substring(colnames(SJandNFE2L2andKEAP1mutationinLUSC),1,2) == "SJ"])
#plot(tsne$Y, col = labelling$color[match(SJandNFE2L2andKEAP1mutationinLUSC$mutlabels,labelling$Gene)], pch = "X", main = "NFEL2 and KEAP1 in LUSC")
#legend("bottomright", legend = labelling$Gene,col = labelling$color , pch = "X")

 fwrite(SJandNFE2L2andKEAP1andCUL3mutationinLUSC, "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUSC.txt", sep = "\t")

#dev.off()


