#pdf("C:/Users/Raul/Documents/NRF2/tsne_LUAD_LUAD_bothmutations_exon2NORMALIZED.pdf",height = 12, width =12)

SJselection = fread("C:/Users/Raul/Documents/NRF2//Features added version 201912/SJselectionLUADallgenes.txt",check.names = F, stringsAsFactors = F)
SJselection = tableSJselection
colnames(SJselection)[1]<- "TCGA_ID"
SJselection$TCGA_ID = apply(str_split_fixed(SJselection$TCGA_ID, "-", 5)[,1:4], 1, paste0 , collapse = "-")
#  SJselection = read.table("C:/Users/Raul/Documents/NRF2/allgenesSJ2/SJselectionLUADallgenes.txt",check.names = F, stringsAsFactors = F)
#SJselection = data.frame(TCGA_ID = colnames(SJselection) , t(SJselection), stringsAsFactors = F)
labelling = data.frame(numero = c(0,1,2,3,4,5), color = c("#333333", "Red", "Orange" , "#339933" , "Magenta" , "#3300FF"), Gene = c("None","NFE2L2_Exon_2", "NFE2L2_Other_Exon", "KEAP1", "CUL3", "Multihit"), stringsAsFactors = F)



LUAD = read.table("C:/Users/Raul/Documents/NRF2/Pancancer_LUAD.txt", header = T,  stringsAsFactors = F, sep = "\t")
#LUAD = LUAD[LUAD$Variant_Classification != "Silent",]
LUAD = LUAD[LUAD$Variant_Classification != "3UTR",]
LUAD = LUAD[LUAD$Variant_Classification != "5UTR",]
LUAD = LUAD[LUAD$Variant_Classification != "Intron",]
LUAD$TCGA_ID = apply(str_split_fixed(LUAD$Tumor_Sample_Barcode, "-", 5)[,1:4], 1, paste0, collapse = "-")
#LUAD NFE2L2

NFE2L2inLUAD = unique(data.frame(TCGA_ID = LUAD$TCGA_ID , isNFE2L2mut = LUAD$Hugo_Symbol == "NFE2L2", stringsAsFactors = F))
NFE2L2inLUAD = data.frame(TCGA_ID  = names((table(NFE2L2inLUAD)>0)[,2]),isNFE2L2mut =(table(NFE2L2inLUAD)>0)[,2],stringsAsFactors = F)
NFE2L2inLUAD$ismutexon2 <- "None"
#this deals with multiple NFE2L2 mutations in the same sample
mutationsinNFE2L2 = data.frame(LUAD$TCGA_ID[LUAD$TCGA_ID %in% NFE2L2inLUAD$TCGA_ID[NFE2L2inLUAD$isNFE2L2mut] & LUAD$Hugo_Symbol == "NFE2L2"],LUAD$Exon_Number[LUAD$TCGA_ID %in% NFE2L2inLUAD$TCGA_ID[NFE2L2inLUAD$isNFE2L2mut] & LUAD$Hugo_Symbol == "NFE2L2"])
mutationsinNFE2L2 = data.frame(TCGA_ID = rownames(table(mutationsinNFE2L2)),ismutexon2 = data.frame(table(mutationsinNFE2L2)>0)$X2.5)
NFE2L2inLUAD$ismutexon2[match(mutationsinNFE2L2$TCGA_ID , NFE2L2inLUAD$TCGA_ID)] <-  mutationsinNFE2L2$ismutexon2
####
NFE2L2inLUAD$ismutexon2[NFE2L2inLUAD$ismutexon2=="FALSE"]<- "NFE2L2_Other_Exon"
NFE2L2inLUAD$ismutexon2[NFE2L2inLUAD$ismutexon2=="TRUE"]<- "NFE2L2_Exon_2"
NFE2L2inLUAD$ismutexon2 <- factor(NFE2L2inLUAD$ismutexon2, levels = c("None", "NFE2L2_Exon_2", "NFE2L2_Other_Exon" ))
genesandNFE2L2mutationinLUADNORMALIZED = inner_join(SJselection, NFE2L2inLUAD)


#NFE2L2 and KEAP1 IN LUAD

KEAP1inLUAD = unique(data.frame(TCGA_ID = LUAD$TCGA_ID , isKEAP1mut = LUAD$Hugo_Symbol == "KEAP1",stringsAsFactors = F))
KEAP1inLUAD = data.frame(TCGA_ID  = names((table(KEAP1inLUAD)>0)[,2]),isKEAP1mut =(table(KEAP1inLUAD)>0)[,2],stringsAsFactors = F)
SJandNFE2L2andKEAP1mutationinLUAD = inner_join(genesandNFE2L2mutationinLUADNORMALIZED, KEAP1inLUAD )



#NFE2L2 and KEAP1 and CUL3 IN LUAD

CUL3inLUAD = unique(data.frame(TCGA_ID = LUAD$TCGA_ID , isCUL3mut = LUAD$Hugo_Symbol == "CUL3",stringsAsFactors = F)) 
CUL3inLUAD = data.frame(TCGA_ID  = names((table(CUL3inLUAD)>0)[,2]),isCUL3mut =(table(CUL3inLUAD)>0)[,2],stringsAsFactors = F)
SJandNFE2L2andKEAP1andCUL3mutationinLUAD = inner_join(SJandNFE2L2andKEAP1mutationinLUAD, CUL3inLUAD )



SJandNFE2L2andKEAP1andCUL3mutationinLUAD$mutlabels <- SJandNFE2L2andKEAP1andCUL3mutationinLUAD$ismutexon2
levels(SJandNFE2L2andKEAP1andCUL3mutationinLUAD$mutlabels) <- c("None"  , "NFE2L2_Exon_2" , "NFE2L2_Other_Exon", "KEAP1", "CUL3", "Multihit")
SJandNFE2L2andKEAP1andCUL3mutationinLUAD$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isKEAP1mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isNFE2L2mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isCUL3mut] <- "KEAP1"
SJandNFE2L2andKEAP1andCUL3mutationinLUAD$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isCUL3mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isNFE2L2mut & !SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isKEAP1mut] <- "CUL3"
SJandNFE2L2andKEAP1andCUL3mutationinLUAD$mutlabels[(SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isKEAP1mut + SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isNFE2L2mut + SJandNFE2L2andKEAP1andCUL3mutationinLUAD$isCUL3mut)>1] <- "Multihit"
#tsne = Rtsne(SJandNFE2L2andKEAP1mutationinLUAD[,substring(colnames(SJandNFE2L2andKEAP1mutationinLUAD),1,2) == "SJ"])
#plot(tsne$Y, col = labelling$color[match(SJandNFE2L2andKEAP1mutationinLUAD$mutlabels,labelling$Gene)], pch = "X", main = "NFEL2 and KEAP1 in LUAD")
#legend("bottomright", legend = labelling$Gene,col = labelling$color , pch = "X")

fwrite(SJandNFE2L2andKEAP1andCUL3mutationinLUAD, "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUAD.txt", sep = "\t")

#dev.off()


