SJandNFE2L2andKEAP1andCUL3mutationinLUSC = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUSC.txt")
CNAinNFE2L2_LUSC = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/CNAinNFE2L2_LUSC.txt")
exskipinNFE2L2_LUSC = fread("C:/Users/Raul/Documents/NRF2/Features added version 201912/exonskippingLUSC.txt")
exskipinNFE2L2_LUSC_shorten = unique(data.frame(TCGA_ID = exskipinNFE2L2_LUSC$SJ, exon2_skipping = TRUE))
CNAinNFE2L2_LUSC_shorten = unique(data.frame(TCGA_ID =CNAinNFE2L2_LUSC$TCGA_ID, NRE2L2_CNA = TRUE ))
TCGA_ID_CNA_exskip = SJandNFE2L2andKEAP1andCUL3mutationinLUSC[,1] %>%
   left_join(CNAinNFE2L2_LUSC_shorten, by = "TCGA_ID") %>%
  left_join(exskipinNFE2L2_LUSC_shorten , by = "TCGA_ID")
TCGA_ID_CNA_exskip$NRE2L2_CNA[is.na(TCGA_ID_CNA_exskip$NRE2L2_CNA)] <- FALSE
TCGA_ID_CNA_exskip$exon2_skipping[is.na(TCGA_ID_CNA_exskip$exon2_skipping)] <- FALSE
SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip = data.frame(SJandNFE2L2andKEAP1andCUL3mutationinLUSC,
                                                              NRE2L2_CNA = TCGA_ID_CNA_exskip$NRE2L2_CNA ,
                                                              exon2_skipping  = TCGA_ID_CNA_exskip$exon2_skipping)
SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$NRE2L2_CNA] <- "CNA"
SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$mutlabels[SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$exon2_skipping] <- "Exon2_Skipping"
SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$mutlabels[rowSums(data.frame(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$isNFE2L2mut,
    SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$isKEAP1mut,
    SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$isCUL3mut,
    SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$NRE2L2_CNA,
    SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip$exon2_skipping))> 1] <- "Multihit"
fwrite(SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip, "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJandNFE2L2andKEAP1andCUL3mutationinLUSC_CNA_ex2skip.txt", sep = "\t")
