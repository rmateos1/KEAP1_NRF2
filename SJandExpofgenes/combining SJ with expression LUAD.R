#combining SJ with expression

SJselectedforexpression = psignificantisthereanymutinLUAD_forplot
SJselectedforexpression$sampleTCGA_ID = substring(SJselectedforexpression$TCGA_ID, 1,12)
expressionLUAD =  read.table("C:/Users/Raul/Documents/NRF2/Features added version 201912/expressiondatanormalizedofgeneswithsignificantSJLUAD.txt",rownames(T), header = T, stringsAsFactors = F, check.names = F)

expressionLUAD$ensembl_gene_id_version = rownames(expressionLUAD)
IDs = colnames(expressionLUAD)
expressionLUAD = data.frame(t(expressionLUAD), stringsAsFactors = F)
expressionLUAD$sampleTCGA_ID =IDs

#EXPhavingsamplesinSJ = semi_join(expressionLUAD , SJselectedforexpression , by = "sampleTCGA_ID")


EXPhavingsamplesinSJlabeled = inner_join(expressionLUAD, SJselectedforexpression[,c("sampleTCGA_ID", "mutlabels")])

EXPhavingsamplesinSJlabeled$mutornot = EXPhavingsamplesinSJlabeled$mutlabels != "None"
EXPhavingsamplesinSJlabeledmeanbymutornot = EXPhavingsamplesinSJlabeled %>% 
  group_by(mutornot) %>%
  dplyr::select(starts_with("ENS")) %>%
  mutate_all(as.numeric)  %>%
  dplyr::summarize_all(funs(median))

SJselectedforexpression$mutornot = SJselectedforexpression$mutlabels != "None"
SJselectedforexpressionlabeledmeanbymutornot = SJselectedforexpression %>% 
  group_by(mutornot) %>%
  dplyr::select(starts_with("SJ")) %>%
  mutate_all(as.numeric)  %>%
  dplyr::summarize_all(funs(median))

#SJselectedforexpressionlabeledmeanbymutornot
#SJandENSTandexternamegenenameLUAD

differenceSJmutornot = t(SJselectedforexpressionlabeledmeanbymutornot[SJselectedforexpressionlabeledmeanbymutornot$mutornot == TRUE,-1]  -  SJselectedforexpressionlabeledmeanbymutornot[SJselectedforexpressionlabeledmeanbymutornot$mutornot == FALSE,-1])
differenceEXPmutornot = t(EXPhavingsamplesinSJlabeledmeanbymutornot[EXPhavingsamplesinSJlabeledmeanbymutornot$mutornot == TRUE,-1]  -  EXPhavingsamplesinSJlabeledmeanbymutornot[EXPhavingsamplesinSJlabeledmeanbymutornot$mutornot == FALSE,-1] )
#for each SJ we add the expression
#some SJ will have more than one gene
#some SJs will share the same gene
SJandENSTandexternamegenenameLUAD$ensembl_gene_id_version<-as.character(SJandENSTandexternamegenenameLUAD$ensembl_gene_id_version)

SjwithitsexpressionLUAD = data.frame()
#this is going to be a loop
caso = rownames(differenceSJmutornot)[1]
for (cont in 1: length(differenceSJmutornot)){
  caso = rownames(differenceSJmutornot)[cont]
  GeneassociatedwiththeSJ  = unique(unlist(str_split(SJandENSTandexternamegenenameLUAD[grep(caso, SJandENSTandexternamegenenameLUAD$SJlist),"ensembl_gene_id_version"], ",")))
  ExpressionassociatedwiththeSJ = differenceEXPmutornot[rownames(differenceEXPmutornot) %in% GeneassociatedwiththeSJ]
  ExpressionassociatedwiththeSJGenesID = rownames(differenceEXPmutornot)[rownames(differenceEXPmutornot) %in% GeneassociatedwiththeSJ]
  if( length(ExpressionassociatedwiththeSJ) >0){
    SjwithitsexpressionLUAD = rbind( SjwithitsexpressionLUAD , cbind(caso , SJdiffvalue = differenceSJmutornot[cont],ensembl_gene_id_version =ExpressionassociatedwiththeSJGenesID,  Expdiffvalue = ExpressionassociatedwiththeSJ))
  }
  
}





pdf("C:/Users/Raul/Documents/NRF2/Features added version 201912/LUADSJEXP.pdf")
SjwithitsexpressionLUAD$SJdiffvalue = as.numeric(as.character(SjwithitsexpressionLUAD$SJdiffvalue ))
SjwithitsexpressionLUAD$Expdiffvalue = as.numeric(as.character(SjwithitsexpressionLUAD$Expdiffvalue ))

gLUAD = ggplot(SjwithitsexpressionLUAD , aes(x = SJdiffvalue , y = Expdiffvalue)) +  
  geom_point()+ ggtitle("LUAD:  Exp vs SJ") + theme_minimal()
(gLUAD)
dev.off()



externalnamefromgeneidLUAD = getBM(
  values = SjwithitsexpressionLUAD$ensembl_gene_id_version,
  filters = c("ensembl_gene_id_version"),
  attributes = c("external_gene_name", "ensembl_gene_id_version"),
  mart = mart
) %>% inner_join(SjwithitsexpressionLUAD)
externalnamefromgeneidLUAD[externalnamefromgeneidLUAD$SJdiffvalue>1 & externalnamefromgeneidLUAD$Expdiffvalue< 2e-6,]
