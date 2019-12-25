#combining SJ with expression

SJselectedforexpression = psignificantisthereanymutinLUSC_forplot
SJselectedforexpression$sampleTCGA_ID = substring(SJselectedforexpression$TCGA_ID, 1,12)
expressionLUSC =  read.table("C:/Users/Raul/Documents/NRF2/Features added version 201912/expressiondatanormalizedofgeneswithsignificantSJLUSC.txt",rownames(T), header = T, stringsAsFactors = F, check.names = F)

expressionLUSC$ensembl_gene_id_version = rownames(expressionLUSC)
IDs = colnames(expressionLUSC)
expressionLUSC = data.frame(t(expressionLUSC), stringsAsFactors = F)
expressionLUSC$sampleTCGA_ID =IDs

#EXPhavingsamplesinSJ = semi_join(expressionLUSC , SJselectedforexpression , by = "sampleTCGA_ID")


EXPhavingsamplesinSJlabeled = inner_join(expressionLUSC, SJselectedforexpression[,c("sampleTCGA_ID", "mutlabels")])

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
#SJandENSTandexternamegenenameLUSC

differenceSJmutornot = t(SJselectedforexpressionlabeledmeanbymutornot[SJselectedforexpressionlabeledmeanbymutornot$mutornot == TRUE,-1]  -  SJselectedforexpressionlabeledmeanbymutornot[SJselectedforexpressionlabeledmeanbymutornot$mutornot == FALSE,-1])
differenceEXPmutornot = t(EXPhavingsamplesinSJlabeledmeanbymutornot[EXPhavingsamplesinSJlabeledmeanbymutornot$mutornot == TRUE,-1]  -  EXPhavingsamplesinSJlabeledmeanbymutornot[EXPhavingsamplesinSJlabeledmeanbymutornot$mutornot == FALSE,-1] )
#for each SJ we add the expression
#some SJ will have more than one gene
#some SJs will share the same gene
SJandENSTandexternamegenenameLUSC$ensembl_gene_id_version<-as.character(SJandENSTandexternamegenenameLUSC$ensembl_gene_id_version)

SjwithitsexpressionLUSC = data.frame()
#this is going to be a loop
caso = rownames(differenceSJmutornot)[1]
for (cont in 1: length(differenceSJmutornot)){
  caso = rownames(differenceSJmutornot)[cont]
  GeneassociatedwiththeSJ  = unique(unlist(str_split(SJandENSTandexternamegenenameLUSC[grep(caso, SJandENSTandexternamegenenameLUSC$SJlist),"ensembl_gene_id_version"], ",")))
  ExpressionassociatedwiththeSJ = differenceEXPmutornot[rownames(differenceEXPmutornot) %in% GeneassociatedwiththeSJ]
  ExpressionassociatedwiththeSJGenesID = rownames(differenceEXPmutornot)[rownames(differenceEXPmutornot) %in% GeneassociatedwiththeSJ]
  if( length(ExpressionassociatedwiththeSJ) >0){
    SjwithitsexpressionLUSC = rbind( SjwithitsexpressionLUSC , cbind(caso , SJdiffvalue = differenceSJmutornot[cont],ensembl_gene_id_version =ExpressionassociatedwiththeSJGenesID,  Expdiffvalue = ExpressionassociatedwiththeSJ))
  }
  
}



#pdf("C:/Users/Raul/Documents/NRF2/Features added version 201912/LUSCSJEXP.pdf")
SjwithitsexpressionLUSC$SJdiffvalue = as.numeric(as.character(SjwithitsexpressionLUSC$SJdiffvalue ))
SjwithitsexpressionLUSC$Expdiffvalue = as.numeric(as.character(SjwithitsexpressionLUSC$Expdiffvalue ))

gLUSC = ggplot(SjwithitsexpressionLUSC , aes(x = SJdiffvalue , y = Expdiffvalue)) +  
  geom_point()+ ggtitle("LUSC:  Exp vs SJ") + theme_minimal()
(gLUSC)
#dev.off()



externalnamefromgeneidLUSC = getBM(
  values = SjwithitsexpressionLUSC$ensembl_gene_id_version,
  filters = c("ensembl_gene_id_version"),
  attributes = c("external_gene_name", "ensembl_gene_id_version"),
  mart = mart
) %>% inner_join(SjwithitsexpressionLUSC)
externalnamefromgeneidLUSC[externalnamefromgeneidLUSC$SJdiffvalue>1 & externalnamefromgeneidLUSC$Expdiffvalue< 4e-6,]

