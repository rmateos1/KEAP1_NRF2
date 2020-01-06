
folderpath = "C:/Users/Raul/Documents/NRF2/CNV LUAD/samples/"

require("dplyr")
require("stringr")
require("tidyr")
require("data.table")
rbindSJ = function(CNV , folderpath = "C:/Users/Raul/Documents/NRF2/CNV LUAD/samples/"){
  require(data.table)
  return(cbind(CNV, fread(paste0(folderpath, CNV))))
  
}
aliquot = fread("C:/Users/Raul/Documents/NRF2/CNV LUAD/aliquot.tsv")
CNV = list.files(folderpath, recursive = T)
fullCNVwithlabelsLUAD =rbindlist(lapply(CNV, FUN = rbindSJ))
#to get rid of the data.table format
fullCNVwithlabelsLUAD = data.frame(fullCNVwithlabelsLUAD, stringsAsFactors = F)

fullCNVwithlabelswithNRF2LUAD = fullCNVwithlabelsLUAD %>% 
  filter(Chromosome == 2) %>%
  filter(Start <= 177230308) %>% 
  filter(End >= 177264727) %>%
  filter(Segment_Mean > log2(3/2) ) %>%
  rename(aliquot_id =  GDC_Aliquot ) %>%
  left_join(y= aliquot, by="aliquot_id") %>%
  rename(TCGA_ID = sample_submitter_id) %>%
  select(TCGA_ID, Chromosome, Start ,End ,Num_Probes ,Segment_Mean  )


fwrite(fullCNVwithlabelswithNRF2LUAD, "C:/Users/Raul/Documents/NRF2/Features added version 201912/CNAinNFE2L2_LUAD.txt", sep = "\t")
