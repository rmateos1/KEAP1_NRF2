
folderpath = "C:/Users/Raul/Documents/NRF2/allgenesSJ2withNAD/LUAD/"
outputpath = "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJselectionLUADallgenes.txt"
keepSJfile = FALSE 
read_threshold = 2
howcommon = 10
heatmap_pdf = FALSE

#SJfilter <-
# function(folderpath,
#           outputpath,
#           read_threshold = 2,
#           howcommon = 102,
#           keepSJfile = FALSE,
#           heatmap_pdf = FALSE) {
require("dplyr")
require("stringr")
require("tidyr")
require("data.table")
rbindSJ = function(SJ , folderpath = "C:/Users/Raul/Documents/NRF2/allgenesSJ2withNAD/LUAD/"){
  require(data.table)
  return(cbind(SJ, fread(paste0(folderpath, SJ))))
  
}

#2	177233339	177264532
#2	177233339	177263403
#2  177233339 177263437
#2  177233339 177263529
SJ = list.files(folderpath, recursive = T)
fullSJwithlabels =rbindlist(lapply(SJ, FUN = rbindSJ))
#to get rid of the data.table format
fullSJwithlabels = data.frame(fullSJwithlabels, stringsAsFactors = F)

NFE2L2SJsetLUAD = fullSJwithlabels[grep(fullSJwithlabels$genes, pattern = "NFE2L2"),]
#exonskippingLUAD = NFE2L2SJsetLUAD[NFE2L2SJsetLUAD$start < 177233249 & NFE2L2SJsetLUAD$start >=  177230308 &  NFE2L2SJsetLUAD$end  >  177233339 & NFE2L2SJsetLUAD$end < 177264727 ,]  
#exonskippingLUAD = NFE2L2SJsetLUAD[NFE2L2SJsetLUAD$start < 177234004	 & NFE2L2SJsetLUAD$start >=  177230308 &  NFE2L2SJsetLUAD$end  >  177234271	 & NFE2L2SJsetLUAD$end < 177264727 ,]  
exonskippingLUAD = NFE2L2SJsetLUAD[( NFE2L2SJsetLUAD$start == 177233339 | NFE2L2SJsetLUAD$start == 177232583 ) &  
                                     (NFE2L2SJsetLUAD$end  ==  177264532	|
                                        NFE2L2SJsetLUAD$end  ==  177263403 |
                                        NFE2L2SJsetLUAD$end  ==  177263437 |
                                        NFE2L2SJsetLUAD$end  ==  177263529	 )
                                   ,]  
#fullSJwithlabels =  fullSJwithlabels[fullSJwithlabels$V11  != "DA",]
exonskippingLUAD = exonskippingLUAD[exonskippingLUAD$score>2,]
exonskippingLUAD$SJ  = str_split_fixed(exonskippingLUAD$SJ, "_", 4)[,3]
exonskippingLUAD$SJ  = apply(data.frame(str_split_fixed(exonskippingLUAD$SJ, "-", 5))[,1:4], 1, FUN = paste0, collapse = "-")
fwrite(exonskippingLUAD, "C:/Users/Raul/Documents/NRF2/Features added version 201912/exonskippingLUAD.txt", sep = "\t")

