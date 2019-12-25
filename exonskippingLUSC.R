
folderpath = "C:/Users/Raul/Documents/NRF2/allgenesSJ2withNAD/LUSC/"
outputpath = "C:/Users/Raul/Documents/NRF2/Features added version 201912/SJselectionLUSCallgenes.txt"
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
rbindSJ = function(SJ , folderpath = "C:/Users/Raul/Documents/NRF2/allgenesSJ2withNAD/LUSC/"){
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

NFE2L2SJsetLUSC = fullSJwithlabels[grep(fullSJwithlabels$genes, pattern = "NFE2L2"),]
#exonskippingLUSC = NFE2L2SJsetLUSC[NFE2L2SJsetLUSC$start < 177233249 & NFE2L2SJsetLUSC$start >=  177230308 &  NFE2L2SJsetLUSC$end  >  177233339 & NFE2L2SJsetLUSC$end < 177264727 ,]  
#exonskippingLUSC = NFE2L2SJsetLUSC[NFE2L2SJsetLUSC$start < 177234004	 & NFE2L2SJsetLUSC$start >=  177230308 &  NFE2L2SJsetLUSC$end  >  177234271	 & NFE2L2SJsetLUSC$end < 177264727 ,]  
exonskippingLUSC = NFE2L2SJsetLUSC[( NFE2L2SJsetLUSC$start == 177233339 | NFE2L2SJsetLUSC$start == 177232583 ) &  
                                     (NFE2L2SJsetLUSC$end  ==  177264532	|
                                     NFE2L2SJsetLUSC$end  ==  177263403 |
                                     NFE2L2SJsetLUSC$end  ==  177263437 |
                                     NFE2L2SJsetLUSC$end  ==  177263529	 )
                                   ,]  
#fullSJwithlabels =  fullSJwithlabels[fullSJwithlabels$V11  != "DA",]
exonskippingLUSC = exonskippingLUSC[exonskippingLUSC$score>2,]
exonskippingLUSC$SJ  = str_split_fixed(exonskippingLUSC$SJ, "_", 4)[,3]
exonskippingLUSC$SJ  = apply(data.frame(str_split_fixed(exonskippingLUSC$SJ, "-", 5))[,1:4], 1, FUN = paste0, collapse = "-")
fwrite(exonskippingLUSC, "C:/Users/Raul/Documents/NRF2/Features added version 201912/exonskippingLUSC.txt", sep = "\t")

