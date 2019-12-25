
#SJ = list.files("Splicing_Project/TP53/TCGA-TP53/splicing/", recursive = T)


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
    
    SJ = list.files(folderpath, recursive = T)
    fullSJwithlabels =rbindlist(lapply(SJ, FUN = rbindSJ))
    #to get rid of the data.table format
    fullSJwithlabels = data.frame(fullSJwithlabels, stringsAsFactors = F)
    
    
    
    #fullSJwithlabels =  fullSJwithlabels[fullSJwithlabels$V11  != "DA",]
    
    
    
    #colnames(fullSJwithlabels) <- c("labels" , "chrom" , "start" , "end" , "strand")
    ## POR SI SE QUIERE SE PUEDE DEJAR ESTO PARA GUARDAR LOS SJ SOLOS
    
    if (keepSJfile) {
      write.table(
        fullSJwithlabels,
        paste0(folderpath, "/SJLUAD.txt"),
        sep = "\t",
        row.names = F,
        col.names = T,
        quote = F
      )
    }
    
    ###
    
    #FILTERING OF WHICH ONES MORE THAN 2 READS
    
    fullSJwithlabelsreads = fullSJwithlabels[fullSJwithlabels$score > read_threshold ,]
    fullSJwithlabelsreadscollapsed = apply(
      cbind(
        fullSJwithlabelsreads[, 2],
        ":",
        fullSJwithlabelsreads[, 3],
        "-",
        fullSJwithlabelsreads[, 4]
      ),
      1,
      paste0 ,
      collapse = ""
    )
    selectedSJnames = names(table(fullSJwithlabelsreadscollapsed)[table(fullSJwithlabelsreadscollapsed) > howcommon])
    selectedSJ = table(fullSJwithlabelsreadscollapsed)[table(fullSJwithlabelsreadscollapsed) > howcommon]
    gc()
    collapsedlabel = apply(
      cbind(
        fullSJwithlabels[, 2],
        ":",
        fullSJwithlabels[, 3],
        "-",
        fullSJwithlabels[, 4]
      ),
      1,
      paste0 ,
      collapse = ""
    )
    gc()
    fullsjwithlabelsandcollapsed = cbind(fullSJwithlabels, collapsedlabel)
    fullsjwithlabelsandcollapsed$collapsedlabel <-
      as.character(fullsjwithlabelsandcollapsed$collapsedlabel)
    selectedSJnames = data.frame(collapsedlabel = as.character(selectedSJnames))
    selectionSJamongsamples = inner_join(fullsjwithlabelsandcollapsed, selectedSJnames)
    #selectionSJamongsamples$labels tiene todos los samples en forma de levels, los cuales son los nombres de todas las muestras
    #esto nos permite hacer table con las labels y los reads, y conseguir los ceros, que tambien nos interesan
    
    #little fix here with "labels" to keep using the script as it was
    colnames(selectionSJamongsamples)[1] <- "labels"
    selectionshorten = selectionSJamongsamples[c("labels", "collapsedlabel" , "score")]
    selectionshorten = selectionshorten[!duplicated(selectionshorten[1:2]),]
    tableSJselection = spread(selectionshorten, collapsedlabel , score, fill = 0)
    tableSJselection$labels <-
      str_split_fixed(tableSJselection$labels, "\\.", 4)[, 1]
    tableSJselection$labels <-
      str_split_fixed(tableSJselection$labels, "_", 4)[, 3]
    colnames(tableSJselection) = gsub(":", ".", colnames(tableSJselection))
    colnames(tableSJselection) = gsub("-", ".", colnames(tableSJselection))
    colnames(tableSJselection) = gsub("chr", "", colnames(tableSJselection))
    colnames(tableSJselection)[-1] = paste0("SJ_", colnames(tableSJselection)[-1])
    gc()
    fwrite(
      tableSJselection,
      paste0(outputpath),
      sep = "\t",
      row.names = F,
      col.names = T,
      quote = F
    )
    
    if (heatmap_pdf) {
      pdf(paste0(outputpath, ".tableSJselection.pdf"))
      heatmap(apply(scale(tableSJselection[-1]), 2, as.numeric),
              cexCol =  0.5,
              Colv = NA)
      dev.off()
    }
  #}
