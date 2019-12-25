SJIdatapath =  "C:/Users/Raul/Documents/TP53/anygenepipeline/SJIdata.txt"
SJImetadatapath = "C:/Users/Raul/Documents/TP53/anygenepipeline/SJImetadata.txt"
outputpath = "C:/Users/Raul/Documents/TP53/anygenepipeline/"
nrounds = 2
psignificantisthereanymutinLUSC_forplot

dataforNBD = psignificantisthereanymutinLUSC_forplot

dataforNBD = dataforNBD[dataforNBD$mutlabels != "NFE2L2_Other_Exon", ]
dataforNBD$mutornot =  dataforNBD$mutlabels != "None"
dataforNBD$mutornot =  dataforNBD$mutlabels != "None"
dataforNBD$mutornot[dataforNBD$mutlabels != "None"] <- "Mut"
dataforNBD$mutornot[dataforNBD$mutlabels == "None"] <- "Not_Mut"
dataforNBD$mutornot= factor(dataforNBD$mutornot)

dataforNBD[,substring(colnames(dataforNBD),1,3) == "SJ_"]  = dataforNBD[,substring(colnames(dataforNBD),1,3) == "SJ_"] -  min(dataforNBD[,substring(colnames(dataforNBD),1,3) == "SJ_"])

           SJImetadatapath
           outputpath
           nrounds = 10 
           ncores = 6
           newWTsize = FALSE
           WTsize = c()


    require(doParallel)
    require(foreach)
    registerDoParallel(cores = ncores)
    clasificaciones = c()
    tableclassifications = c()
    my_palette = colorRampPalette(c("white", "blue"))(256)
    dataforNBD$mutornot <- relevel(dataforNBD$mutornot, ref = "Not_Mut")
    dataforNBD$mutornot = droplevels(dataforNBD$mutornot)
    if (!newWTsize)  {
      WTsize =  sum(dataforNBD$mutornot == "Not_Mut")
    }
   # pdf(paste0(outputpath, "feature_relevance_heatmap.pdf"))
    #for (rounds in 1:nrounds) {
      print(rounds)
      tooptimize = function(r , k = dmfeature) {
        sum(digamma(k + r)) - length(k) * digamma(r) + length(k) * log(r / (r + sum(k /
                                                                                      length(k))))
      }
      
      ##################################################
      ######thetrainingset and test set generation######
      ##################################################
      
      
      
      part1 = dataforNBD[sample(which(dataforNBD$mutornot == "Not_Mut"), length(which(dataforNBD$mutornot == "Not_Mut")))[1:WTsize],]
      whichsampletestWT = sample(1:dim(part1)[1], (dim(part1)[1]) / 2)
      sampletestWT = part1[whichsampletestWT,]
      sampletrainWT = part1[-whichsampletestWT,]
      part2 =  dataforNBD[dataforNBD$mutornot != "Not_Mut",]
      
      sampletestmut = c()
      sampletrainmut = c()
      motifs = as.character(unique(as.character(part2$mutornot)))
      for (mot in 1:length(motifs)) {
        casosmotivo = part2[part2$mutornot == motifs[mot],]
        whichcasosmotivo =  sample(1:dim(casosmotivo)[1], (dim(casosmotivo)[1]) /
                                     2)
        sampletrainmut = rbind(sampletrainmut, casosmotivo[whichcasosmotivo,])
        sampletestmut = rbind(sampletestmut, casosmotivo[-whichcasosmotivo,])
        
      }
      
      
      SJItrainset = rbind(sampletrainWT, sampletrainmut)
      SJItrainset$mutornot <- droplevels(SJItrainset$mutornot)
      testset = rbind(sampletestWT, sampletestmut)
  #    testset$mutornot <- droplevels(testset$mutornot)
      
      
      
      
      Motif = SJItrainset$mutornot
      
      
      SJItrainset = SJItrainset[, substring(colnames(SJItrainset), 1, 2) == "SJ" |
                                  substring(colnames(SJItrainset), 1, 2) == "I_"]
      #trainTP53bpcount = SJItrainset$TP53bpcount
      WTmedians = c()
      motif2medians = c()
      frequencies = data.frame()
      
      uniqueselectedMotifs = names(table(Motif)[table(Motif) > 4])
      SJItrainsetselected = SJItrainset[Motif %in% uniqueselectedMotifs ,]
      selectedMotifs = Motif[Motif %in% uniqueselectedMotifs]
      ropt = matrix(
        0,
        nrow = length(uniqueselectedMotifs),
        ncol = dim(SJItrainsetselected)[2]
      )
      popt = ropt
      
      py = table(selectedMotifs) / sum(table(selectedMotifs))
      
      for (cont in 1:length(uniqueselectedMotifs)) {
        themotif = uniqueselectedMotifs[cont]
        dm = SJItrainsetselected[selectedMotifs == uniqueselectedMotifs[cont],]
        dm[dm == 0] <- exp(-10)
        for (feature in 1:dim(dm)[2]) {
          dmfeature  = dm[, feature]
          opt =  optim(
            10,
            tooptimize,
            method = "BFGS",
            control = list(maxit = 100000, abstol = 0)
          )
          ropt[cont, feature] =  opt$par
          popt[cont, feature] = sum(dmfeature) / (length(dmfeature) * ropt[cont, feature] + sum(dmfeature))
        }
        
      }
      
      
      lamoda = ropt * (popt / (1 - popt))
      colnames(lamoda) <- colnames(SJItrainset)
      rownames(lamoda) <- uniqueselectedMotifs
      heatmap.2(
        lamoda,
        na.rm = T,
        Rowv = NULL,
        Colv = NULL,
        cexCol = 0.5,
        cexRow = 0.5,
        margins = c(8, 8),
        trace = "none",
        col = my_palette
      )
      
      
      
      
      SJItestset = testset[, substring(colnames(testset), 1, 2) == "SJ" |
                             substring(colnames(testset), 1, 2) == "I_"]
      SJItestset[SJItestset == 0] <- exp(-10)
      Motif = testset$mutornot
      sides = str_split_fixed(colnames(SJItestset), "\\.", 3)
      lengthies = as.numeric(sides[, 3]) - as.numeric(sides[, 2])
      lengthies[substring(colnames(SJItestset), 1, 2) == "SJ"] <-  1
      uniqueselectedMotifs = names(table(Motif)[table(Motif) > 4])
      SJItestsetselected = SJItestset[Motif %in% uniqueselectedMotifs ,]
      selectedMotifs = Motif[Motif %in% uniqueselectedMotifs]
      
      coleccionscorestesting = c()
      coleccionscorestestingexpminusmin = c(rep(0, 10))
      
      
      motifprobability = 0
      namesselectedMotifs  = names(table(selectedMotifs))
      clasificacion = c()
      SJItestset[SJItestset == 0] <- exp(-10)
      ropt = mpfr(ropt, 100)
      
      logfactoriald = log(factorial(SJItestset))
      loggammaropt = log(gamma(ropt))
      logpopt = log(popt)
      ropttimeslogoneminuspopt = ropt * log(1 - popt)
      ropttimeslogoneminuspoptMINUSloggammaropt = ropttimeslogoneminuspopt - loggammaropt
      
      selection = c()
      selection <- foreach(caso = 1:dim(SJItestset)[1]) %dopar%
        {
          library(Rmpfr)
          start_time <- Sys.time()
          scores = c()
          scores = mpfr(scores, 100)
          SJItestsetloop = mpfr(as.numeric(SJItestset[caso,]) , 100)
          logfactorialdd = logfactoriald[caso,]
          for (motivo in 1:(length(uniqueselectedMotifs))) {
            motifscore = 0
            for (feature in 1:length(SJItestsetloop)) {
              motifscore = motifscore + log(gamma(ropt[motivo, feature] + SJItestsetloop[feature])) -
                logfactorialdd[feature] +
                SJItestsetloop[feature] * logpopt[motivo, feature]  +
                ropttimeslogoneminuspoptMINUSloggammaropt[motivo, feature]
            }
            scores[motivo] = log(py[motivo]) +   motifscore
            
          }
          print(scores)
          
          which.max(as.numeric(exp(scores) / sum(exp(scores))))
          
          
        }
      selection = unlist(selection)
      tableclasificationsperround = table(factor(selection, levels = 1:length(unique(Motif))), factor(as.numeric(Motif[1:length(Motif)]), 1:length(unique(Motif))))
      clasificaciones = rbind(clasificaciones,
                              tableclasificationsperround[cbind(1:length(unique(Motif)), 1:length(unique(Motif)))] / table(as.numeric(Motif[1:length(Motif)])))
      motifnumero = unique(data.frame(as.numeric(Motif), Motif))
      colnames(tableclasificationsperround) <-
        motifnumero[match(colnames(tableclasificationsperround),
                          motifnumero$as.numeric.Motif.),]$Motif
      rownames(tableclasificationsperround) <-
        motifnumero[match(rownames(tableclasificationsperround),
                          motifnumero$as.numeric.Motif.),]$Motif
      
      tableclassifications  = abind(tableclassifications,
                                    tableclasificationsperround,
                                    along = 3)
      


colnames(clasificaciones) <-
  motifnumero[match(colnames(clasificaciones), motifnumero$as.numeric.Motif.),]$Motif

dev.off()


pdf("heatmapsclassificationLUSC_SJ.pdf")
for (cont in 1:dim(tableclassifications)[3]) {
  totalTP = sum(tableclassifications[cbind(2:length(unique(Motif)), 2:length(unique(Motif)), cont)]) / sum(table(as.numeric(Motif[1:length(Motif)]))[-1])
  heatmap.2(
    t(t(tableclassifications[, , cont]) / (colSums(
      tableclassifications[, , cont]
    ))),
    na.rm = T,
    Rowv = NULL,
    Colv = NULL,
    cexCol = 1,
    cexRow = 1,
    margins = c(12, 12),
    trace = "none",
    col = my_palette,
    ylab = "Predicted",
    xlab = "Observed",
    main = paste0("LUSC: \n Total ratio of TP = ", round(totalTP, digits = 2), "%")
  )
  
}

dev.off() 
write.table(
  tableclassifications,
  "tableclassificationsmutbySJLUSC_SJ.txt",
  sep = "\t",
  row.names = T,
  col.names = T,
  quote = F
)
