lungs = TCGAquery_recount2('TCGA', tissue = "lung")
classgenesoflungs = rownames(lungs[[1]])
datos = assay(lungs[[1]])
#CUIDADO CON ESTO!!! HA SIDO MUY A LO LOCO PERO ENCAJA
colnames(datos) <-  lungs$TCGA_lung$xml_bcr_patient_barcode
to_normalize <-  lungs$TCGA_lung$mapped_read_count
bplength = rowData(lungs[[1]])$bp_length
datosnormalizados = t(t(datos)/to_normalize)
datosnormalizados = datosnormalizados/bplength

ENSTLUAD = read.table("SJandENSTandexternamegenenameLUAD.txt", sep = "\t", stringsAsFactors = F, header = T)
uniqueENSTLUAD = unique(unlist(strsplit(ENSTLUAD$ensembl_gene_id_version, ",")))
selecciondatos = datosnormalizados[(classgenesoflungs %in% uniqueENSTLUAD), ]
write.table(selecciondatos,  "expressiondatanormalizedofgeneswithsignificantSJLUSC.txt", sep = "\t", col.names =  T, row.names =  T)

ENSTLUSC = read.table("SJandENSTandexternamegenenameLUSC.txt", sep = "\t", stringsAsFactors = F, header = T)
uniqueENSTLUSC = unique(unlist(strsplit(ENSTLUSC$ensembl_gene_id_version, ",")))
selecciondatos = datosnormalizados[(classgenesoflungs %in% uniqueENSTLUSC), ]
write.table(selecciondatos,  "expressiondatanormalizedofgeneswithsignificantSJLUAD.txt", sep = "\t", col.names =  T, row.names =  T)
