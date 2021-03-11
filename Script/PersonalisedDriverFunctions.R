#================================================================
#' This function allows you to convert miRNA names to version 21.
#' @param miRs miRNA names
#' @return Converted miRNA names
#================================================================
convertmiRs <- function(miRs) {
  # Get version of miRs
  version = checkMiRNAVersion(miRs, verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature = miRNA_PrecursorToMature(miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miRNA list
  # If a converted item is "NA", get the corresponding value from the mature list
  miRs <- miRNameList_Converted[, 2]
  naList <- which(is.na(miRNameList_Converted[, 2]))
  for(i in 1:length(naList)){
    k <- naList[i]
    miRs[k] <- miRMature[k, 2]
  }
  
  return(miRs)
}

#================================================================
#' miRNA target prediction with the Pearson correlation coefficient method, returning p-value
#' Calculate the Pearson correlation coefficient of each pair of miRNA-mRNA,and return a matrix of p-value of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' 
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @return A  matrix that includes the p-value of Pearson correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' 
#' @references
#' Pearson, K. (1920) Notes on the history of correlation. Biometrika, 13, 25 - 45.
#================================================================
Pearson_pValue=function(datacsv, cause, effect){
  data=Read(datacsv)
  data=scale(data) #standardise the data
  
  header<-readHeader(datacsv)
  num_miRNA<-length(cause)
  num_mRNA<-length(effect)
  miR<-header[1:num_miRNA]
  mR<-header[-(1:num_miRNA)]
  
  miRNA=data[,cause]
  mRNA=data[,effect]
  
  r <- matrix(data = NA, nrow = num_mRNA, ncol = num_miRNA)
  row.names(r) <- mR
  colnames(r) <- miR
  for (i in 1:num_mRNA) {
    for (j in 1:num_miRNA) {
      r[i,j] <- cor.test(mRNA[,i], miRNA[,j], method="pearson")$p.value
    }
  }
  
  return(r)
}

#================================================================
#' Adjust p-values
#================================================================
adjustpValues = function(results) {
  
  r <- results
  
  nR <- nrow(r)
  nC <- ncol(r)
  t <- as.vector(r)
  t <- p.adjust(t, method="fdr")
  r <- matrix(t, nrow = nR, ncol = nC)
  
  row.names(r) <- row.names(results)
  colnames(r) <- colnames(results)
  
  return(r)
}

#================================================================
#' Identify links among nodes from their expression data
#' @param data Expression data with rows being samples and columns being biological features
#' @param cause Range of cause
#' @param effect Range of effect
#' @param rootDir Root folder
#' @return Edges from cause to effect
#================================================================
identifyEdges = function(data, cause, effect, rootDir) {
  # Use Pearson to evaluate the relationship among cause and effect
  dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
  write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
  results = Pearson(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
  results <- t(results)
  
  # Get links which have the absolute coefficients more than
  # the average of all absolute coefficients
  a <- mean(abs(results))
  results <- abs(results) >= a
  ind <- which(results, arr.ind=TRUE)
  edges <- ind
  edges[, 1] <- row.names(ind)
  edges[, 2] <- colnames(results)[ind[, 2]]
  
  # Remove dataset.csv file
  file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  
  return(edges)
}

#================================================================
#' Identify links among nodes from their expression data
#' @param data Expression data with rows being samples and columns being biological features
#' @param cause Range of cause
#' @param effect Range of effect
#' @param rootDir Root folder
#' @param usingpVal TRUE if using p-value
#' @param cutoff FDR cutoff
#' @return Edges from cause to effect
#================================================================
identifyWeightedEdges = function(data, cause, effect, rootDir, usingpVal = FALSE, cutoff = 0.05) {
  if(usingpVal == TRUE) {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson_pValue(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Get links which have p-value < cutoff
    results <- adjustpValues(results)
    results <- results < cutoff
    ind <- which(results, arr.ind=TRUE)
    edges <- matrix(nrow = nrow(ind), ncol = 3)
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    n <- nrow(edges)
    for (i in 1:n) {
      edges[i, 3] <- results[ind[i, 1], ind[i, 2]]
    }
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  } else {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Get links which have the absolute coefficients more than
    # the average of all absolute coefficients
    a <- mean(abs(results))
    results[which(abs(results) < a)] <- 0
    ind <- which(results != 0, arr.ind=TRUE)
    edges <- matrix(nrow = nrow(ind), ncol = 3)
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    n <- nrow(edges)
    for (i in 1:n) {
      edges[i, 3] <- results[ind[i, 1], ind[i, 2]]
    }
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))  
  }
  
  return(edges)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomiR Number of miRNAs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order miRNAs, mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkWithmiRs = function(interactions, nomiR, nomR, noTF, data, rootDir){
  # Build network
  # miRNA => TF
  edges_non_cod <- identifyEdges(data, 1:nomiR, (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir)
  # miRNA => mRNA
  edges_non_cod <- rbind(edges_non_cod,
                         identifyEdges(data, 1:nomiR, (nomiR+1):(nomiR + nomR), rootDir))
  edges_non_cod[, 1] <- gsub("\\.", "-", edges_non_cod[, 1])
  # TF => miRNA
  edges_cod_non <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF), 1:nomiR, rootDir)
  edges_cod_non[, 2] <- gsub("\\.", "-", edges_cod_non[, 2])
  # TF => mRNA
  edges_cod_cod <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                                 (nomiR+1):(nomiR + nomR), rootDir)
  # mRNA => mRNA
  temp <- identifyEdges(data, (nomiR+1):(nomiR + nomR), (nomiR+1):(nomiR + nomR), rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  temp <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                        (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_non_cod) <- c("cause", "effect")
  colnames(edges_cod_non) <- c("cause", "effect")
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  # miRNA => TF & miRNA => mRNA: Filter with
  # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
  confirmedList <- read.csv(paste(rootDir, "/Data/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                  sep = ""), header = FALSE)
  colnames(confirmedList) <- c("cause", "effect")
  predictedList <- read.csv(paste(rootDir, "/Data/TargetScan_7.0.csv", sep = ""), header = TRUE)
  predictedList <- predictedList[, -2]
  colnames(predictedList) <- c("cause", "effect")
  confirmedList[, 1] <- convertmiRs(confirmedList[, 1])
  predictedList[, 1] <- convertmiRs(predictedList[, 1])
  targetbinding <- rbind(confirmedList, predictedList)
  targetbinding <- unique(targetbinding)
  targetbinding$effect <- as.character(targetbinding$effect)
  targetbinding <- targetbinding[which(targetbinding$cause %in% colnames(data)[1:nomiR]),]
  targetbinding <- targetbinding[
    which(targetbinding$effect %in% colnames(data)[(nomiR + 1):(nomiR+nomR+noTF)]),]
  targetbinding <- merge(targetbinding, edges_non_cod)
  # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
  interacts <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
                          sep = '\t', header = FALSE)
  interacts <- interacts[, 1:2]
  colnames(interacts) <- c("cause", "effect")
  interacts[, 2] <- convertmiRs(interacts[, 2])
  interacts <- unique(interacts)
  interacts$cause <- as.character(interacts$cause)
  interacts <- interacts[
    which(interacts$cause %in% colnames(data)[(nomiR+nomR+1):(nomiR+nomR+noTF)]),]
  interacts <- interacts[which(interacts$effect %in% colnames(data)[1:nomiR]),]
  interacts <- merge(interacts, edges_cod_non)
  
  # Combine
  interactions <- rbind(interactions, targetbinding)
  interactions <- rbind(interactions, interacts)
  
  return(interactions)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomiR Number of miRNAs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order miRNAs, mRNAs, TFs
#' @param rootDir Root folder
#' @param usingpVal TRUE if using p-value
#' @param cutoff FDR cutoff
#' @return Edges of the network from cause to effect
#================================================================
buildWeightedNetworkWithmiRs = function(interactions, nomiR, nomR, noTF, data, rootDir, usingpVal = FALSE, cutoff = 0.05){
  # Build network
  # miRNA => TF
  edges_non_cod <- identifyWeightedEdges(data, 1:nomiR, (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir, usingpVal, cutoff)
  # miRNA => mRNA
  edges_non_cod <- rbind(edges_non_cod,
                         identifyWeightedEdges(data, 1:nomiR, (nomiR+1):(nomiR + nomR), rootDir, usingpVal, cutoff))
  edges_non_cod[, 1] <- gsub("\\.", "-", edges_non_cod[, 1])
  # TF => miRNA
  edges_cod_non <- identifyWeightedEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF), 1:nomiR, rootDir, usingpVal, cutoff)
  edges_cod_non[, 2] <- gsub("\\.", "-", edges_cod_non[, 2])
  # TF => mRNA
  edges_cod_cod <- identifyWeightedEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                                 (nomiR+1):(nomiR + nomR), rootDir, usingpVal, cutoff)
  # mRNA => mRNA
  temp <- identifyWeightedEdges(data, (nomiR+1):(nomiR + nomR), (nomiR+1):(nomiR + nomR), rootDir, usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  temp <- identifyWeightedEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                        (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir, usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_non_cod) <- c("cause", "effect", "weight")
  colnames(edges_cod_non) <- c("cause", "effect", "weight")
  colnames(edges_cod_cod) <- c("cause", "effect", "weight")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  # miRNA => TF & miRNA => mRNA: Filter with
  # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
  confirmedList <- read.csv(paste(rootDir, "/Data/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                  sep = ""), header = FALSE)
  colnames(confirmedList) <- c("cause", "effect")
  predictedList <- read.csv(paste(rootDir, "/Data/TargetScan_7.0.csv", sep = ""), header = TRUE)
  predictedList <- predictedList[, -2]
  colnames(predictedList) <- c("cause", "effect")
  confirmedList[, 1] <- convertmiRs(confirmedList[, 1])
  predictedList[, 1] <- convertmiRs(predictedList[, 1])
  targetbinding <- rbind(confirmedList, predictedList)
  targetbinding <- unique(targetbinding)
  targetbinding$effect <- as.character(targetbinding$effect)
  targetbinding <- targetbinding[which(targetbinding$cause %in% colnames(data)[1:nomiR]),]
  targetbinding <- targetbinding[
    which(targetbinding$effect %in% colnames(data)[(nomiR + 1):(nomiR+nomR+noTF)]),]
  targetbinding <- merge(targetbinding, edges_non_cod)
  # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
  interacts <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
                          sep = '\t', header = FALSE)
  interacts <- interacts[, 1:2]
  colnames(interacts) <- c("cause", "effect")
  interacts[, 2] <- convertmiRs(interacts[, 2])
  interacts <- unique(interacts)
  interacts$cause <- as.character(interacts$cause)
  interacts <- interacts[
    which(interacts$cause %in% colnames(data)[(nomiR+nomR+1):(nomiR+nomR+noTF)]),]
  interacts <- interacts[which(interacts$effect %in% colnames(data)[1:nomiR]),]
  interacts <- merge(interacts, edges_cod_non)
  
  # Combine
  interactions <- rbind(interactions, targetbinding)
  interactions <- rbind(interactions, interacts)
  
  return(interactions)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkForPersonalised = function(interactions, nomR, noTF, data, rootDir){
  # Build network
  # TF => mRNA
  edges_cod_cod <- identifyEdges(data, (nomR+1):(nomR + noTF),
                                 1:nomR, rootDir)
  # mRNA => mRNA
  temp <- identifyEdges(data, 1:nomR, 1:nomR, rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  temp <- identifyEdges(data, (nomR+1):(nomR + noTF),
                        (nomR+1):(nomR + noTF), rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  
  return(interactions)
}

#================================================================
#' Analyse a network
#' @param nomiR Number of miRNAs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param network Edges of the network
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order miRNAs, mRNAs, TFs
#' @param fileName File to be written
#================================================================
analyseNetwork = function(nomiR, nomR, noTF, network, data, fileName){
  # Prepare the data
  dat <- ""
  dat <- paste(dat, "Number of nodes:", nomiR + nomR + noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNAs:", nomiR, "\n", sep = "\t")
  dat <- paste(dat, "Number of TFs:", noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNAs:", nomR, "\n", sep = "\t")
  dat <- paste(dat, "\n", sep = "\t")
  miRs <- colnames(data)[1:nomiR]
  mRs <- colnames(data)[(nomiR+1):(nomiR+nomR)]
  TFs <- colnames(data)[(nomiR+nomR+1):(nomiR+nomR+noTF)]
  edge_type1 <- network[which(network$cause %in% miRs),]
  edge_type1 <- edge_type1[which(edge_type1$effect %in% TFs),]
  edge_type1 <- nrow(edge_type1)
  edge_type2 <- network[which(network$cause %in% miRs),]
  edge_type2 <- edge_type2[which(edge_type2$effect %in% mRs),]
  edge_type2 <- nrow(edge_type2)
  edge_type3 <- network[which(network$cause %in% TFs),]
  edge_type3 <- edge_type3[which(edge_type3$effect %in% miRs),]
  edge_type3 <- nrow(edge_type3)
  edge_type4 <- network[which(network$cause %in% TFs),]
  edge_type4 <- edge_type4[which(edge_type4$effect %in% TFs),]
  edge_type4 <- nrow(edge_type4)
  edge_type5 <- network[which(network$cause %in% TFs),]
  edge_type5 <- edge_type5[which(edge_type5$effect %in% mRs),]
  edge_type5 <- nrow(edge_type5)
  edge_type6 <- network[which(network$cause %in% mRs),]
  edge_type6 <- edge_type6[which(edge_type6$effect %in% mRs),]
  edge_type6 <- nrow(edge_type6)
  edge <- edge_type1 + edge_type2 + edge_type3 + edge_type4 + edge_type5 + edge_type6
  dat <- paste(dat, "Number of edges:", edge, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNA-TF edges:", edge_type1, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNA-mRNA edges:", edge_type2, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-miRNA edges:", edge_type3, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-TF edges:", edge_type4, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-mRNA edges:", edge_type5, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNA-mRNA edges:", edge_type6, "\n", sep = "\t")
  
  # Write file
  writeLines(dat, fileName)
}

#================================================================
#' Analyse a network
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param network Edges of the network
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param fileName File to be written
#================================================================
analyseNetworkForPersonalised = function(nomR, noTF, network, data, fileName){
  # Prepare the data
  dat <- ""
  dat <- paste(dat, "Number of nodes:", nomR + noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of TFs:", noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNAs:", nomR, "\n", sep = "\t")
  dat <- paste(dat, "\n", sep = "\t")
  mRs <- colnames(data)[1:nomR]
  TFs <- colnames(data)[(nomR+1):(nomR+noTF)]
  edge_type4 <- network[which(network$cause %in% TFs),]
  edge_type4 <- edge_type4[which(edge_type4$effect %in% TFs),]
  edge_type4 <- nrow(edge_type4)
  edge_type5 <- network[which(network$cause %in% TFs),]
  edge_type5 <- edge_type5[which(edge_type5$effect %in% mRs),]
  edge_type5 <- nrow(edge_type5)
  edge_type6 <- network[which(network$cause %in% mRs),]
  edge_type6 <- edge_type6[which(edge_type6$effect %in% mRs),]
  edge_type6 <- nrow(edge_type6)
  edge <- edge_type4 + edge_type5 + edge_type6
  dat <- paste(dat, "Number of edges:", edge, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-TF edges:", edge_type4, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-mRNA edges:", edge_type5, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNA-mRNA edges:", edge_type6, "\n", sep = "\t")
  
  # Write file
  writeLines(dat, fileName)
}

#================================================================
#' Analyse controllability of the network
#' @param result Result file of analysing controllability of the network
#' @param outFile File to be written
#================================================================
analyseControllability = function(result, outFile){
  # Read files
  result_output <- read.csv(result, header = FALSE)
  # Prepare the data
  items <- c("1: # of nodes having links", "2: # of edges", "3: average degree",
             "4: # of driver nodes", "5: fraction of driver nodes = Nd/N",
             "6: fraction of type-I critical nodes",
             "7: fraction of type-I redundant nodes", "8: fraction of type-I ordinary nodes",
             "9: fraction of type-II critical nodes",
             "10: fraction of type-II redundant nodes", "11: fraction of type-II ordinary nodes",
             "12: fraction of critical links",
             "13: fraction of redundant links", "14: fraction of ordinary links")
  dat <- ""
  for (i in 1:14) {
    dat <- paste(dat, items[i], result_output[1,i], "\n", sep = "\t")  
    if(i %in% c(3,5,8,11)) {
      dat <- paste(dat, "\n")
    }
  }
  
  # Write file
  writeLines(dat, outFile)
}

#================================================================
#' Convert a list to a string
#' @param aList A list
#' @return A string contains items of the list
#================================================================
getList = function(aList){
  len <- length(aList)
  l <- aList[1]
  if(len > 1) {
    for(i in 2:len) {
      l <- paste(l, aList[i], sep=",")
    }
  }
  
  return(l)
}

#================================================================
#' Print gene lists for finding overlaps among methods
#' @param OncodriveCLUST A list of genes from OncodriveCLUST method
#' @param OncodriveFM A list of genes from OncodriveFM method
#' @param DawnRank A list of genes from DawnRank method
#' @param CBNA A list of genes from CBNA method
#================================================================
printGeneList = function(OncodriveCLUST, OncodriveFM,
                         DawnRank, CBNA) {
  print("OncodriveCLUST")
  print(getList(OncodriveCLUST))
  print("OncodriveFM")
  print(getList(OncodriveFM))
  print("DawnRank")
  print(getList(DawnRank))
  print("CBNA")
  print(getList(CBNA))  
}

#================================================================
#' Process data by filtering out genes with constant expression values
#' @param matchedData Matched expression data, rownames are samples and colnames are biological features
#' @return Processed data
#================================================================
processData <- function(matchedData) {
  # Convert input matched expression data to rownames being genes and colnames being samples
  miR <- matchedData$miRs
  mRNA <- matchedData$mRNAs
  Exp = t(cbind(miR, mRNA))
  
  # Filter out genes with constant expression values
  sdGenes <- apply(Exp, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    Exp <- Exp[sdGenes > 0 & !is.na(sdGenes), , drop=FALSE]
  }
  
  return(Exp)
}

#================================================================
#' Classify samples into Epi, Mes, Intermediate Epi, or Intermediate Mes
#' @param Generic_EMT_signature EMT signatures
#' @param Exp Expression data
#' @return Classified samples
#================================================================
scoreEMT <- function(Generic_EMT_signature, Exp) {
  # Get EMT genes, Epi genes, and Mes genes
  EMT_signature = lapply(1:dim(Generic_EMT_signature)[1], function(i) Generic_EMT_signature[i,1])
  Epi = Generic_EMT_signature[which(Generic_EMT_signature[, 2] %in% "Epi"), 1]
  Mes = Generic_EMT_signature[which(Generic_EMT_signature[, 2] %in% "Mes"), 1]
  
  # Only get genes which have expression data
  Epi_Update = Epi[which(Epi %in% rownames(Exp))]
  Mes_Update = Mes[which(Mes %in% rownames(Exp))]
  
  # Estimate GSVA enrichment scores (ES) for two lists of EMT signatures: Epi and Mes.
  # The method is gsva (Hanzelmann et al., 2013), ES is calculated as the maximum distance of the random walk from 0.
  # gsva_es_list = gsva(Exp, list(Epi_Update, Mes_Update), mx.diff=0)$es.obs
  gsva_es_list = gsva(Exp, list(Epi_Update, Mes_Update), mx.diff=0)
  
  # The difference of GSVA enrichment scores between Epi and Mes
  EMT_score = abs(gsva_es_list[2,]) - abs(gsva_es_list[1,])
  
  # Estimate GSVA enrichment score for each gene of EMT signatures.
  # The method is gsva (Hanzelmann et al., 2013), ES is calculated as the maximum distance of the random walk from 0.
  # gsva_es_signature = gsva(Exp, EMT_signature, mx.diff=0)$es.obs
  gsva_es_signature = gsva(Exp, EMT_signature, mx.diff=0)
  
  # Using two-sample KS test, compute the p-values of the difference of GSVA enrichment scores between Epi and Mes
  Epi_index = 1:length(Epi_Update)
  Mes_index = (length(Epi_Update)+1):(length(Epi_Update)+length(Mes_Update))
  EMT_p.value = unlist(lapply(1:dim(gsva_es_signature)[2], function(i) ks.test(gsva_es_signature[Epi_index,i], gsva_es_signature[Mes_index,i])$p.value))
  
  # We divide the samples to four types:
  # Mes (EMT_p.value<0.05 & EMT_score>0),
  # Intermediate Mes (EMT_p.value>=0.05 & EMT_score>0),
  # Epi (EMT_p.value<0.05 & EMT_score<0),
  # and Intermediate Epi (EMT_p.value>=0.05 & EMT_score<0)
  Sample_Phenotype=c()
  for(i in 1:dim(gsva_es_signature)[2]){
    if(EMT_p.value[i]<0.05 & EMT_score[i]>0){
      Sample_Phenotype[i]="Mes"}
    else if(EMT_p.value[i]>=0.05 & EMT_score[i]>0){
      Sample_Phenotype[i]="Intermediate Mes"}
    else if(EMT_p.value[i]<0.05 & EMT_score[i]<0){
      Sample_Phenotype[i]="Epi"}
    else if(EMT_p.value[i]>=0.05 & EMT_score[i]<0){
      Sample_Phenotype[i]="Intermediate Epi"}
  }
  Sample_EMT_Score=cbind(EMT_score,EMT_p.value,Sample_Phenotype)
  
  return(Sample_EMT_Score)
}

#================================================================
#' This function allows you to prepare data for classifying cancer subtypes
#' @param d a matrix of expression of 50 mRNAs in Pam50 with rows being samples and columns being mRNA names.
#' @param directoryPath the directory path to save result file.
#' @param dataDir the directory for saving the result data.
#' @param inputFileName the file name of the input data which is used to classify cancer subtypes.
#================================================================
prepareData=function(d, directoryPath, dataDir = "bioclassifier_data", inputFileName = "inputFile.txt"){
  d <- t(d)
  noOfRow <- nrow(d)
  noOfCol <- ncol(d)
  temp <- matrix(0, nrow = (noOfRow + 1), ncol = noOfCol)
  row.names(temp) <- c("T", row.names(d))
  colnames(temp) <- colnames(d)
  temp[2:(noOfRow + 1), 1:noOfCol] <- d[,]
  write.table(temp, paste(directoryPath, "/", dataDir, "/", inputFileName, sep = ""), col.names=NA, sep = "\t")
}

#================================================================
#' This function allows you to compute measures (Precision, Recall, and F1 Score)
#' @param noOfConfirmed the number of found genes in CGC.
#' @param noOfNotConfirmed the number of found genes not in CGC.
#' @param noOfCGC the number of genes in CGC.
#' @return A list including Precision, Recall, and F1 Score
#================================================================
computeMeasures=function(noOfConfirmed, noOfNotConfirmed, noOfCGC){
  p <- noOfConfirmed/(noOfConfirmed + noOfNotConfirmed)
  r <- noOfConfirmed/noOfCGC
  f1 <- 2*((p*r)/(p+r))
  
  l <- list(Precision = p, Recall = r, F1Score = f1)
  
  return(l)
}

#================================================================
#' This function allows you to prepare data for comparing methods
#' @param Found_Gene_List_Validated the validated found gene list.
#' @param Found_Gene_List the found gene list.
#' @param gold_standard the CGC list.
#' @return A list of genes with corresponding Precision, Recal, and F1 Score for top found genes
#================================================================
prepareDataForComparison=function(Found_Gene_List_Validated, Found_Gene_List, gold_standard){
  noOfCGC <- length(gold_standard)
  noOfFoundGenes <- length(Found_Gene_List)
  r <- matrix(nrow = noOfFoundGenes, ncol = 5)
  colnames(r) <- c("Gene", "InCGC", "Precision", "Recall", "F1Score")
  r[,1] <- Found_Gene_List
  r[,2] <- "No"
  r[which(r[,1] %in% Found_Gene_List_Validated),2] <- "Yes"
  for (i in 1:noOfFoundGenes) {
    topGenes <- r[1:i, , drop=FALSE]
    noOfConfirmed <- sum(topGenes[,2] == "Yes")
    noOfNotConfirmed <- i - noOfConfirmed
    t <- computeMeasures(noOfConfirmed, noOfNotConfirmed, noOfCGC)
    r[i,3] <- t$Precision
    r[i,4] <- t$Recall
    r[i,5] <- t$F1Score
  }
  
  return(r)
}

#================================================================
#' This function allows you to draw a line chart
#' @param OncodriveCLUSTData the data for OncodriveCLUST.
#' @param OncodriveFMData the data for OncodriveFM.
#' @param DawnRankData the data for DawnRank.
#' @param CBNAData the data for CBNA.
#' @param Type it can be "Precision", "Recall", or "F1Score.
#================================================================
drawLineChart=function(OncodriveCLUSTData, OncodriveFMData, DawnRankData, CBNAData, Type="Precision"){
  # Parameters for each chart type
  if(Type == "Precision") {
    yText <- "Precision according to CGC"
    t <- "Precision Comparison"
    pos <- "topright"
  } else if (Type == "Recall") {
    yText <- "Recall according to CGC"
    t <- "Recall Comparison"
    pos <- "topleft"
  } else { # "F1Score"
    yText <- "F1 Score according to CGC"
    t <- "F1 Score Comparison"
    pos <- "topleft"
  }
  
  # Prepare data
  Method=c(rep("OncodriveCLUST", 200), rep("OncodriveFM", 200), rep("DawnRank", 200), rep("CBNA", 200))
  TopNGenes=c(1:200, 1:200, 1:200, 1:200)
  Val=c(OncodriveCLUSTData, OncodriveFMData, DawnRankData, CBNAData)
  indexes<-seq(5,800,5)
  Method <- Method[indexes]
  TopNGenes <- TopNGenes[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
            TopNGenes=TopNGenes,
            Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("OncodriveCLUST", "OncodriveFM", "DawnRank", "CBNA")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$TopNGenes) 
  yrange <- range(data$Val) 
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Top N genes",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], TopNGenes=data$TopNGenes[ind], Val=data$Val[ind])
    lines(method$TopNGenes, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
         pch=plotchar, lty=linetype, title="Method")
}

#================================================================
#' This function allows you to draw a line chart
#================================================================
drawLineChart3=function(pr, im, nc, Type="Precision"){
  # Parameters for each chart type
  if(Type == "Precision") {
    yText <- "Precision"
    t <- ""
    pos <- "topright"
  } else if (Type == "Recall") {
    yText <- "Recall"
    t <- ""
    pos <- "topleft"
  } else { # "F1Score"
    yText <- "F1 Score"
    t <- ""
    pos <- "topleft"
  }
  
  # Prepare data
  Method=c(rep("PageRank", 200), rep("Influence_Maximisation", 200), rep("Network_Control", 200))
  TopNGenes=c(1:200, 1:200, 1:200)
  Val=c(pr, im, nc)
  indexes<-seq(5,600,5)
  Method <- Method[indexes]
  TopNGenes <- TopNGenes[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
               TopNGenes=TopNGenes,
               Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("PageRank", "Influence_Maximisation", "Network_Control")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$TopNGenes) 
  yrange <- range(data$Val) 
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Top N genes",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], TopNGenes=data$TopNGenes[ind], Val=data$Val[ind])
    lines(method$TopNGenes, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  # title(t, cex.main=2.5)
  
  # Add a legend 
  legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
         pch=plotchar, lty=linetype, title="", bty = "n")
}

#================================================================
#' This function allows you to draw a line chart
#' @param dat the data for the chart.
#' @param ind index of the network.
#' @param hasyText if display yText.
#' @param hasLegend if display legend.
#================================================================
drawLineChartForEachNetwork=function(dat, ind, hasyText, hasLegend){
  # Parameters for each chart type
  if(hasyText) {
    yText <- "No. of validated genes"
  } else {
    yText <- ""
  }
  t <- paste("Network ", ind, sep = "")
  pos <- "topleft"
  
  # Prepare data
  Method=c(rep("GoldStandard", 20), rep("Random", 20), rep("jointIDA", 20), rep("GroupDriver", 20))
  NCases=c(1:20, 1:20, 1:20, 1:20)
  Val=c(dat[,2], dat[,3], dat[,4], dat[,5])
  indexes<-seq(1,80,1)
  Method <- Method[indexes]
  NCases <- NCases[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
               NCases=NCases,
               Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("GoldStandard", "Random", "jointIDA", "GroupDriver")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$NCases) 
  # yrange <- range(data$Val) 
  yrange <- c(0,40)
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Case",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], NCases=data$NCases[ind], Val=data$Val[ind])
    lines(method$NCases, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  if(hasLegend) {
    legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
           pch=plotchar, lty=linetype, title="", bty = "n")
  }
}

#================================================================
#' This function allows you to draw a line chart
#' @param dat the data for the chart.
#' @param ind index of the network.
#' @param hasyText if display yText.
#' @param hasLegend if display legend.
#' @param yRange Range of y.
#================================================================
drawLineChartForEachNetwork2=function(dat, ind, hasyText, hasLegend, yRange=40){
  # Parameters for each chart type
  if(hasyText) {
    yText <- "No. of validated genes"
  } else {
    yText <- ""
  }
  t <- paste("Network ", ind, sep = "")
  pos <- "topleft"
  
  # Prepare data
  Method=c(rep("GoldStandard", 20), rep("Random", 20), rep("jointIDA", 20), rep("GroupDriver", 20))
  NCases=c(1:20, 1:20, 1:20, 1:20)
  Val=c(dat[,2], dat[,3], dat[,4], dat[,5])
  indexes<-seq(1,80,1)
  Method <- Method[indexes]
  NCases <- NCases[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
               NCases=NCases,
               Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("GoldStandard", "Random", "jointIDA", "GroupDriver")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$NCases) 
  # yrange <- range(data$Val) 
  yrange <- c(0,yRange)
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Number of cases",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], NCases=data$NCases[ind], Val=data$Val[ind])
    lines(method$NCases, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  if(hasLegend) {
    legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
           pch=plotchar, lty=linetype, title="", bty = "n")
  }
}

#================================================================
#' This function allows you to prepare data for a clustergram
#' @param dat the data to draw.
#' @param termTop the number of top terms to draw.
#' @param geneTop the number of top genes to draw.
#' @return Data for draw clustergram
#================================================================
prepareDataForClustergram=function(dat, termTop = 10, geneTop  = 20){
  d <- dat
  d <- d[1:termTop,]
  d[,1] <- substr(d[,1], nchar(d[,1])-10, nchar(d[,1])-1)
  genes <- NULL
  geneList <- NULL
  for (i in 1:termTop) {
    genes <- c(genes, unlist(strsplit(d[i,4], ";")))
    geneList <- c(geneList, strsplit(d[i,4], ";"))
  }
  genes <- unique(genes)
  nGenes <- length(genes)
  r <- matrix(0, nrow = termTop+1, ncol = nGenes)
  row.names(r) <- c(d[,1], "total")
  colnames(r) <- genes
  for (i in 1:termTop) {
    ind <- which(colnames(r)[] %in% unlist(geneList[i]))
    r[i,ind] <- 1
  }
  for (i in 1:nGenes) {
    r[termTop+1,i] <- sum(r[1:termTop,i])
  }
  r <- r[,order(r[termTop+1,], decreasing = FALSE)]
  r <- r[1:termTop, (nGenes-geneTop + 1):nGenes]
  r <- t(r)
  newData <- cbind("", r)
  newData[,1] <- rownames(r)
  colnames(newData)[1] <- "Gene"
  
  return(newData)
}

#================================================================
#' This function allows you to draw a clustergram
#' @param dat the data to draw.
#' @param t title
#================================================================
drawClustergram=function(dat, t){
  d <- dat
  d.m <- melt(d)
  d.m <- ddply(d.m, .(variable), transform, rescale = rescale(value))
  p <- ggplot(d.m, aes(variable, Gene)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient(low = "#CCCCCC", high = "#E85642")
  base_size <- 15
  p + theme_grey(base_size = base_size) + 
    labs(title= t, x = "Enriched Terms", y = "Predicted Cancer Drivers") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(legend.position = "none", axis.ticks = element_blank(), 
          axis.text.x = element_text(size = base_size*1.5, 
                                     angle = 285, hjust = 0,
                                     colour = "black"),
          axis.text.y = element_text(size = base_size*1.5,
                                     colour = "black"),
          axis.title = element_text(size=base_size*2),
          plot.title = element_text(size=base_size*2.5))
}

#================================================================
#' This function allows you to read results of exploring drivers for cancer subtypes
#' @param resultDir the directory to get data.
#' @param type the cancer subtype.
#' @return Results of cancer subtype drivers
#================================================================
readResult = function(resultDir, type) {
  # Get the folder name
  outDir <- paste(resultDir, "/", type, sep = "")
  
  # Coding with mutations
  coding_mutations = read.csv(
    file = paste(outDir, "/coding_candidate_cancer_drivers_mutations.csv", sep = ""), as.is = TRUE)
  
  # Coding without mutations
  coding_no_mutations = read.csv(
    file = paste(outDir, "/coding_candidate_cancer_drivers_no_mutations.csv", sep = ""), as.is = TRUE)
  
  # Noncoding
  noncoding = read.csv(
    file = paste(outDir, "/noncoding_candidate_cancer_drivers.csv", sep = ""), as.is = TRUE)
  
  # Combine
  all <- rbind(coding_mutations, coding_no_mutations, noncoding)
  l = list(subType = type, coding_mutations = coding_mutations, coding_no_mutations = coding_no_mutations, noncoding = noncoding, all = all)
  
  return(l)
}

#================================================================
#' This function allows you to get drivers which are specific for a subtype
#' @param basalList Basal drivers.
#' @param her2List Her2 drivers.
#' @param lumAList LumA drivers.
#' @param lumBList LumB drivers.
#' @param normalLikeList Normal-like drivers.
#' @return Different drivers of subtypes
#================================================================
findDiff = function(basalList, her2List, lumAList, lumBList, normalLikeList) {
  basalDiff <- setdiff(basalList$Name, c(her2List$Name, lumAList$Name, lumBList$Name, normalLikeList$Name))
  her2Diff <- setdiff(her2List$Name, c(basalList$Name, lumAList$Name, lumBList$Name, normalLikeList$Name))
  lumADiff <- setdiff(lumAList$Name, c(basalList$Name, her2List$Name, lumBList$Name, normalLikeList$Name))
  lumBDiff <- setdiff(lumBList$Name, c(basalList$Name, her2List$Name, lumAList$Name, normalLikeList$Name))
  normalLikeDiff <- setdiff(normalLikeList$Name, c(basalList$Name, her2List$Name, lumAList$Name, lumBList$Name))
  
  l = list(basalDiff = basalDiff, her2Diff = her2Diff, lumADiff = lumADiff, lumBDiff = lumBDiff, normalLikeDiff = normalLikeDiff)
  
  return(l)
}

#================================================================
#' This function allows you to validate with the CGC
#' @param l List of drivers.
#' @param gold_standard The CGC.
#' @return Validated cancer drivers
#================================================================
validateCGC=function(l, gold_standard) {
  Top50 <- l$coding_mutations[1:50,]
  Top100 <- l$coding_mutations[1:100,]
  Top150 <- l$coding_mutations[1:150,]
  Top200 <- l$coding_mutations[1:200,]
  Top50_Validated <- intersect(Top50[, "Name"], gold_standard)
  Top100_Validated <- intersect(Top100[, "Name"], gold_standard)
  Top150_Validated <- intersect(Top150[, "Name"], gold_standard)
  Top200_Validated <- intersect(Top200[, "Name"], gold_standard)
  
  r <- list(Top50_Validated=Top50_Validated, Top100_Validated=Top100_Validated,
            Top150_Validated=Top150_Validated, Top200_Validated=Top200_Validated)
  
  return(r)
}

#================================================================
#' This function allows you to process miRs of Mes
#' Refer to http://www.mirbase.org/help/nomenclature.shtml for more information
#' @param noncoding_gold_standard_Mes List of miRs.
#' @return Processed miRs
#================================================================
processmiR=function(noncoding_gold_standard_Mes) {
  # Put data to returned variable
  r <- noncoding_gold_standard_Mes
  
  # Get miRs which do not start with "miR"
  ind <- which(substr(r[],1,3) != "miR")
  
  # Add prefix
  r[-ind] <- paste("hsa-", r[-ind], sep = "")
  
  # Get version of miRs
  miRs <- r
  version = checkMiRNAVersion(miRs, verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature = miRNA_PrecursorToMature(miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miRNA list
  # If a converted item is "NA", get the corresponding value from the original
  miRs <- miRNameList_Converted[, 2]
  naList <- which(is.na(miRNameList_Converted[, 2]))
  for(i in 1:length(naList)){
    k <- naList[i]
    miRs[k] <- miRMature[k, 1]
  }
  
  return(miRs)
}

#================================================================
#' This function allows you to compute node weight
#' @param cancer_data Cancer expression data.
#' @param normal_data Normal expression data.
#' @return A list of weighted nodes
#================================================================
computeNodeWeight=function(cancer_data, normal_data) {
  cancer <- colMeans(cancer_data)
  normal <- colMeans(normal_data)
  nodeList <- cancer
  nodeList[] <- abs(normal[] - cancer[])
  
  return(nodeList)
}

#================================================================
#' This function allows you to evaluate influence of a node set on another node set
#' @param cause Node set of cause.
#' @param effect Node set of effect.
#' @param network Network data.
#' @param nodes Nodes with weight.
#' @param getGeneName If TRUE, return gene list, if FALSE, return number of genes.
#' @return Influence of cause on effect
#================================================================
evaluateInfluence=function(cause, effect, network, nodes, getGeneName = FALSE) {
  n_nodes <- nrow(nodes)
  n_activeNodes <- 0
  activeNodes <- cause[[1]]
  n_newActiveNodes <- length(activeNodes)
  while ((n_newActiveNodes > n_activeNodes) & (n_newActiveNodes < n_nodes)) {
    n_activeNodes <- length(activeNodes)
    newActiveNodes <- NULL
    inactiveNodes <- nodes[which(!(nodes[,1] %in% activeNodes)),1]
    n_inactiveNodes <- length(inactiveNodes)
    for (i in 1:n_inactiveNodes) {
      curNode <- inactiveNodes[i]
      edges <- network[which(network$effect %in% curNode),]
      edges <- edges[which(edges$cause %in% activeNodes),]
      n_edges <- nrow(edges)
      if(n_edges > 0) {
        if(abs(sum(edges[,3])) >= nodes[which(nodes[,1] == curNode),2]) {
          newActiveNodes <- union(newActiveNodes, curNode)
        }
      }
    }
    activeNodes <- union(activeNodes, newActiveNodes)
    n_newActiveNodes <- length(activeNodes)
  }
  if (getGeneName == TRUE) {
    r <- intersect(activeNodes, effect[,1])
  } else {
    r <- length(intersect(activeNodes, effect[,1]))  
  }
  
  return(r)
}

#================================================================
#' This function allows you to get influence of a node set on another node set
#' @param cause Node set of cause.
#' @param effect Node set of effect.
#' @param network Network data.
#' @param nodes Nodes with weight.
#' @return Influence of cause on effect
#================================================================
getInfluence=function(cause, effect, network, nodes) {
  n_nodes <- nrow(nodes)
  n_activeNodes <- 0
  activeNodes <- cause[[1]]
  n_newActiveNodes <- length(activeNodes)
  while ((n_newActiveNodes > n_activeNodes) & (n_newActiveNodes < n_nodes)) {
    n_activeNodes <- length(activeNodes)
    newActiveNodes <- NULL
    inactiveNodes <- nodes[which(!(nodes[,1] %in% activeNodes)),1]
    n_inactiveNodes <- length(inactiveNodes)
    for (i in 1:n_inactiveNodes) {
      curNode <- inactiveNodes[i]
      edges <- network[which(network$effect %in% curNode),]
      edges <- edges[which(edges$cause %in% activeNodes),]
      n_edges <- nrow(edges)
      if(n_edges > 0) {
        if(abs(sum(edges[,3])) >= nodes[which(nodes[,1] == curNode),2]) {
          newActiveNodes <- union(newActiveNodes, curNode)
        }
      }
    }
    activeNodes <- union(activeNodes, newActiveNodes)
    n_newActiveNodes <- length(activeNodes)
  }
  r <- intersect(activeNodes, effect[,1])
  
  return(r)
}

#================================================================
#' getInfluenceForDREAMData
#================================================================
getInfluenceForDREAMData=function(rootDir, i_net, nodeWeight) {
  
  # cause
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net,
                    "_dualknockouts_indexes.tsv", sep = "")
  c <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  n_cause <- nrow(c)
  
  # effect
  effect <- matrix(1:100, nrow = 100, ncol = 1)
  
  # network
  filename <- paste(rootDir,
                    "/Data/DREAM4_InSilicoNetworks_GoldStandard/DREAM4_Challenge2_GoldStandards/Size 100/DREAM4_GoldStandard_InSilico_Size100_",
                    i_net, ".tsv", sep = "")
  net <- read.table(filename, header=FALSE, sep="\t", as.is = TRUE)
  net <- net[which(net$V3 == 1),]
  net <- net[,c(1:2)]
  net[,1] <- gsub("[G]", "", net[,1])
  net[,2] <- gsub("[G]", "", net[,2])
  colnames(net) <- c("cause", "effect")
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net, "_knockouts.tsv", sep = "")
  net_data <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  net_data <- as.matrix(net_data)
  dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
  write.csv(net_data[, c(1:100, 1:100)], dataset, row.names = FALSE)
  results = Pearson(dataset, 1:100, 101:200)
  results <- t(results)
  ind <- which(results != 0, arr.ind=TRUE)
  edges <- matrix(nrow = nrow(ind), ncol = 3)
  edges[, 1] <- row.names(ind)
  edges[, 2] <- colnames(results)[ind[, 2]]
  n <- nrow(edges)
  for (i in 1:n) {
    edges[i, 3] <- results[ind[i, 1], ind[i, 2]]
  }
  file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  colnames(edges) <- c("cause", "effect", "weight")
  edges[,1] <- gsub("[G]", "", edges[,1])
  edges[,2] <- gsub("[G]", "", edges[,2])
  edges[,2] <- gsub("\\..*", "", edges[,2])
  network <- merge(net, edges)
  network$weight <- as.numeric(levels(network$weight))[network$weight]
  network <- normaliseEdgeWeight(network)
  
  # nodes
  nodes <- matrix(nrow = 100, ncol = 2)
  nodes[,1] <- 1:100
  nodes[,2] <- nodeWeight
  
  # get influence
  r <- list()
  for (i in 1:n_cause) {
    cause <- list(c(c[i,1], c[i,2]))
    inf <- getInfluence(cause, effect, network, nodes)
    r[[i]] <- inf
  }
  
  return(r)
}

#================================================================
#' This function allows you to normalise edge weight
#' @param network Network data.
#' @return Network with edge weight normalised
#================================================================
normaliseEdgeWeight=function(network) {
  targets <- network[,2]
  targets <- unique(targets)
  n_targets <- length(targets)
  maxWeight <- 0
  for (i in 1:n_targets) {
    weight <- sum(network[which(network$effect == targets[i]),3])
    if(abs(weight) > maxWeight) {
      maxWeight <- abs(weight)
    }
  }
  network[,3] <- network[,3]/maxWeight
  
  return(network)
}

#================================================================
#' This function allows you to normalise edge weight
#' @param network Network data.
#' @return Network with edge weight normalised
#================================================================
normaliseEdges=function(network) {
  targets <- network[,2]
  targets <- unique(targets)
  n_targets <- length(targets)
  for (i in 1:n_targets) {
    weight <- sum(network[which(network$effect == targets[i]),3])
    network[which(network$effect == targets[i]),3] <- network[which(network$effect == targets[i]),3]/abs(weight)
  }
  
  return(network)
}

#================================================================
#' This function allows you to select candidate groups
#' @param groups Candidate groups.
#' @return Indexes of selected groups
#================================================================
selectCandidates=function(groups) {
  groups <- groups[,1]
  n_groups <- length(groups)
  l <- NULL
  for (i in 1:n_groups) {
    l <- c(l, c(strsplit(groups[i], split=' ', fixed=TRUE)))
  }
  ind <- c(rep("Validated", n_groups))
  i <- n_groups
  while (i > 1) {
    for (j in 1:(i-1)) {
      if(all(l[[j]] %in% l[[i]])) {
        ind[i] <- "Invalidated"
        break
      }
    }
    i <- i-1
  }
  
  return(ind)
}

#================================================================
#' This function allows you to build the network for jointIda function
#' @param net Network data.
#' @return Network used for joinIda
#================================================================
prepareNetworkForJointIda <- function(net) {
  net <- net[,1:2]
  n_net <- nrow(net)
  nodes <- unique(union(net$cause, net$effect))
  n_nodes <- length(nodes)
  A <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  colnames(A) <- nodes
  row.names(A) <- nodes
  for (i in 1:n_net) {
    A[net$cause[i], net$effect[i]] <- 1
  }
  graph <- getGraph(A)
  # saveRDS(A, paste(outDir, "/network_jointIda.rds", sep = ""))
  # library(pcalg)
  # A <- readRDS(paste("R/pro/network_jointIda.rds", sep = ""))
  # isValidGraph(amat = A, type = "dag", verbose = TRUE)
  # isValidGraph(amat = A, type = "cpdag", verbose = TRUE)
  # isValidGraph(amat = A, type = "pdag", verbose = TRUE)
  return(graph)
}

#================================================================
#' This function allows you to build the network for jointIda function
#' @param net Network data.
#' @param nod Nodes.
#' @return Network used for joinIda
#================================================================
getNetworkForJointIda <- function(net, nod) {
  net <- net[,1:2]
  n_net <- nrow(net)
  n_nod <- length(nod)
  A <- matrix(0, nrow = n_nod, ncol = n_nod)
  colnames(A) <- nodes
  rownames(A) <- nodes
  for (i in 1:n_net) {
    A[net$cause[i], net$effect[i]] <- 1
  }
  graph <- getGraph(A)
  # isValidGraph(amat = A, type = "dag", verbose = TRUE)
  # isValidGraph(amat = A, type = "cpdag", verbose = TRUE)
  # isValidGraph(amat = A, type = "pdag", verbose = TRUE)
  return(graph)
}

#=======================================
#' This function allows you to select mRNAs
#' which have the most variant Median Absolute Deviation (MAD).
#' @param matchedData Expression data of mRNAs.
#' @param topk_mR Number of mRNAs to be cut off.
#' @return A list containing selected mRNAs. The list of elements includes:
#' \item{d} {The data with rows being samples and columns being mRNAs.}
#' \item{mRs} {The names of mRNAs.}
#=======================================
getDatabyMAD <- function(matchedData, topk_mR) {
  # Identify significant mRNAs by using function FSbyMAD in CancerSubtypes package
  mRNAsData = FSbyMAD(t(matchedData$mRNAs), value=topk_mR)
  mRNAsData <- t(mRNAsData)
  
  # Remove duplicated data
  mRNAsData <- mRNAsData[,!duplicated(colnames(mRNAsData))]
  
  # Prepare the result
  mRs <- colnames(mRNAsData)
  l = list(d = mRNAsData, mRs = mRs)
  
  return(l)
}

#================================================================
#' This function allows you to build the network for jointIda function
#' @param d Expression data.
#' @param pcmethod Method to build the network.
#' @param num.cores Number fo cores to be used.
#' @param mem.efficient Flag if using efficient memory.
#' @return Network used for joinIda
#================================================================
buildNetwork <- function(d, pcmethod = "parallel",
                         num.cores = 4, mem.efficient = FALSE) {
  data <- d
  data <- scale(data) # standardise the data
  suffStat <- list(C=cor(data), n = nrow(data))
  if ((pcmethod == "stable") || (pcmethod == "original")){
    pcFit <- pc_stable(suffStat, indepTest = gaussCItest, p = ncol(data), skel.method = pcmethod, alpha = 0.01)
  }else {
    pcFit <- pc_parallel(suffStat, indepTest = gaussCItest, p = ncol(data), skel.method = pcmethod, alpha = 0.01, num.cores = num.cores, mem.efficient = mem.efficient)
  }
  
  return(pcFit@graph)
}

#================================================================
#' This function allows you to find jointIda effect
#' @param dat Expression data.
#' @param cause Cause indexes.
#' @param effect Effect indexes.
#' @param method Method to get result.
#' @param num.cores Number of cores to be used.
#' @param network Network data.
#' @param tecnique Technique to be used.
#' @return joinIda effect
#================================================================
jointIdaGroup <- function(dat, cause, effect, method, num.cores, network, technique = "RRC") {
  data <- dat
  data <- scale(data) # Standardise the data
  result <- matrix(nrow = length(effect), ncol = length(cause))
  
  # Get number of cores to run
  cl<-makeCluster(num.cores)
  registerDoParallel(cl)
  
  # jointIda
  temp = rep(NA, length(cause))
  result.tmp <- foreach(k = 1:length(effect)) %dopar% {
    caef <- pcalg::jointIda(cause, effect[k], cov(data), network, technique=technique, type="cpdag")
    caefabs <- abs(caef)
    mat.tmp = temp
    for(l in 1:length(cause)){
      if(method=="min"||method=="max"){
        if(method=="min"){
          index <- which(caefabs==min(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }else{
          index <- which(caefabs==max(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }
        pos <- index[1, 2]
        mat.tmp[l] <- caef[l, pos]
      }else if(method=="median"){
        mat.tmp[l] <- median(caef[l, ], na.rm = TRUE)
      }
    }
    mat.tmp
  }
  
  # Shut down the workers
  stopCluster(cl)
  stopImplicitCluster()
  
  # Create a matrix from a list
  result = matrix(unlist(result.tmp), ncol = length(cause), byrow = T)
  
  # Calculate joint effect
  r <- NULL
  for (i in 1:length(effect)) {
    temp <- 0
    for (j in 1:length(cause)) {
      temp <- temp + result[i,j]*mean(data[,cause[j]])
    }
    r <- c(r, temp)
  }
  
  return(r)
}

#================================================================
#' This function allows you to find TFs regulating miRs
#' @param rootDir Folder to get TF-miRNA interaction database.
#' @param miRs miRNAs.
#' @return TFs regulating miRs
#================================================================
getTFsformiRs <- function(rootDir, miRs) {
  # TransmiR: TF => miRNA
  # TransmiR from http://www.cuilab.cn/transmir 
  interactions <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
                             sep = '\t', header = FALSE)
  interactions <- interactions[, 1:2]
  interactions <- unique(interactions)
  colnames(interactions) <- c("TF", "Target")
  interactions[,2] <- convertmiRs(interactions[,2])
  TF_Targets_TransmiR <- interactions[which(interactions$Target %in% miRs),]
  r <- as.character(unique(TF_Targets_TransmiR[,1]))
  return(r)
}

#================================================================
#' getGoldStandard
#================================================================
getGoldStandard <- function(rootDir, i_net, topk=5) {
  
  # bij is the noise-free steady state expression level of gene j in the ith double gene knockout experiment, where we let (c1(i), c2(i)) denote the pair of genes that is knocked out in the ith experiment.
  filename <- paste(rootDir, "/Data/DREAM4_InSilicoNetworks_GoldStandard/DREAM4_Challenge2_GoldStandards/Size 100/Size 100 bonus round/insilico_size100_", i_net, "_nonoise_dualknockouts.tsv", sep = "")
  b <- read.table(filename, header=FALSE, sep="\t", as.is = TRUE)
  
  # The noise-free steady state expression level
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net, "_wildtype.tsv", sep = "")
  w <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  
  # The gold-standard total increase in the expression value of gene j due to knocking out the gene pair (c1(i), c2(i))
  # delta <- b - w
  delta <- b
  delta <- as.matrix(delta)
  w <- as.matrix(w)
  for (j in 1:nrow(delta)) {
    delta[j,] <- delta[j,] - w[1,]
  }
  
  # Define our target set as the top k% of |deltaij| values for the 20 x 98 triples (c1(i), c2(i), j), i = 1, . . . , 20 and j belongs to {1, . . . , 100} \  {c1(i), c2(i)}.
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net, "_dualknockouts_indexes.tsv", sep = "")
  c <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  r <- matrix(nrow = 20, ncol = 2+topk)
  r[,1:2] <- as.matrix(c)
  delta <- abs(delta)
  for (j in 1:nrow(delta)) {
    delta[j, c(r[j,1], r[j,2])] <- 0
    r[j, 3:(2+topk)] <- head(sort(delta[j,], decreasing = TRUE, index.return=TRUE)$ix, topk)
  }
  
  return(r)
}

#================================================================
#' getEdges
#================================================================
getEdges <- function(net, cause) {
  cause <- paste("G", cause, sep = "")
  edgs <- net[which(net$V1 == cause),]
  if(nrow(edgs) == 0) {
    r <- NA
  } else {
    effects <- edgs$V2
    effects <- gsub("[G]", "", effects)
    effects <- as.double(effects)
    r <- list(edges=effects, weights=rep(1, nrow(edgs)))
  }
  
  return(r)
}

#================================================================
#' getNetwork
#================================================================
getNetwork <- function(rootDir, i_net) {
  
  # Network
  filename <- paste(rootDir, "/Data/DREAM4_InSilicoNetworks_GoldStandard/DREAM4_Challenge2_GoldStandards/Size 100/DREAM4_GoldStandard_InSilico_Size100_", i_net, ".tsv", sep = "")
  net <- read.table(filename, header=FALSE, sep="\t", as.is = TRUE)
  net <- net[which(net$V3 == 1),]
  
  # Build the network for jointIda
  p <- 100
  V <- as.character(1:p)
  my_edges <- vector("list", length = p)
  names(my_edges) <- c(1:p)
  for (i in 1:p) {
    my_edges[[as.character(i)]] <- getEdges(net, i)
  }
  my_edges <- lapply(my_edges, function(x) if(is.na(x[1])) NULL else x)  #changing NA to NULL 
  myDAG <- new("graphNEL", nodes=V, edgeL=my_edges, edgemode="directed") ## true DAG
  
  return(myDAG)
}

#================================================================
#' jointIdaForGroup
#================================================================
jointIdaForGroup <- function(rootDir, i_net, topk=5, technique = "RRC") {
  
  # Expression data
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net, "_knockouts.tsv", sep = "")
  net_data <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  net_data <- as.matrix(net_data)
  
  # The average expression level 
  net_data_aver <- colMeans(x=net_data, na.rm = TRUE)
  
  # Estimate PDAG
  p <- 100
  suffStat <- list(C = cor(net_data), n = nrow(net_data))
  pc.fit <- pc(suffStat, indepTest = gaussCItest, p = p, alpha = 0.01, u2pd="relaxed")
  pc.fit.pdag <- addBgKnowledge(pc.fit@graph)
  
  # jointIda
  filename <- paste(rootDir, "/Data/DREAM4_InSilico_Size100/insilico_size100_", i_net, "/insilico_size100_", i_net,
                    "_dualknockouts_indexes.tsv", sep = "")
  c <- read.table(filename, header=TRUE, sep="\t", as.is = TRUE)
  r_jointIda <- matrix(nrow = 20, ncol = p)
  r <- matrix(nrow = 20, ncol = 2+topk)
  r[, 1:2] <- as.matrix(c)
  for (i in 1:nrow(c)) {
    for (j in 1:p) {
      x_pos <- c(c[i,1], c[i,2])
      if(j %in% x_pos) {
        temp <- NULL
      } else {
        temp <- jointIda(x.pos=x_pos, y.pos=j, cov(net_data), graphEst=pc.fit.pdag, technique=technique, type = "pdag")    
      }
      if(!is.null(temp)) {
        temp <- -net_data_aver[c[i,1]]*(mean(x=temp[1,], na.rm = TRUE))-net_data_aver[c[i,2]]*(mean(x=temp[2,], na.rm = TRUE))
      } else {
        temp <- 0
      }
      r_jointIda[i, j] <- abs(temp)
    }
    r[i, 3:(2+topk)] <- head(sort(r_jointIda[i,], decreasing = TRUE, index.return=TRUE)$ix, topk)
  }
  
  return(r)
}

#================================================================
#' getRandomValidated
#================================================================
getRandomValidated <- function(goldRecord, noGene) {
  k <- 100
  res <- 0
  x <- goldRecord[3:100]
  x <- x[order(x, decreasing = FALSE)]
  x <- as.numeric(x)
  for (i in 1:k) {
    r <- sample(x, noGene)
    r_validated <- length(intersect(goldRecord[3:(2+noGene)], r))
    res <- res + r_validated
  }
  res <- res/100
  
  return(res)
}

#================================================================
#' Survival analysis (Survival curves, Log-rank test) and compute Silhouette information for cancer subtypes
#' 
#' Survival analysis is a very common tool to explain and validate the cancer subtype identification result. It provides the significance testing and 
#' graphical display for the verification of the survival patterns between the identified cancer subtypes.
#' 
#' @param mainTitle A character will display in the result plot.
#' @param time A numeric vector representing the survival time (days) of a set of samples.
#' @param status A numeric vector representing the survival status of a set of samples. 0=alive/censored, 1=dead.
#' @param group A vector represent the cluster label for a set of samples.
#' @param distanceMatrix A data matrix represents the similarity matrix or dissimilarity matrix between samples.\cr
#' If NULL, it will not compute silhouette width and draw the plot.
#' @param similarity A logical value. If TRUE, the distanceMatrix is a similarity distance matrix between samples. Otherwise a dissimilarity distance matrix between samples
#'
#' @return
#' The log-rank test p-value
#' 
#' @references
#' Xu T, Le TD, Liu L, Su N, Wang R, Sun B, Colaprico A, Bontempi G, Li J (2017). "CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation, and visualization." Bioinformatics. https://doi.org/10.1093/bioinformatics/btx378.
#================================================================
survivalAnalysis <- function(mainTitle = "Survival Analysis", time, status, group, distanceMatrix = NULL, similarity = TRUE){
  
  clusterNum = length(unique(group))
  dataset = list(time, status, x = group)  
  surv = survfit(Surv(time, status) ~ x, dataset)
  if(clusterNum>1){
    sdf=NULL
    sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
    cat("                                                     \n")
    cat("*****************************************************\n")
    cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
    print(sdf)
    p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  }else{
    cat("There is only one cluster in the group")
    p_value=1
  }
  
  if(!is.null(distanceMatrix[1,1])){
    layout(matrix(c(1,2,3,3), 2, 2, byrow = FALSE), widths=c(2.2,2), heights=c(2,2))
  }
  
  myCol <- wes_palette("Zissou1")
  
  # Graph 1
  title=paste(mainTitle, " (Number of clusters: ", clusterNum, ")", sep="")
  plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
       main = title, font.main=2, cex.lab=1.1, cex.axis=1, cex.main=1.2, cex.sub=1)
  legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
         lty=1, lwd=3, cex=1, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
         seg.len = 0.3)
  digit=ceiling(-log10(p_value)+2)
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1)   
  
  # Graph 2 & 3
  if(!is.null(distanceMatrix[1,1])){
    if(class(distanceMatrix)=="Similarity"){
      si=silhouette_SimilarityMatrix(group,distanceMatrix)
    }else{
      si=silhouette(group,distanceMatrix)
    }
    
    attr(distanceMatrix,'class')=NULL
    
    ind=order(group,-si[, "sil_width"])
    
    num=length(unique(group))
    annotation=data.frame(group=as.factor(group))
    # Var1 = c(palette()[2:(num+1)])
    Var1 = myCol[2:(num+1)]
    names(Var1) = sort(unique(group))
    ann_colors =  list(group=Var1)
    
    consensusmap(distanceMatrix,Rowv=ind,Colv=ind,main = "Clustering display",
                 annCol = annotation,annColors=ann_colors,
                 labRow ="Sample", labCol = "Sample",scale="none")
    
    plot(si,col =myCol[2:(clusterNum+1)], border=NA, cex.lab=1.1, cex.axis=1, cex.main=1, cex.sub=1)
  }
  
  # par(mfrow=c(1,1))
  return(p_value)
}

#================================================================
#' This function allows you to get type of a group: Additional group (ADD) or Enhanced group (ENH)
#' @param group The group to identify type.
#' @param singles Single genes.
#' @return Type of group
#================================================================
getType <- function(group, singles) {
  
  # Get cause genes
  cause_genes <- strsplit(group[1,1], split=' ', fixed=TRUE)
  cause_genes <- cause_genes[[1]]
  
  # Get influenced genes
  influenced_genes <- group$Influenced_Genes
  influenced_genes <- strsplit(influenced_genes, split=' ', fixed=TRUE)
  influenced_genes <- influenced_genes[[1]]
  
  # Get genes influenced by individual genes
  l <- NULL
  for (i in 1:group$Group_Size) {
    t <- singles[which(singles[,1] == cause_genes[i]),3]
    if (!isEmpty(t)) {
      t <- strsplit(t, split=' ', fixed=TRUE)
      t <- t[[1]]
      l <- union(l, t)  
    }
  }
  
  # Evaluate
  if (setequal(influenced_genes, l)) {
    r <- "ADD"
  } else {
    r <- "ENH"
  }
  
  return(r)
}

#================================================================
#' This function allows you to get regulators of a gene.
#' @param gene The gene to identify regulators
#' @param type 1 for coding, 2 for non-coding
#' @param dbDir Directoy which contains databases
#' @return Regulators of gene
#================================================================
getRegulators <- function(gene, type=1, dbDir) {
  
  if(type==1) { # coding
    # Get PPI network
    edges <- read_excel(paste(dbDir, "/PPI.xls",
                              sep = ""), sheet = 1)
    interactions <- edges[, c(1, 3)]
    colnames(interactions) <- c("cause", "effect")
    
    # miRNA => TF & miRNA => mRNA: Filter with
    # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
    confirmedList <- read.csv(paste(dbDir, "/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                    sep = ""), header = FALSE)
    colnames(confirmedList) <- c("cause", "effect")
    predictedList <- read.csv(paste(dbDir, "/TargetScan_7.0.csv", sep = ""), header = TRUE)
    predictedList <- predictedList[, -2]
    colnames(predictedList) <- c("cause", "effect")
    confirmedList[, 1] <- convertmiRs(confirmedList[, 1])
    predictedList[, 1] <- convertmiRs(predictedList[, 1])
    targetbinding <- rbind(confirmedList, predictedList)
    targetbinding <- unique(targetbinding)
    targetbinding$effect <- as.character(targetbinding$effect)
    
    # Combine
    interactions <- rbind(interactions, targetbinding)
    colnames(interactions) <- c("cause", "effect")
  } else { # non-coding
    # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
    interactions <- read.table(file = paste(dbDir, "/hsa.tsv", sep = ""),
                            sep = '\t', header = FALSE)
    interactions <- interactions[, 1:2]
    colnames(interactions) <- c("cause", "effect")
    interactions[, 2] <- convertmiRs(interactions[, 2])
    interactions <- unique(interactions)
    interactions$cause <- as.character(interactions$cause)
  }
  
  interactions <- interactions[which(interactions$effect == gene),]
  nodes <- interactions$cause
  
  return(nodes)
}

#================================================================
#' This function allows you to get common regulators of a gene group.
#' @param geneGroup The gene group to identify common regulators
#' @param type 1 for coding, 2 for non-coding
#' @param dbDir Directoy which contains databases
#' @return Regulators of gene group
#================================================================
getCommonRegulators <- function(geneGroup, type=1, dbDir) {
  
  n <- length(geneGroup)
  r <- getRegulators(geneGroup[1], type, dbDir)
  for (i in 2:n) {
    t <- getRegulators(geneGroup[i], type, dbDir)
    r <- intersect(r,t)
  }
  
  return(r)
}

#================================================================
#' This function allows you to get children of a gene.
#' @param gene The gene to identify children
#' @param type 1 for coding, 2 for non-coding
#' @param dbDir Directoy which contains databases
#' @return Children of gene
#================================================================
getChildren <- function(gene, type=1, dbDir) {
  
  if(type==1) { # coding
    # Get PPI network
    edges <- read_excel(paste(dbDir, "/PPI.xls",
                              sep = ""), sheet = 1)
    interactions <- edges[, c(1, 3)]
    colnames(interactions) <- c("cause", "effect")
    
    # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
    targetbinding <- read.table(file = paste(dbDir, "/hsa.tsv", sep = ""),
                               sep = '\t', header = FALSE)
    targetbinding <- targetbinding[, 1:2]
    colnames(targetbinding) <- c("cause", "effect")
    targetbinding[, 2] <- convertmiRs(targetbinding[, 2])
    targetbinding <- unique(targetbinding)
    targetbinding$cause <- as.character(targetbinding$cause)
    
    # Combine
    interactions <- rbind(interactions, targetbinding)
    colnames(interactions) <- c("cause", "effect")
  } else { # non-coding
    # miRNA => TF & miRNA => mRNA: Filter with
    # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
    confirmedList <- read.csv(paste(dbDir, "/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                    sep = ""), header = FALSE)
    colnames(confirmedList) <- c("cause", "effect")
    predictedList <- read.csv(paste(dbDir, "/TargetScan_7.0.csv", sep = ""), header = TRUE)
    predictedList <- predictedList[, -2]
    colnames(predictedList) <- c("cause", "effect")
    confirmedList[, 1] <- convertmiRs(confirmedList[, 1])
    predictedList[, 1] <- convertmiRs(predictedList[, 1])
    targetbinding <- rbind(confirmedList, predictedList)
    targetbinding <- unique(targetbinding)
    targetbinding$effect <- as.character(targetbinding$effect)
    interactions <- targetbinding
  }
  
  interactions <- interactions[which(interactions$cause == gene),]
  nodes <- interactions$effect
  
  return(nodes)
}

#================================================================
#' This function allows you to get common children of a gene group.
#' @param geneGroup The gene group to identify common children
#' @param type 1 for coding, 2 for non-coding
#' @param dbDir Directoy which contains databases
#' @return Children of gene group
#================================================================
getCommonChildren <- function(geneGroup, type=1, dbDir) {
  
  n <- length(geneGroup)
  r <- getChildren(geneGroup[1], type, dbDir)
  for (i in 2:n) {
    t <- getChildren(geneGroup[i], type, dbDir)
    r <- intersect(r,t)
  }
  
  return(r)
}

#================================================================
#' This function allows you to select the most variably expressed genes.
#' @param exp Gene expression: rows - genes, cols - samples
#' @param nsel Number of selected genes
#' @return Expression of selected genes
#================================================================
selectVariablyExpressedGenes <- function(exp, nsel = 500) {
  # subset this dataset to the nsel most variably expressed genes
  cvar <- apply(as.array(as.matrix(exp)), 1, sd)
  dat <- cbind(cvar, exp)
  dat <- dat[order(dat[,1], decreasing=T),]
  dat <- dat[1:nsel, -1]
  dat <- as.matrix(dat)
  
  return(dat)
}

#================================================================
#' This function allows you to construct networks for patients.
#' @param data Gene expression: rows - genes, cols - samples
#' @param method Method to build networks
#' @param rootDir Root folder
#' @return Networks
#================================================================
constructNetwork <- function(data, method, rootDir) {
  
  cancer_network <- NULL
  
  if (method == "LIONESS") {
    # model the single-sample networks based on co-expression using lionessR
    cormat <- lioness(data, netFun)
    row.names(cormat) <- paste(cormat[,1], cormat[,2], sep="_")
    cancer_network <- refineNetwork(cormat, rootDir)
    cancer_network <- as.data.table(cancer_network)
  }
  
  return(cancer_network)
}

#================================================================
#' This function allows you to construct networks for patients, using the new lioness function (i.e. return SummarizedExperiment).
#' @param data Gene expression: rows - genes, cols - samples
#' @param method Method to build networks
#' @param rootDir Root folder
#' @return Networks
#================================================================
constructNetwork2 <- function(data, method, rootDir) {
  
  cancer_network <- NULL
  
  if (method == "LIONESS") {
    # model the single-sample networks based on co-expression using lionessR
    cormat <- lioness(data, netFun)
    cormat <- assays(cormat)$lioness
    cormat <- cbind(reg = "", tar = "", cormat)
    n <- nrow(cormat)
    for (i in 1:n) {
      cormat[i,1] <- strsplit(row.names(cormat)[i],split='_', fixed=TRUE)[[1]][1]
      cormat[i,2] <- strsplit(row.names(cormat)[i],split='_', fixed=TRUE)[[1]][2]
    }
    cancer_network <- refineNetwork(cormat, rootDir)
    cancer_network <- as.data.table(cancer_network)
  }
  
  return(cancer_network)
}

#================================================================
#' This function allows you to construct networks for patients.
#' @param data Gene expression: rows - genes, cols - samples
#' @param method Method to build networks
#' @param rootDir Root folder
#' @return Networks
#================================================================
constructNetwork8 <- function(data, method, rootDir) {
  
  cancer_network <- NULL
  
  if (method == "LIONESS") {
    # model the single-sample networks based on co-expression using lionessR
    cormat <- lioness8(data, netFun8)
    row.names(cormat) <- paste(cormat[,1], cormat[,2], sep="_")
    cancer_network <- refineNetwork(cormat, rootDir)
    cancer_network <- as.data.table(cancer_network)
  }
  
  return(cancer_network)
}

#' LIONESS
#'
#' This function uses the LIONESS equation to estimate single-sample networks.
#' @param x Numeric matrix with samples in columns.
#' @param f Network reconstruction function. Defaults to Pearson correlation.
#' @keywords lioness
#' @export
#' @examples
#' exp <- matrix(sample(1000,1000)/1000, 100, 10)
#' colnames(exp) <- paste("sample", c(1:ncol(exp)), sep="_")
#' row.names(exp) <- paste("gene", c(1:nrow(exp)), sep="_")
#' lionessResults <- lioness(exp, netFun)

lioness8 <- function(x, f=netFun, ...){
  
  if(!is.function(f)){ stop("please use a function") }
  if(!is.matrix(x)) { print("please use a numeric matrix as input") }
  
  nrsamples <- ncol(x)
  samples <- colnames(x)
  
  # this applies netFun and extracts the aggregate network
  net <- f(x)
  agg <- c(net)
  
  # prepare the lioness output
  lionessOutput <- matrix(NA, nrow(net)*ncol(net), nrsamples+2)
  colnames(lionessOutput) <- c("reg", "tar", samples)
  lionessOutput[,1] <- rep(row.names(net), ncol(net))
  lionessOutput[,2] <- rep(colnames(net), each=nrow(net))
  lionessOutput <- as.data.frame(lionessOutput, stringsAsFactors=F)
  lionessOutput[,3:ncol(lionessOutput)] <- sapply(lionessOutput[,3:ncol(lionessOutput)], as.numeric)
  
  # run function f and the LIONESS equation
  for(i in 1:nrsamples){
    ss <- c(f(x[,-i])) # apply netFun on all samples minus one
    lionessOutput[,i+2] <- nrsamples*(agg-ss)+ss # apply LIONESS equation
  }
  return(lionessOutput)  
}

#' netFun
#'
#' This is the network reconstruction function that will be used to build aggregate networks.
#' @param x Numeric matrix with samples in columns.
#' @keywords netFun
#' @export
#' @examples
#' netFun()


netFun8 <- function(x, ...) {
  stats::cor(t(x), method="pearson") 
}

#================================================================
#' Refine networks based on the existing datasets
#' @param data List of edges, rows - edges, cols - start, end, samples
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
refineNetwork = function(data, rootDir){
  
  # nodes
  nodes <- unique(union(data[,1], data[,2]))
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with PPI interactions
  # Get PPI network
  edges <- read_excel(paste(rootDir, "/Data/PPI.xls",
                            sep = ""), sheet = 1)
  interactions <- edges[, c(1, 3)]
  colnames(interactions) <- c("cause", "effect")
  interactions <- interactions[which(interactions$cause %in% nodes),]
  interactions <- interactions[which(interactions$effect %in% nodes),]
  # miRNA => TF & miRNA => mRNA: Filter with
  # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
  confirmedList <- read.csv(paste(rootDir, "/Data/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                  sep = ""), header = FALSE)
  colnames(confirmedList) <- c("cause", "effect")
  predictedList <- read.csv(paste(rootDir, "/Data/TargetScan_7.0.csv", sep = ""), header = TRUE)
  predictedList <- predictedList[, -2]
  colnames(predictedList) <- c("cause", "effect")
  confirmedList[, 1] <- miRNAVersionConvert(confirmedList[, 1], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
  predictedList[, 1] <- miRNAVersionConvert(predictedList[, 1], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
  targetbinding <- rbind(confirmedList, predictedList)
  targetbinding <- unique(targetbinding)
  targetbinding$effect <- as.character(targetbinding$effect)
  # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
  interacts <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
                          sep = '\t', header = FALSE)
  interacts <- interacts[, 1:2]
  colnames(interacts) <- c("cause", "effect")
  interacts[, 2] <- miRNAVersionConvert(interacts[, 2], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
  interacts <- unique(interacts)
  interacts$cause <- as.character(interacts$cause)
  interacts <- interacts[which(!(is.na(interacts$effect))),]
  interacts$effect <- gsub("hsa-mir", "hsa-miR", interacts$effect)
  
  # Combine
  interactions <- rbind(interactions, targetbinding)
  interactions <- rbind(interactions, interacts)
  interactions <- unique(interactions)
  interactions <- interactions[which(!(is.na(interactions$cause))),]
  interactions <- interactions[which(!(is.na(interactions$effect))),]
  interactions <- as.data.frame(interactions)
  row.names(interactions) <- paste(interactions[,1], interactions[,2], sep="_")
  
  # refine
  r <- data[which(row.names(data) %in% row.names(interactions)),]
  r <- as.matrix(r)
  
  return(r)
}

#' #================================================================
#' #' Refine networks based on the existing datasets
#' #' @param data List of edges, rows - edges, cols - samples
#' #' @param rootDir Root folder
#' #' @return Edges of the network from cause to effect
#' #================================================================
#' refineNetwork2 = function(data, rootDir){
#'   
#'   # nodes
#'   reg_tar <- row.names(data)
#'   nodes <- strsplit(reg_tar, split = "_")
#'   nodes <- unlist(nodes)
#'   nodes <- unique(nodes)
#'   
#'   # Filter with existing datasets
#'   # TF => mRNA, mRNA => mRNA, TF => TF: Filter with PPI interactions
#'   # Get PPI network
#'   edges <- read_excel(paste(rootDir, "/Data/PPI.xls",
#'                             sep = ""), sheet = 1)
#'   interactions <- edges[, c(1, 3)]
#'   colnames(interactions) <- c("cause", "effect")
#'   interactions <- interactions[which(interactions$cause %in% nodes),]
#'   interactions <- interactions[which(interactions$effect %in% nodes),]
#'   # miRNA => TF & miRNA => mRNA: Filter with
#'   # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
#'   confirmedList <- read.csv(paste(rootDir, "/Data/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
#'                                   sep = ""), header = FALSE)
#'   colnames(confirmedList) <- c("cause", "effect")
#'   predictedList <- read.csv(paste(rootDir, "/Data/TargetScan_7.0.csv", sep = ""), header = TRUE)
#'   predictedList <- predictedList[, -2]
#'   colnames(predictedList) <- c("cause", "effect")
#'   confirmedList[, 1] <- miRNAVersionConvert(confirmedList[, 1], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
#'   predictedList[, 1] <- miRNAVersionConvert(predictedList[, 1], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
#'   targetbinding <- rbind(confirmedList, predictedList)
#'   targetbinding <- unique(targetbinding)
#'   targetbinding$effect <- as.character(targetbinding$effect)
#'   # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
#'   interacts <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
#'                           sep = '\t', header = FALSE)
#'   interacts <- interacts[, 1:2]
#'   colnames(interacts) <- c("cause", "effect")
#'   interacts[, 2] <- miRNAVersionConvert(interacts[, 2], targetVersion = "v21", exact = TRUE, verbose = FALSE)[,2]
#'   interacts <- unique(interacts)
#'   interacts$cause <- as.character(interacts$cause)
#'   interacts <- interacts[which(!(is.na(interacts$effect))),]
#'   interacts$effect <- gsub("hsa-mir", "hsa-miR", interacts$effect)
#'   
#'   # Combine
#'   interactions <- rbind(interactions, targetbinding)
#'   interactions <- rbind(interactions, interacts)
#'   interactions <- unique(interactions)
#'   interactions <- interactions[which(!(is.na(interactions$cause))),]
#'   interactions <- interactions[which(!(is.na(interactions$effect))),]
#'   interactions <- as.data.frame(interactions)
#'   row.names(interactions) <- paste(interactions[,1], interactions[,2], sep="_")
#'   
#'   # refine
#'   r <- data[which(row.names(data) %in% row.names(interactions)),]
#'   r <- as.matrix(r)
#'   
#'   return(r)
#' }

#================================================================
#' Identify cancer drivers for individuals by network control
#' @param network Network: rows - edges, cols - start, end, sample
#' @param nodes Nodes of networks
#' @param outDir Output directory
#' @param controlDir Directory to run network control
#' @param env Windows or Linux
#' @return Predicted drivers
#================================================================
identifyIndDriversbyControl <- function(network, nodes, outDir, controlDir, env) {
  
  drivers <- nodes
  
  # Analyse controllability of the network
  interactions <- network
  threshold <- mean(abs(network[,3]))
  interactions <- interactions[which(abs(interactions[,3]) >= threshold),1:2]
  # Write the edges of the network for analysing controllability
  mainDir <- paste(outDir, "/Controllability", sep = "")
  subDir <- gsub(".", "-", colnames(network)[3], fixed = TRUE)
  ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  write.table(interactions[,1:2], paste(mainDir, "/", subDir, "/edges.dat", sep = ""),
              row.names = FALSE, col.names=FALSE, quote=FALSE)
  # Run the controllability analysis
  if (env == "windows") {
    cmd <- paste(controlDir, "/parse.exe ", mainDir, "/", subDir, 
                 "/edges.dat", sep = "")
    system(cmd)
    cmd <- paste(controlDir, "/controllability_analysis.exe ", mainDir, "/", subDir,
                 "/edges.dat", sep = "")
    system(cmd)
  } else {
    cmd <- paste(controlDir, "/Parse ", mainDir, "/", subDir, 
                 "/edges.dat", sep = "")
    system(cmd)
    cmd <- paste(controlDir, "/ControllabilityAnalysis ", mainDir, "/", subDir,
                 "/edges.dat", sep = "")
    system(cmd)
  }
  
  # Identify critical nodes in the network
  # Read the result
  nodetype <- read.table(paste(mainDir, "/", subDir, "/edges.dat.nodetype", sep = ""))
  colnames(nodetype) <- c("Name", "K", "Kin", "Kout", "TypeI", "TypeII")
  # Critical nodes of the network
  critical_nodes <- nodetype[which(nodetype$TypeI == 0),]
  
  # output
  drivers[,2] <- ifelse(drivers[,1] %in% critical_nodes$Name, 1, 0)
  colnames(drivers) <- c("Gene", subDir)
  
  return(drivers)
  
}

#' #================================================================
#' #' Identify cancer drivers for individuals by Maximisation Influence
#' #' @param network Network: rows - edges, cols - start, end, sample
#' #' @param nodes Nodes of networks
#' #' @param outDir Output directory
#' #' @return Predicted drivers
#' #================================================================
#' identifyIndDriversbyMI <- function(network, nodes, outDir) {
#' 
#'   drivers <- nodes
#' 
#'   # Influence of active sets
#'   activeSets <- unique(union(network$reg, network$tar))
#'   n_activeSets <- length(activeSets)
#'   inf <- matrix(nrow = n_activeSets, ncol = 2)
#'   colnames(inf) <- c("Gene", "Influence")
#' 
#'   # Normalise node weight so that it is in [0,1]
#'   maxNodeWeight <- max(abs(nodes[,2]))
#'   nodes[,2] <- nodes[,2]/maxNodeWeight
#' 
#'   # Normalise edge weight
#'   colnames(network)[1:2] <- c("cause", "effect")
#'   network <- normaliseEdgeWeight(network)
#' 
#'   # Evaluate influence
#'   registerDoParallel(4)  # Use multi cores, set to the number of cores
#'   geneData <- activeSets
#'   geneData <- matrix(geneData, ncol = 1)
#'   r <- foreach (i=1:n_activeSets, .combine = c) %dopar% {
#'     evaluateInfluence(activeSets[i], geneData, network, nodes)
#'   }
#'   for (i in 1:n_activeSets) {
#'     inf[i,1] <- activeSets[i]
#'     inf[i,2] <- r[i]
#'   }
#'   stopImplicitCluster() # Clean up the cluster
#' 
#'   # output
#'   subDir <- gsub(".", "-", colnames(network)[3], fixed = TRUE)
#'   l <- nrow(drivers)
#'   for (i in 1:l) {
#'     ind <- which(inf[,1] == drivers[i,1])
#'     if(!is_empty(ind)) {
#'       v <- inf[ind,2]
#'     } else {
#'       v <- 0
#'     }
#'     drivers[i,2] <- v
#'   }
#'   colnames(drivers) <- c("Gene", subDir)
#' 
#'   return(drivers)
#' 
#' }
#' 
#' 

#================================================================
#' Get top drivers
#' @param d Influence of drivers: rows - genes, cols - gene, samples
#' @param n Number of selected tops
#' @return Top selected drivers
#================================================================
selectTopDrivers <- function(d, n = 200) {
  r <- d
  row.names(r) <- r[,1]
  
  nCol <- ncol(r)
  for (i in 2:nCol) {
    t <- r[, c(1,i)]
    t <- t[order(t[,2], decreasing = TRUE),]
    t <- t[1:n,1]
    r[,i] <- 0
    r[t,i] <- 1
  }
  
  return(r)
}

#================================================================
#' Convert edge list to network matrix
#' @param net Network: rows - edges, cols - start, end, sample
#' @param nodes Nodes of networks
#' @param hasWeight TRUE if the network has edge weights
#' @return Network matrix
#================================================================
convertToNetworkMatrix <- function(net, nodes, hasWeight = TRUE) {
  # Convert
  n <- nrow(nodes)
  m <- matrix(0, nrow = n, ncol = n)
  row.names(m) <- nodes[,1]
  colnames(m) <- nodes[,1]
  l <- nrow(net)
  for (i in 1:l) {
    m[net[i,1],net[i,2]] <- net[i,3]
  }
  
  # Remove weights
  if(!hasWeight) {
    m[which(m != 0)] <- 1
  }
  
  return(m)
}

#================================================================
#' Identify cancer drivers for individuals by PageRank
#' @param net Network: rows - edges, cols - start, end, sample
#' @param nodes Nodes of networks
#' @return Predicted drivers
#================================================================
identifyIndDriversbyPageRank <- function(net, nodes) {
  
  drivers <- nodes
  
  # network: row - regulators, col - targets
  threshold <- mean(abs(net[,3]))
  net <- net[which(abs(net[,3]) >= threshold),]
  # network <- convertToNetworkMatrix(net, nodes, hasWeight = FALSE)
  
  # set.seed(100)
  # g <- graph_from_data_frame(net, directed=TRUE, vertices=nodes)
  # ranking <- page_rank(g, algo="arpack")$vector
  
  set.seed(1000)
  g <- graph_from_data_frame(net, directed=TRUE, vertices=nodes)
  ranking <- page_rank(g, algo="arpack")$vector
  
  # output
  subDir <- gsub(".", "-", colnames(net)[3], fixed = TRUE)
  drivers[,2] <- ranking
  colnames(drivers) <- c("Gene", subDir)
  
  return(drivers)
  
}

#================================================================
#' Identify cancer drivers
#' @param network Network: rows - edges, cols - start, end, samples
#' @param nodes Nodes of networks
#' @param method Method to be used
#' @param outDir Output directory
#' @param controlDir Directory to run network control
#' @param env Windows or Linux
#' @return Predicted drivers
#================================================================
identifyDrivers <- function(network, nodes, method="control", outDir, controlDir, env = "Windows") {
  
  drivers <- nodes[,1]
  drivers <- as.matrix(drivers, ncol=1)
  colnames(drivers) <- "Gene"
  
  if(method == "control") {
    for (i in 3:ncol(network)) {
      temp <- identifyIndDriversbyControl(network[,c(1:2, i)], nodes, outDir, controlDir, env)
      drivers <- cbind(drivers, temp[,2])
      colnames(drivers)[i-1] <- colnames(temp)[2]
    }
  }
  
  # if(method == "MI") {
  #   for (i in 3:ncol(network)) {
  #     temp <- identifyIndDriversbyMI(network[,c(1:2, i)], nodes, outDir)
  #     drivers <- cbind(drivers, temp[,2])
  #     colnames(drivers)[i-1] <- colnames(temp)[2]
  #   }
  # }
  
  if(method == "PageRank") {
    for (i in 3:ncol(network)) {
      temp <- identifyIndDriversbyPageRank(network[,c(1:2, i)], nodes)
      drivers <- cbind(drivers, temp[,2])
      colnames(drivers)[i-1] <- colnames(temp)[2]
    }
  }
  
  return(drivers)
}

#================================================================
#' Remove genes whose expression is == 0 in more than 50% of the samples
#' @param x Gene expression: rows - genes, cols - samples or cells
#' @param filterParam Percentage to remove
#' @return Indexes to be removed
#================================================================
rem <- function(x, filterParam = 0.5){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*filterParam)
  return(remove)
}

#================================================================
#' Create table data for a box plot
#' @param cancer_type Cancer type, e.g. "BRCA"
#' @param drivers Matrix for personalised cancer drivers, rows - genes, columns - patients
#' @param gold_standard Gold standard (e.g. CGC)
#' @param n_genes Number of genes of interest
#' @return Table with 3 columns: "p_value", "Cancer_type", and "Patient"
#================================================================
createBoxPlotTable <- function(cancer_type, drivers, gold_standard, n_genes) {
  
  # Process data
  row.names(drivers) <- drivers[,1]
  drivers <- drivers[,-1]
  
  # Compute p-value for each patient
  n <- ncol(drivers) # Number of patients
  p_values <- matrix(NA, nrow = n, ncol = 3)
  colnames(p_values) <- c("p_value", "Cancer_type", "Patient")
  p_values[,2] <- cancer_type
  p_values[,3] <- colnames(drivers)
  # n_A: number of estimated cancer drivers
  # n_B: number of confirmed cancer drivers (get from gold standard)
  # n_A_B: number of cancer drivers validated
  # n_C: number of genes (get from TCGA)
  n_B <- length(gold_standard)
  n_C <- n_genes
  for (i in 1:n) {
    p_drivers <- rownames(drivers)[which(drivers[,i] == 1)]
    n_A <- length(p_drivers)
    n_A_B <- length(intersect(p_drivers, gold_standard))
    p_values[i,1] <-  1 - phyper(n_A_B, n_B, n_C-n_B, n_A)
  }
  tmp <- as.data.frame(p_values)
  tmp$p_value <- as.numeric(as.character(tmp$p_value))
  
  return(tmp)
}

#================================================================
#' Find proportion of patients driven by rare drivers in CGC
#' @param drivers Personalised drivers
#' @param rare_drivers Rare drivers in CGC
#' @return The proportion
#================================================================
computeRareProportion <- function(drivers, rare_drivers) {
  n_samples <- ncol(drivers)
  n_rare <- 0
  for (i in 1:n_samples) {
    p_drivers <- row.names(drivers)[which(drivers[,i] == 1)]
    if (all(rare_drivers %in% p_drivers)) {
      n_rare <- n_rare + 1
    }
  }
  return(n_rare/n_samples)
}

#================================================================
#' Find proportion of patients driven by rare drivers in CGC (4 rare lists)
#' @param drivers Personalised drivers
#' @param patients Patients
#' @param rare_drivers_1 Rare list 1
#' @param rare_drivers_2 Rare list 2
#' @param rare_drivers_3 Rare list 3
#' @param rare_drivers_4 Rare list 4
#' @param gold_standard Gold standard
#' @return The proportions
#================================================================
computeAllRareProportions <- function(drivers, patients, rare_drivers_1, rare_drivers_2, rare_drivers_3, rare_drivers_4, gold_standard) {
  
  drivers <- drivers[, which(colnames(drivers) %in% patients[,1])]
  # rare 1
  rare_1 <- computeRareProportion(drivers, intersect(rare_drivers_1$Gene, gold_standard))
  # rare 2
  rare_2 <- computeRareProportion(drivers, intersect(rare_drivers_2$Gene, gold_standard))
  # rare 3
  rare_3 <- computeRareProportion(drivers, intersect(rare_drivers_3$Gene, gold_standard))
  # rare 4
  rare_4 <- computeRareProportion(drivers, intersect(rare_drivers_4$Gene, gold_standard))
  
  return(c(rare_1, rare_2, rare_3, rare_4))
}
