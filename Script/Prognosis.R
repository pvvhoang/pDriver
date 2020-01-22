#=====================================
# This function is to format a list of miRs and mRNAs
# Input parameters:
#   miRmRNAList: a list with four columns - Source (miR), Target (mRNA), Effect, Group
# The output data is a list of below columns:
#   miRs
#   mRNA
#   Effect
#=====================================
formatList<-function(miRmRNAList){
  miRmRNAList <- miRmRNAList[,-4]

  # Format the list
  noOfmRNAs <- length(unique(as.character(miRmRNAList[,2])))
  l = matrix(nrow=noOfmRNAs, ncol=3)
  k <- 0
  startPoint <- 1
  while(startPoint < (nrow(miRmRNAList)+1)) {
    curmR <- as.character(miRmRNAList[startPoint,2])
    endPoint <- startPoint
    while((endPoint+1) < (nrow(miRmRNAList)+1)) {
      if (curmR == as.character(miRmRNAList[endPoint+1,2])) {
        endPoint <- endPoint + 1
      } else {
        break
      }
    }

    temp = as.character(miRmRNAList[startPoint,1])
    if(startPoint != endPoint) {
      for(i in (startPoint+1):endPoint) {
        temp <- paste(temp, as.character(miRmRNAList[i,1]), sep = "/")
      }
    }
    k <- k+1
    l[k,] <- c(temp, as.character(miRmRNAList[startPoint,2]), miRmRNAList[startPoint,3])

    # Move next
    startPoint <- endPoint + 1
  }
  row.names(l) <- l[,2]
  colnames(l) <- c("miRs", "mRNA", "Effect")

  return(l)
}

#=====================================
# This function is to merge 2 lists based on row.names
# Input parameters:
#   x: a list which has mRNAs in row.names
#   y: a list which has mRNAs in row.names
# The output data is a merged list
#=====================================
mergeLists<-function(x, y){
  result <- merge(x, l, by = "row.names", all.x = TRUE)
  colnames(result)[1] <- "mRNA"
  for (k in 1:ncol(result)){
    result[[k]] <- as.character(result[[k]])
  }
  result[is.na(result)] <- ""

  return(result)
}

#=====================================
# This function is to get related miRs
# Input parameters:
#   bindingListFinal: binding list between miRs and mRNAs
#   mRNAList: a list of mRNAs
#   delta: remove miRNAs that regulate less than delta of number of genes in a signature
# The output data is a list of miRs & mRNAs
#=====================================
getmiRsmRs<-function(bindingListFinal, mRNAList, delta){
  temp1 <- bindingListFinal
  colnames(temp1) <- c("miR", "mR")
  temp2 <- data.frame(mR=mRNAList)
  result <- merge(temp1, temp2, by.x = "mR", by.y = "mR", all.y = TRUE)
  result <- calEffectLevel(result)
  # result <- result[complete.cases(result),]
  result <- result[which(result[,3] > delta),]
  miRList <- result[,2]
  miRList <- unique(miRList)
  miRList <- as.character(miRList)
  mRNAList <- result[,1]
  mRNAList <- unique(mRNAList)
  mRNAList <- as.character(mRNAList)
  result = list(left=miRList, right=mRNAList)
  return(result)
}

#=====================================
# This function is to calculate the effect level of each miR
# Input parameters:
#   mRmiRList: a list of mRs and miRs
# The output data is a list of miRs & mRNAs & effect levels
#=====================================
calEffectLevel<-function(mRmiRList){
  result <- mRmiRList
  # Add one column
  newCol <- c("Effect Level")
  result[ , newCol] <- 0
  # Get number of mRNAs
  mRNAList <- result[,1]
  mRNAList <- unique(mRNAList)
  noOfmRNAs <- length(mRNAList)
  # Calculate effect levels
  for(i in 1:nrow(result)) {
    if(!is.na(result[i,2])) {
      curmiR <- as.character(result[i,2])
      count <- length(which(result[,2] == curmiR))
      result[i,3] <- count/noOfmRNAs
    }
  }
  # Return result
  return(result)
}

#=====================================
# This function is to predict the risk and calculate concordance index
# Input parameters:
#   trainingData: a training data set
#   testData: a test data set
#   survData: the survival data of testData
#   type: algorithm
# The output data is a list including an object of class coxph representing the fit, prediction and concordance index
#=====================================
predictRisk<-function(trainingData, testData, survData, type = "coxph"){
  if(type == "coxph") {
    ## Identify the fit
    fit <- coxph(Surv(time, event) ~ ., data = data.frame(trainingData))
    ## Predict
    pred <- predict(fit, newdata = data.frame(testData), type = "risk")
    ## Calculate concordance index
    perf <- concordance.index(x = pred, surv.time = survData[, "time"], surv.event = survData[, "event"], method = "noether", na.rm = TRUE)
    cIndex <- perf[1]
  # } else if(type == "KM") {
  #   fit <- survfit(Surv(time, event) ~ 1, data=data.frame(trainingData))
  # } else if(type == "aareg") {
  #   fit <- aareg(Surv(time, event) ~ . , data = data.frame(trainingData))
  # } else if(type == "ranger") {
  #   fit <- ranger(Surv(time, event) ~ ., data = data.frame(trainingData),
  #                   importance = "permutation",
  #                   splitrule = "extratrees",
  #                   verbose = TRUE)
  } else if(type == "rf") {
    temp <- data.frame(trainingData)[data.frame(trainingData)[,1] > 0,]
    fit <- rfsrc(Surv(time, event) ~ ., temp, ntree = 100)
    pred <- predict(fit, data.frame(testData))
    perf <- concordance.index(x = pred$predicted, surv.time = survData[, "time"], surv.event = survData[, "event"], method = "noether", na.rm = TRUE)
    cIndex <- perf[1]
  }

  # Return result
  result = list(fit = fit, prediction = pred, cIndex = cIndex)
  return(result)
}

#=====================================
# This function is to convert miR names to matured ones
# Input parameters:
#   testData: a data set with columns being miRs
#   startPoint: start column to convert
#   endPoint: end column to convert
# The output data is the data set with matured miR names
#=====================================
convertmiRNames <- function(testData, startPoint, endPoint){
  library(miRBaseConverter)

  miRNameList <- colnames(testData)[startPoint:endPoint]
  version=checkMiRNAVersion(miRNameList,verbose=FALSE)
  # Convert non mature miRNAs' names to mature names
  miRMature=miRNA_PrecursorToMature(miRNameList, version=version)
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")

  # Process miRNAs' names
  # Update the converted items to the miR list
  # if a converted item is "NA", get the corresponding value from the mature list
  miRNameList <- miRNameList_Converted[,2]
  naList <- which(is.na(miRNameList_Converted[,2]))
  for(i in 1:length(naList)) {
    k <- naList[i]
    miRNameList[k] <- miRMature[k,2]
  }
  colnames(testData)[startPoint:endPoint] <- miRNameList
  return(testData)
}

#=====================================
# This function is to process test data
# Input parameters:
#   testData: a data set with columns being miRs/ mRNAs
#   matchedResult: mirs/ mRNAs are needed to get, in case there are many records for a miR/ mRNA in test data, get the average value
# The output data is the test data
#=====================================
processTestData <- function(testData, matchedResult){
  result = matrix(nrow = nrow(testData), ncol = length(matchedResult))
  row.names(result) <- row.names(testData)
  colnames(result) <- matchedResult

  for(i in 1:length(matchedResult)) {
    t <- testData[,which(colnames(testData) %in% matchedResult[i])]
    if(!is.null(ncol(t)) && ncol(t) > 1) {
      for(j in 1:nrow(result)) {
        result[j,i] <- mean(t[j,])
      }
    } else {
      result[,i] <- t
    }
  }

  return(result)
}

#=====================================
# This function is to get unique miRs
# Input parameters:
#   miRs: a list of miRs
# The output data is the miR list
#=====================================
getmiRs <- function(miRs){
  result = NULL
  for(i in 1:length(miRs)) {
    t <- as.character(miRs[i])
    if(t != "") {
      result <- c(result, unlist(strsplit(t, "/")))
    }
  }
  result <- unique(result)
  return(result)
}

#=====================================
# This function is to predict risks using mRs
# Input parameters:
#   RNASeq_data: training data set
#   survival_data: survival data of training data set
#   dataMETABRIC_exprs: test data set
#   mRNAs: a list of mRs
#   type: algorithm
# The output data is a list including an object of class coxph representing the fit, prediction and concordance index
#=====================================
predictRiskUsingmR <- function(RNASeq_data, survival_data, dataMETABRIC_exprs, mRNAs, type = "coxph"){
  # Get survival info
  survData <- survival_data[,c("new_death", "death_event")]
  colnames(survData) <- c("time", "event")

  # Append gene expressions
  dd <- RNASeq_data[which(row.names(RNASeq_data) %in% mRNAs), ]
  if(is.null(dim(dd))) {# there is only one record
    dd <- cbind(survData, dd)
    colnames(dd)[ncol(dd)] <- mRNAs[1]
  } else {
    dd <- t(dd)
    dd <- cbind(survData, dd)
  }

  # Using METABRIC as test data
  testData <- t(exprs(dataMETABRIC_exprs)) # col 1-823: miRs, col 824-25191: mRNAs

  # Just get the intersection between training data set and test data set
  l = list()
  k <- 1
  l[[k]] <- colnames(dd)[3:(ncol(dd))]
  k <- k+1
  l[[k]] <- colnames(testData)
  result = Reduce(intersect, l)
  if(length(result) == 0) {
    print("There is no intersection between training data set and test data set!")
    return(NA)
  } else {
    print(result)
  }

  # update training data
  match = match(result, colnames(dd))
  dd <- dd[,c(1,2,match)]

  # update test data set
  testData <- processTestData(testData, result)

  # Get the survival data of test data
  survData <- phenoData(dataMETABRIC_exprs)@data [, c("t.rfs", "e.rfs")]
  colnames(survData) <- c("time", "event")

  # Run prediction
  result <- predictRisk(dd, testData, survData, type)

  return(result)
}

#=====================================
# This function is to predict risks using miRs
# Input parameters:
#   survival_data: survival data of training data set
#   dataMETABRIC_exprs: test data set
#   miRs: a list of miRs
#   directoryPath: directory
#   type: algorithm
# The output data is a list including an object of class coxph representing the fit, prediction and concordance index
#=====================================
predictRiskUsingmiR <- function(survival_data, dataMETABRIC_exprs, miRs, directoryPath, type = "coxph"){
  # Get survival info
  survData <- survival_data[,c("new_death", "death_event")]
  colnames(survData) <- c("time", "event")

  # Append gene expressions
  filename <- paste(directoryPath, "datacsv.csv", sep="")
  datacsv <- read.csv(filename, header=TRUE, sep=",")
  colnames(datacsv) <- gsub("[.]", "-", colnames(datacsv))
  dd <- datacsv[,which(colnames(datacsv) %in% miRs)]
  if(is.null(dim(dd))) {# there is only one record
    dd <- cbind(survData, dd)
    colnames(dd)[ncol(dd)] <- miRs[1]
  } else {
    dd <- cbind(survData, dd)
  }

  # Using METABRIC as test data
  testData <- t(exprs(dataMETABRIC_exprs)) # col 1-823: miRs, col 824-25191: mRNAs
  startPoint <- 1
  endPoint <- 823
  testData <- convertmiRNames(testData, startPoint, endPoint)
  miRs <- testData[, startPoint:endPoint]
  miRs <- sapply(split(seq_len(ncol(miRs)), colnames(miRs)), function(cis) rowMeans(miRs[, cis, drop=F]))
  testData <- testData[,-(startPoint:endPoint)]
  testData <- cbind(miRs, testData)

  # Just get the intersection between training data set and test data set
  l = list()
  k <- 1
  l[[k]] <- colnames(dd)[3:(ncol(dd))]
  k <- k+1
  l[[k]] <- colnames(testData)
  result = Reduce(intersect, l)
  if(length(result) == 0) {
    print("There is no intersection between training data set and test data set!")
    return(NA)
  } else {
    print(result)
  }

  # update training data
  match = match(result, colnames(dd))
  dd <- dd[,c(1,2,match)]

  # update test data set
  testData <- processTestData(testData, result)

  # Get the survival data of test data
  survData <- phenoData(dataMETABRIC_exprs)@data [, c("t.rfs", "e.rfs")]
  colnames(survData) <- c("time", "event")

  # Run prediction
  result <- predictRisk(dd, testData, survData, type)

  return(result)
}

#=====================================
# This function is to predict risks using both mRs and miRs
# Input parameters:
#   survival_data: survival data of training data set
#   dataMETABRIC_exprs: test data set
#   mRNAs: a list of mRs
#   miRs: a list of miRs
#   directoryPath: directory
#   type: algorithm
# The output data is a list including an object of class coxph representing the fit, prediction and concordance index
#=====================================
predictRiskUsingmRAndmiR <- function(survival_data, dataMETABRIC_exprs, mRNAs, miRs, directoryPath, type = "coxph"){
  # Get survival info
  survData <- survival_data[,c("new_death", "death_event")]
  colnames(survData) <- c("time", "event")

  # Append gene expressions
  filename <- paste(directoryPath, "datacsv.csv", sep="")
  datacsv <- read.csv(filename, header=TRUE, sep=",")
  colnames(datacsv) <- gsub("[.]", "-", colnames(datacsv))
  dd <- datacsv[,which(colnames(datacsv) %in% miRs)]
  if(is.null(dim(dd))) {# there is only one record
    dd <- cbind(survData, dd)
    colnames(dd)[ncol(dd)] <- miRs[1]
  } else {
    dd <- cbind(survData, dd)
  }
  dd2 <- datacsv[,which(colnames(datacsv) %in% mRNAs), ]
  if(is.null(dim(dd2))) {# there is only one record
    dd <- cbind(dd, dd2)
    colnames(dd)[ncol(dd)] <- mRNAs[1]
  } else {
    dd <- cbind(dd, dd2)
  }

  # Using METABRIC as test data
  testData <- t(exprs(dataMETABRIC_exprs)) # col 1-823: miRs, col 824-25191: mRNAs
  startPoint <- 1
  endPoint <- 823
  testData <- convertmiRNames(testData, startPoint, endPoint)
  miRs <- testData[, startPoint:endPoint]
  miRs <- sapply(split(seq_len(ncol(miRs)), colnames(miRs)), function(cis) rowMeans(miRs[, cis, drop=F]))
  testData <- testData[,-(startPoint:endPoint)]
  testData <- cbind(miRs, testData)

  # Just get the intersection between training data set and test data set
  l = list()
  k <- 1
  l[[k]] <- colnames(dd)[3:(ncol(dd))]
  k <- k+1
  l[[k]] <- colnames(testData)
  result = Reduce(intersect, l)
  if(length(result) == 0) {
    print("There is no intersection between training data set and test data set!")
    return(NA)
  } else {
    print(paste("Number of miRNAs & mRNAs: ", length(result), sep = ""))
  }

  # Update traing data set
  match = match(result, colnames(dd))
  dd <- dd[,c(1,2,match)]

  # Update test data set
  testData <- processTestData(testData, result)

  # Get the survival data of test data
  survData <- phenoData(dataMETABRIC_exprs)@data [, c("t.rfs", "e.rfs")]
  colnames(survData) <- c("time", "event")

  # Run prediction
  result <- predictRisk(dd, testData, survData, type)

  return(result)
}

# makeDataSets <- function(train.dataset = NULL, test.dataset = NULL, list = NULL){
#   temp=exprs(train.dataset)
#   rownames(temp)=miRs2
#   temp = apply(temp, 2, function(x) tapply(x, rownames(temp), mean))
#   idx = which(rownames(temp) %in% list)
#   if (length(idx) == 0) return (list("train.task" = NULL, "test.task" = NULL))
#   temp=temp[idx,]
#   expd = t(temp)
#   clind = pData(train.dataset)[,c('t.rfs',"e.rfs")]
#   t1 = as.numeric(clind[,1])
#   t1[t1<0] <- 0
#   t2 = as.numeric(clind[,2])
#   clind = data.frame(cbind(t1,t2))
#   colnames(clind) = c("time", "status")
#   data = cbind(expd, clind)
#   colnames(data) = gsub("-", "\\.", colnames(data))
#   #train.task = data
#   train.task <- makeSurvTask(data = data, target = c("time", "status"))
#
#   ####Build testing task
#   temp = t(exprs(test.dataset))
#   colnames(temp)=miRs1
#   temp = apply(temp, 1, function(x) tapply(x, colnames(temp), mean))
#   train.genes = colnames(expd)
#   expd = NULL
#   for (i in train.genes) {
#     idx = which(colnames(temp) == i)
#     if (length(idx)) expd = cbind(expd, temp[,idx])
#     else expd = cbind(expd, rep(0, ncol(temp)))
#   }
#   colnames(expd) = train.genes
#   clind = pData(test.dataset)[,c('t.rfs',"e.rfs")]
#   t1 = as.numeric(clind[,1])
#   t2 = as.numeric(clind[,2])
#   clind = data.frame(cbind(t1,t2))
#   colnames(clind) = c("time", "status")
#   data = cbind(expd, clind)
#   colnames(data) = gsub("-", "\\.", colnames(data))
#   #test.task = data
#   test.task <- makeSurvTask(data = data, target = c("time", "status"))
#
#   return(list("train.task" = train.task, "test.task" = test.task))
# }
