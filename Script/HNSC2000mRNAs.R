#=========================================================================
#=========================================================================
# pDriver: A novel method for unravelling personalised coding and non-coding cancer drivers
#=========================================================================
#=========================================================================
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
library(lionessR)
library(readxl)
library(data.table)
library(miRBaseConverter)

# # Please install the R package lionessR via the devtools package from CRAN:
# install.packages("devtools")
# library(devtools)
# devtools::install_github("mararie/lionessR")

# if using lionessR from devtools, use constructNetwork
# if using lionessR from Bioconductor, use constructNetwork2

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
rootDir <- "/home/users/phamvv/pDriver" # And put the input files
# in "rootDir/Data"
outDir <- "/home/users/phamvv/pDriver/Data/HNSC" # Output folder
controlDir <- "/home/users/phamvv/Control" # Put here the library from control network
numCores <- 4 # Number of cores to be used
nomiR <- 100
nomR <- 2000
#nomR <- 200
#nomR <- 500
#nomR <- 1000
nomRforSC <- 400
outDir <- paste(outDir, "/", nomR , "mRs", sep = "")
env <- "Linux"
#---------------------------------------

# Include the script of functions
source(paste(rootDir, "/Script/PersonalisedDriverFunctions.R", sep=""))
source(paste(rootDir, "/Script/Prognosis.R", sep=""))

# Main script
#================================================================
# (1) Constructing the miRNA-TF-mRNA network for cancer patients
#================================================================

#---------------------------------------
# Common
#---------------------------------------

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("1. Start of (1) Constructing the miRNA-TF-mRNA network for cancer patients", sep="\n", file = log_con)
close(log_con)
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************

# load the dataset
# Load the tumor expression data
load(paste(rootDir, "/Data/HNSC_matchedData_full.RData", sep = "")) # rows - samples, cols - genes
BRCA_matchedData <- HNSC_matchedData

# Get PPI network
edges <- read_excel(paste(rootDir, "/Data/PPI.xls",
                          sep = ""), sheet = 1)
interactions <- edges[, c(1, 3)]
colnames(interactions) <- c("cause", "effect")
interactions <- interactions[which(interactions$cause %in% colnames(BRCA_matchedData$mRNAs)),]
interactions <- interactions[which(interactions$effect %in% colnames(BRCA_matchedData$mRNAs)),]
nodes <- unique(union(interactions$cause, interactions$effect))

# TFs: Download the list from http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19
tfs <- read.csv(paste(rootDir, "/Data/Browse Transcription Factors hg19 - resource_browser.csv",
                      sep = ""))
i <- which(levels(tfs$Symbol) %in% nodes)
tfData <- BRCA_matchedData$mRNAs[, levels(tfs$Symbol)[i]]

# Update cancer data of mRNAs
mRNAsData_Cancer <- BRCA_matchedData$mRNAs[,
                                           nodes[which(!(nodes %in% levels(tfs$Symbol)[i]))]]

# Get the cancer data of miRNAs
miRNAsData_Cancer <-  BRCA_matchedData$miRs

# Select genes
miRNAsData_Cancer <- t(selectVariablyExpressedGenes(t(miRNAsData_Cancer), nomiR))
mRNAsData_Cancer <- t(selectVariablyExpressedGenes(t(mRNAsData_Cancer),nomR))

# Covert miRNAs
colnames(miRNAsData_Cancer) <- miRNAVersionConvert(colnames(miRNAsData_Cancer),
                                                   targetVersion = "v21",
                                                   exact = TRUE,
                                                   verbose = FALSE)[,2]

# Combine data
noTF <- ncol(tfData)
cancer_data <- cbind(miRNAsData_Cancer, mRNAsData_Cancer, tfData)

# change to: rows - genes, cols - samples
cancer_data <- t(cancer_data)

# Free the memory
gc()

# Save data
saveRDS(cancer_data, file = paste(outDir, "/cancer_data.rds", sep = ""))

#---------------------------------------
# (1.1) LIONESS
#---------------------------------------

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("2. Start of (1.1) LIONESS - construct network", sep="\n", file = summary(log_con)$description, append = TRUE)
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************

# Read data
cancer_data <- readRDS(file = paste(outDir, "/cancer_data.rds", sep = ""))

# construct network
# cancer_network_LIONESS <- constructNetwork8(cancer_data, method="LIONESS", rootDir)

#@@@@@@@@@@@
data <- cancer_data
method="LIONESS"

cancer_network <- NULL

# model the single-sample networks based on co-expression using lionessR
# cormat <- lioness8(data, netFun8)

#@@@@@@@@@@@@@@@@
x <- data
f <- netFun8

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
#lionessOutput[,3:ncol(lionessOutput)] <- sapply(lionessOutput[,3:ncol(lionessOutput)], as.numeric)
nr_lionessOutput <- nrow(lionessOutput)
nr_half <- round(nr_lionessOutput/2)
nr_half_half <- round(nr_half/2)
lionessOutput[1:nr_half_half,3:ncol(lionessOutput)] <- sapply(lionessOutput[1:nr_half_half,3:ncol(lionessOutput)], as.numeric)
lionessOutput[(nr_half_half+1):nr_half,3:ncol(lionessOutput)] <- sapply(lionessOutput[(nr_half_half+1):nr_half,3:ncol(lionessOutput)], as.numeric)
lionessOutput[(nr_half+1):(nr_half+nr_half_half),3:ncol(lionessOutput)] <- sapply(lionessOutput[(nr_half+1):(nr_half+nr_half_half),3:ncol(lionessOutput)], as.numeric)
lionessOutput[(nr_half+nr_half_half+1):nr_lionessOutput,3:ncol(lionessOutput)] <- sapply(lionessOutput[(nr_half+nr_half_half+1):nr_lionessOutput,3:ncol(lionessOutput)], as.numeric)

# run function f and the LIONESS equation
for(i in 1:nrsamples){
  ss <- c(f(x[,-i])) # apply netFun on all samples minus one
  lionessOutput[,i+2] <- nrsamples*(agg-ss)+ss # apply LIONESS equation
}

cormat <- lionessOutput
#@@@@@@@@@@@@@@@@

row.names(cormat) <- paste(cormat[,1], cormat[,2], sep="_")
cancer_network <- refineNetwork(cormat, rootDir)
cancer_network <- as.data.table(cancer_network)

cancer_network_LIONESS <-cancer_network
#@@@@@@@@@@@

# Save the network
fwrite(cancer_network_LIONESS, paste(outDir, "/cancer_network_LIONESS.csv", sep = ""), row.names = FALSE)

#---------------------------------------
# Common (Just run 1 time)
#---------------------------------------

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("3. Common (Just run 1 time)", sep="\n", file = summary(log_con)$description, append = TRUE)
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************

# Compute node weight
load(paste(rootDir, "/Data/HNSC_matchedData_normal_samples_full.RData", sep = "")) # rows - samples, cols - genes
BRCA_matchedData_normal_samples <- HNSC_matchedData_normal_samples
miRNAsData_Normal <- BRCA_matchedData_normal_samples$miRs
# Covert miRNAs
colnames(miRNAsData_Normal) <- miRNAVersionConvert(colnames(miRNAsData_Normal),
                                                   targetVersion = "v21",
                                                   exact = TRUE,
                                                   verbose = FALSE)[,2]
miRNAsData_Normal <-  miRNAsData_Normal[, colnames(miRNAsData_Cancer)]
mRNAsData_Normal <- BRCA_matchedData_normal_samples$mRNAs[, colnames(mRNAsData_Cancer)]
tfData_Normal <- BRCA_matchedData_normal_samples$mRNAs[, colnames(tfData)]
normal_data <- cbind(miRNAsData_Normal, mRNAsData_Normal, tfData_Normal)
nodeList <- computeNodeWeight(t(cancer_data), normal_data)
# Save file
write.csv(nodeList, paste(outDir, "/nodeList.csv", sep = ""), row.names = TRUE)

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("4. End of Common (Just run 1 time)", sep="\n", file = summary(log_con)$description, append = TRUE)
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************

#================================================================
# (2) Identifying cancer drivers
#================================================================

#---------------------------------------
# (2.1) LIONESS - control
#---------------------------------------

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("5. Start of (2) Identifying cancer drivers - (2.1) LIONESS - control", sep="\n", file = summary(log_con)$description, append = TRUE)
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************

network <- read.csv(paste(outDir, "/cancer_network_LIONESS.csv", sep = ""), as.is = TRUE)
nodes <- read.csv(paste(outDir, "/nodeList.csv", sep = ""), as.is = TRUE)
drivers <- identifyDrivers(network, nodes, method="control", outDir, controlDir, env)
# Save drivers
drivers <- as.data.table(drivers)
fwrite(drivers, paste(outDir, "/drivers_LIONESS_control.csv", sep = ""), row.names = FALSE)

#**************************************************
# Write to log
log_con <- file(paste(outDir, "/", nomiR, "miRs", nomR, "mRs.log", sep = ""))
cat("6. End of (2) Identifying cancer drivers - (2.1) LIONESS - control", sep="\n", file = summary(log_con)$description, append = TRUE)
cat(paste0(Sys.time()), sep="\n", file = summary(log_con)$description, append = TRUE)
close(log_con)
#**************************************************