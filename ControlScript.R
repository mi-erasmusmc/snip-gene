# Overview script
#TODO: "Unchanged" expressie status toevoegen
require("yaml")
config = read_yaml("Execution_config.yml")

# Obtain the relevant genes from the Next Generation Sequencing data
calculateAnew = config$calculateAnew
source("getNGSdata.R")

# Get the gene-gene relationships and their annotations
getFromKG = config$getFromKG
source("getGeneTriples.R")

# Create the features that are based on predicates
generatePredicateFeatures = config$generatePredicateFeatures
source("createPredicateFeatures.R")

# Calculate the features that are based on network statistics
recalculateNetwork = config$recalculateNetwork
source("createNetworkStatisticalFeatures.R") # Duurt lang, kijken of ik daar wat aan kan doen

# Create the reference set from the supplemental materials of Farashi et al.
getSNPDataFromWeb = config$getSNPDataFromWeb
source("getDBSNPdata.R") #TODO: Niet elke SNP heeft een locatie
farashi_final = farashi_final[!is.na(farashi_final$BP) & !is.na(farashi_final$Target_Gene), ]
farashi_final = farashi_final[farashi_final$SNP.s.Genomic.Location != "Coding region", ]

# Combine the SNP data with the NGS data to identify the negative cases
diffGenesOnly = config$diffGenesOnly
source("CombineSNPwithDiffGenes.R")
combined = combined[combined$TargetGeneOnChromosome == T, ]

# Create the final feature set
onlySelectNegatives = config$onlySelectNegatives
source("createFeatureSet.R")