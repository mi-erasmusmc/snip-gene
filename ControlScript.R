# Overview script
#TODO: "Unchanged" expressie status toevoegen
require("yaml")
config = read_yaml("Execution_config.yml")

# Obtain the relevant genes from the Next Generation Sequencing data
source("getNGSdata.R")

# Get the gene-gene relationships and their annotations
source("getGeneTriples.R")

# Create the features that are based on predicates
source("createPredicateFeatures.R")

# Calculate the features that are based on network statistics
source("createNetworkStatisticalFeatures.R") # Duurt lang, kijken of ik daar wat aan kan doen

# Create the reference set from the supplemental materials of Farashi et al.
source("getDBSNPdata.R") #TODO: Niet elke SNP heeft een locatie

# Combine the SNP data with the NGS data to identify the negative cases
source("CombineSNPwithDiffGenes.R")
combined = combined[combined$TargetGeneOnChromosome == T, ]

# Create the final feature set
source("createFeatureSet.R")