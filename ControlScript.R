# Overview script
require(doMC)
require("yaml")
config = read_yaml("Execution_config.yml")
registerDoMC(cores = config$coreCount)

# Set the date which should be assigned to all output files
todays_date = format(Sys.Date(), "%d-%m-%Y")

# Obtain the relevant genes from the Next Generation Sequencing data
# TODO: "Unchanged" expressie status toevoegen
source("getNGSdata.R")

# Get the gene-gene relationships and their annotations
source("getGeneTriples.R")

# Create the features that are based on predicates
source("createPredicateFeatures.R")

# Calculate the features that are based on network statistics
source("createNetworkStatisticalFeatures.R") # Duurt lang, kijken of ik daar wat aan kan doen

# Create the reference set from the supplemental materials of Farashi et al.
source("getDbSNPdata.R")

# Combine the SNP data with the NGS data to identify the negative cases
source("CombineSNPwithDiffGenes.R")
combined = combined[combined$TargetGeneOnChromosome == T, ]

# Create the final feature set
source("createFeatureSet.R")