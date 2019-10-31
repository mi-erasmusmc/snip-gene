# Overview script
require(doMC)
require("yaml")
config = read_yaml("Execution_config.yml")
registerDoMC(cores = config$coreCount)

# Set the date which should be assigned to all output files
todays_date = format(Sys.Date(), "%d-%m-%Y")

# Identify which genes are expressed either in normal or in cancer prostate tissue
source("IdentifyExpressedGenes.R")

# Obtain the relevant genes from the Next Generation Sequencing data
source("getNGSdata.R")

# Get the gene-gene relationships and their annotations
source("getGeneTriples.R")

# Create the features that are based on predicates
source("createPredicateFeatures.R")

# Calculate the features that are based on network statistics
source("createNetworkStatisticalFeatures.R") # Duurt lang, kijken of ik daar wat aan kan doen

# Create the reference set from the supplemental materials of Farashi et al.
# Identify other genes in the neighborhood
source("getDbSNPdata.R")
farashi_final = farashi_final[farashi_final$TargetGeneOnChromosome == T, ]

# Create the final feature set
source("createFeatureSet.R")
source("createFeatureSetOfClosestGenes.R")
